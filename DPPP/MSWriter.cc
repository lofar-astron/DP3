//# MSWriter.cc: DPPP step writing to an MS
//# Copyright (C) 2010
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include "MSWriter.h"
#include "MSUpdater.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Version.h"

#include "../Common/VdsMaker.h"
#include "../Common/ParameterSet.h"

#include <casacore/tables/Tables/TableCopy.h>
#include <casacore/tables/DataMan/DataManInfo.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledStManAccessor.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/OS/Path.h>

#include <iostream>
#include <limits>
#include <measures/TableMeasures/TableMeasDesc.h>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    MSWriter::MSWriter (MSReader* reader, const string& outName,
                        const ParameterSet& parset, const string& prefix)
      : itsReader       (reader),
        itsName         (prefix),
        itsOutName      (outName),
        itsParset       (parset),
        itsNrDone       (0)
    {
      // Get tile size (default 1024 KBytes).
      itsTileSize          = parset.getUint (prefix+"tilesize", 1024);
      itsTileNChan         = parset.getUint (prefix+"tilenchan", 0);
      itsOverwrite         = parset.getBool (prefix+"overwrite", false);
      itsNrTimesFlush      = parset.getUint (prefix+"flush", 60);
      itsCopyCorrData      = parset.getBool (prefix+"copycorrecteddata", false);
      itsCopyModelData     = parset.getBool (prefix+"copymodeldata", false);
      itsWriteFullResFlags = parset.getBool (prefix+"writefullresflag", true);
      itsDataColName       = parset.getString (prefix+"datacolumn", "DATA");
      itsWeightColName     = parset.getString (prefix+"weightcolumn",
                               "WEIGHT_SPECTRUM");
      itsVdsDir            = parset.getString (prefix+"vdsdir", string());
      itsClusterDesc       = parset.getString (prefix+"clusterdesc", string());
      if(itsDataColName != "DATA")
        throw Exception("Currently only the DATA column"
                 " can be used as output when writing a new MS");
      if(itsWeightColName != "WEIGHT_SPECTRUM")
        throw Exception("Currently only the "
          "WEIGHT_SPECTRUM column can be used as output when writing a new MS");
      
      itsStManKeys.Set(parset, prefix);
    }

    MSWriter::~MSWriter()
    {}

    bool MSWriter::process (const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      // Form the vector of the output table containing new rows.
      Vector<uint> rownrs(itsNrBl);
      indgen (rownrs, itsMS.nrow());
      // Add the necessary rows to the table.
      itsMS.addRow (itsNrBl);
      // Form the subset of the tables containing the rows.
      // It can happen that a missing slot was inserted. In that case
      // the rownr vector is empty and we use the first itsNrBl input rows.
      // Time related info can only be copied if not averaging and if the
      // the time slot was not missing.
      Table out(itsMS(rownrs));
      // Copy the input columns that do not change.
      writeMeta (out, buf);
      // Now write the data and flags.
      writeData (out, buf);
      // Flush if sufficient time slots are written.
      itsNrDone++;
      if (itsNrTimesFlush > 0  &&  itsNrDone%itsNrTimesFlush == 0) {
        itsMS.flush();
      }
      // Replace the rownrs in the buffer which is needed if in a later
      // step the MS gets updated.
      itsBuffer.setRowNrs (rownrs);
      getNextStep()->process(itsBuffer);
      return true;
    }

    void MSWriter::finish()
    {
      NSTimer::StartStop sstime(itsTimer);
      itsMS.flush();
      ///ROTiledStManAccessor acc1(itsMS, "TiledData");
      ///acc1.showCacheStatistics (cout);
      ///ROTiledStManAccessor acc2(itsMS, "TiledFlag");
      ///acc2.showCacheStatistics (cout);
      ///ROTiledStManAccessor acc3(itsMS, "TiledUVW");
      ///acc3.showCacheStatistics (cout);
      ///ROTiledStManAccessor acc4(itsMS, "TiledFullResFlag");
      ///acc4.showCacheStatistics (cout);
      // Create the VDS file.
      if (! itsClusterDesc.empty()) {
        string vdsName = itsMS.tableName() + ".vds";
        if (! itsVdsDir.empty()) {
          if (itsVdsDir[itsVdsDir.size() - 1] != '/') {
            itsVdsDir.append ("/");
          }
          vdsName = itsVdsDir + string(casacore::Path(vdsName).baseName());
        }
        // Create VDS file without detailed time info.
        DP3::VdsMaker::create (itsMS.tableName(), vdsName,
                                 itsClusterDesc, "", false);
      }
    }

    void MSWriter::addToMS (const string&)
    {
      getPrevStep()->addToMS(itsOutName);
    }

    void MSWriter::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      itsInterval     = info().timeInterval();
      itsNrCorr       = info().ncorr();
      itsNrChan       = info().nchan();
      itsNrBl         = info().nbaselines();
      itsNrTimes      = info().ntime();
      // Input can already be averaged, so take that into account.
      itsNChanAvg     = itsReader->nchanAvgFullRes() * info().nchanAvg();
      itsNTimeAvg     = itsReader->ntimeAvgFullRes() * info().ntimeAvg();
      // Create the MS.
      if (itsTileNChan <= 0  ||  itsTileNChan > getInfo().nchan()) {
        itsTileNChan = getInfo().nchan();
      }
      NSTimer::StartStop sstime(itsTimer);
      createMS (itsOutName, info(), itsTileSize, itsTileNChan);
      // Write the parset info into the history.
      writeHistory (itsMS, itsParset);
      itsMS.flush (true, true);
      DPLOG_INFO ("Finished preparing output MS", false);
      info().clearWrites();
      info().clearMetaChanged();
    }

    void MSWriter::show (std::ostream& os) const
    {
      os << "MSWriter " << itsName << std::endl;
      os << "  output MS:      " << itsMS.tableName() << std::endl;
      os << "  nchan:          " << itsNrChan << std::endl;
      os << "  ncorrelations:  " << itsNrCorr << std::endl;
      os << "  nbaselines:     " << itsNrBl << std::endl;
      os << "  ntimes:         " << itsNrTimes << std::endl;
      os << "  time interval:  " << itsInterval << std::endl;
      os << "  DATA column:    " << itsDataColName << std::endl;
      os << "  WEIGHT column:  " << itsWeightColName << std::endl;
      if (itsStManKeys.stManName == "dysco") {
        os
          << "  Compressed:     yes\n"
          << "  Data bitrate:   " << itsStManKeys.dyscoDataBitRate << '\n'
          << "  Weight bitrate: " << itsStManKeys.dyscoWeightBitRate << '\n'
          << "  Dysco mode:     " << itsStManKeys.dyscoNormalization << ' ' 
            << itsStManKeys.dyscoDistribution << '(' 
            << itsStManKeys.dyscoDistTruncation << ")\n";
      }
      else {
        os << "  Compressed:     no\n";
      }
    }

    void MSWriter::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " MSWriter " << itsName << endl;
    }

    void MSWriter::makeArrayColumn (ColumnDesc desc, const IPosition& ipos,
                                    DataManager* dm, Table& table, bool makeDirectColumn)
    {
      desc.setOptions(0);
      desc.setShape(ipos);
      if (makeDirectColumn) {
        desc.setOptions(ColumnDesc::Direct | ColumnDesc::FixedShape);
      }
      else {
        desc.setOptions(ColumnDesc::FixedShape);
      }
      if (table.tableDesc().isColumn(desc.name())) {
        table.removeColumn(desc.name());
      }
      // Use storage manager if given.
      if (dm == 0) {
        table.addColumn (desc);
      } else {
        table.addColumn (desc, *dm);
      }
    }

    void MSWriter::createMS (const string& outName, const DPInfo& info,
                             uint tileSize, uint tileNChan)
    {
      // Determine the data shape.
      IPosition dataShape(2, itsNrCorr, itsNrChan);
      // Obtain the MS description.
      TableDesc tdesc (itsReader->table().tableDesc());
      // Create the output table without the columns depending
      // on the nr of channels.
      // FLAG_CATEGORY is taken, but ignored when writing.
      Block<String> fixedColumns(20);
      fixedColumns[0]  = "UVW";
      fixedColumns[1]  = "FLAG_CATEGORY";
      fixedColumns[2]  = "WEIGHT";
      fixedColumns[3]  = "SIGMA";
      fixedColumns[4]  = "ANTENNA1";
      fixedColumns[5]  = "ANTENNA2";
      fixedColumns[6]  = "ARRAY_ID";
      fixedColumns[7]  = "DATA_DESC_ID";
      fixedColumns[8]  = "EXPOSURE";
      fixedColumns[9]  = "FEED1";
      fixedColumns[10] = "FEED2";
      fixedColumns[11] = "FIELD_ID";
      fixedColumns[12] = "FLAG_ROW";
      fixedColumns[13] = "INTERVAL";
      fixedColumns[14] = "OBSERVATION_ID";
      fixedColumns[15] = "PROCESSOR_ID";
      fixedColumns[16] = "SCAN_NUMBER";
      fixedColumns[17] = "STATE_ID";
      fixedColumns[18] = "TIME";
      fixedColumns[19] = "TIME_CENTROID";
      Table temptable = itsReader->table().project(fixedColumns);
      TableDesc newdesc = temptable.tableDesc();
      // Now quite some 'magic' is done to get the storage managers right.
      // Most of the columns do not change much and should be stored with
      // the IncrementalStMan. The ANTENNA columns change more often and
      // can best be stored with StandardStMan.
      // The new storage managers are only used for MSs stored with LofarStMan.
      // For 'normal' MSs they won't change.
      //
      // If needed, make WEIGHT and SIGMA a fixed shape, direct column.
      // In this way they are not written if the values are the same.
      {
        ColumnDesc& cdesc = newdesc.rwColumnDesc("WEIGHT");
        if (cdesc.shape().empty()) {
          cdesc.setShape (IPosition(1,itsNrCorr), true);
        }
      }
      {
        ColumnDesc& cdesc = newdesc.rwColumnDesc("SIGMA");
        if (cdesc.shape().empty()) {
          cdesc.setShape (IPosition(1,itsNrCorr), true);
        }
      }
      // Remove possible hypercolumn definitions.
      newdesc.adjustHypercolumns (SimpleOrderedMap<String,String>(String()));
      // Set data manager info.
      Record dminfo = temptable.dataManagerInfo();
      // Determine the DATA tile shape. Use all corrs and the given #channels.
      // The given tile size (in kbytes) determines the nr of rows in a tile .
      IPosition tileShape(3, itsNrCorr, tileNChan, 1);
      tileShape[2] = tileSize * 1024 / (8 * tileShape[0] * tileShape[1]);
      if (tileShape[2] < 1) {
        tileShape[2] = 1;
      }
      // Replace all non-writable storage managers (i.e. LofarStMan) by ISM.
      dminfo = DataManInfo::adjustStMan (dminfo, "IncrementalStMan");
      // Remove ANTENNA1 and ANTENNA2 from the dminfo.
      // Don't remove them if already stored with StandardStMan.
      Vector<String> removeCols(2);
      removeCols[0] = "ANTENNA1";
      removeCols[1] = "ANTENNA2";
      DataManInfo::removeDminfoColumns (dminfo, removeCols, "StandardStMan");
      // Use TiledStMan for UVW.
      // Use as many rows as used for the DATA columns, but minimal 1024.
      int tsmnrow = tileShape[2];
      if (tsmnrow < 1024) {
        tsmnrow = 1024;
      }
      DataManInfo::setTiledStMan (dminfo, Vector<String>(1, "UVW"),
                                  "TiledColumnStMan", "TiledUVW",
                                  IPosition(2, 3, tsmnrow));
      // Test if SSMVar already exists.
      bool hasSSMVar = false;
      for (uint i=0; i<dminfo.nfields(); ++i) {
        if (dminfo.subRecord(i).asString("NAME") == "SSMVar") {
          hasSSMVar = true;
          break;
        }
      }
      // Setup table creation. Exception is thrown if it exists already.
      Table::TableOption opt = itsOverwrite ? Table::New : Table::NewNoReplace;
      SetupNewTable newtab(outName, newdesc, opt);

      // First bind all columns to SSM.
      // For all columns defined in dminfo the bindings will be overwritten.
      // In this way variable columns like ANTENNA1/2 are bound to SSM.
      // Only do it if SSMVar does not exist (otherwise duplicate StMan name).
      if (!hasSSMVar) {
        StandardStMan ssm("SSMVar", 32768);
        newtab.bindAll (ssm);
      }

      // Bind all columns according to dminfo.
      newtab.bindCreate (dminfo);
      DataManagerCtor dyscoConstructor = 0;
      Record dyscoSpec;
      if (itsStManKeys.stManName == "dysco") {
        dyscoSpec = itsStManKeys.GetDyscoSpec();
        dyscoConstructor = DataManager::getCtor("DyscoStMan");
      }
      itsMS = Table(newtab);

      if (itsStManKeys.stManName == "dysco" && itsStManKeys.dyscoDataBitRate != 0) {
        // Add DATA column using Dysco stman.
        std::unique_ptr<DataManager> dyscoStMan(dyscoConstructor("DyscoData", dyscoSpec));
        makeArrayColumn (tdesc["DATA"], dataShape, dyscoStMan.get(), itsMS, true);
      }
      else {
        // Add DATA column using tsm.
        TiledColumnStMan tsm("TiledData", tileShape);
        makeArrayColumn (tdesc["DATA"], dataShape, &tsm, itsMS);
      }

      // Add FLAG column using tsm.
      // Use larger tile shape because flags are stored as bits.
      IPosition tileShapeF(tileShape);
      tileShapeF[2] *= 8;
      TiledColumnStMan tsmf("TiledFlag", tileShapeF);
      makeArrayColumn(tdesc["FLAG"], dataShape, &tsmf, itsMS);

      if (itsWriteFullResFlags) {
        // Add LOFAR_FULL_RES_FLAG column using tsm.
        // The input can already be averaged and averaging can be done in
        // this run, so the full resolution is the combination of both.
        uint orignchan = itsNrChan * itsNChanAvg;
        IPosition dataShapeF(2, (orignchan+7)/8, itsNTimeAvg);
        IPosition tileShapeF(3, (orignchan+7)/8, 1024, tileShape[2]);
        TiledColumnStMan tsmf("TiledFullResFlag", tileShapeF);
        ArrayColumnDesc<uChar> padesc("LOFAR_FULL_RES_FLAG",
                                      "flags in original full resolution",
                                      dataShapeF, ColumnDesc::FixedShape);
        makeArrayColumn(padesc, dataShapeF, &tsmf, itsMS);
      }
      if (itsStManKeys.stManName == "dysco" && itsStManKeys.dyscoWeightBitRate != 0) {
        // Add WEIGHT_SPECTRUM column using Dysco stman.
        std::unique_ptr<DataManager> dyscoStMan(dyscoConstructor(
          "DyscoWeightSpectrum", dyscoSpec) );
        ArrayColumnDesc<float> wsdesc("WEIGHT_SPECTRUM", "weight per corr/chan",
                                   dataShape, 
                                   ColumnDesc::FixedShape | ColumnDesc::Direct);
        makeArrayColumn (wsdesc, dataShape, dyscoStMan.get(), itsMS, true);
      }
      else {
        // Add WEIGHT_SPECTRUM column using tsm.
        TiledColumnStMan tsmw("TiledWeightSpectrum", tileShape);
        ArrayColumnDesc<float> wsdesc("WEIGHT_SPECTRUM", "weight per corr/chan",
                                      dataShape, ColumnDesc::FixedShape);
        makeArrayColumn(wsdesc, dataShape, &tsmw, itsMS);
      }
      // If present handle the CORRECTED_DATA and MODEL_DATA column.
      if (!tdesc.isColumn("CORRECTED_DATA")) {
        itsCopyCorrData = false;
      }
      if (!tdesc.isColumn("MODEL_DATA")) {
        itsCopyModelData = false;
      }
      if (itsCopyCorrData) {
        TiledColumnStMan tsmc("CorrectedData", tileShape);
        makeArrayColumn(tdesc["CORRECTED_DATA"], dataShape, &tsmc,
                        itsMS);

        IPosition iwShape(1, dataShape[1]);
        IPosition iwShapeTile(2, tileShape[1], tileShape[2]);
        TiledColumnStMan tsmw("TiledImagingWeight", iwShapeTile);
        ColumnDesc iwdesc(ArrayColumnDesc<float>("IMAGING_WEIGHT"));
        makeArrayColumn(iwdesc, iwShape, &tsmw, itsMS);
      }
      if (itsCopyModelData) {
        ColumnDesc mdesc = tdesc.columnDesc("MODEL_DATA");
        TableRecord& keyset = mdesc.rwKeywordSet();
        // Redefine possible keywords used by the CASA VisSet classes.
        if (keyset.isDefined("CHANNEL_SELECTION")) {
          keyset.removeField("CHANNEL_SELECTION");
        }
        Matrix<int> selection(2, 1);
        selection(0, 0) = 0;
        selection(1, 0) = itsNrChan;
        keyset.define("CHANNEL_SELECTION", selection);
        TiledColumnStMan tsmm("ModelData", tileShape);
        makeArrayColumn(mdesc, dataShape, &tsmm, itsMS);
      }
      DPLOG_INFO (" copying info and subtables ...", false);
      // Copy the info and subtables.
      TableCopy::copyInfo(itsMS, temptable);
      TableCopy::copySubTables(itsMS, temptable);
      // Adjust the SPECTRAL_WINDOW and DATA_DESCRIPTION table as needed.
      updateSpw (outName, info);
      // Adjust the OBSERVATION table as needed.
      updateObs (outName);
      // Adjust the FIELD table as needed.
      if (! info.phaseCenterIsOriginal()) {
        updatePhaseCentre (outName, info);
      }
      updateBeam (outName, "DATA", info);
    }

    void MSWriter::updateSpw (const string& outName, const DPInfo& info)
    {
      // Fix the SPECTRAL_WINDOW values by updating the values in the subtable.
      IPosition shape(1,itsNrChan);
      Table inSPW  = itsReader->table().keywordSet().asTable("SPECTRAL_WINDOW");
      Table outSPW = Table(outName + "/SPECTRAL_WINDOW", Table::Update);
      Table outDD  = Table(outName + "/DATA_DESCRIPTION", Table::Update);
      assert(outSPW.nrow() == outDD.nrow());
      uint spw = itsReader->spectralWindow();
      // Remove all rows before and after the selected band.
      // Do it from the end, otherwise row numbers change.
      for (uint i=outSPW.nrow(); i>0;) {
        if (--i != spw) {
          outSPW.removeRow (i);
          outDD.removeRow (i);
        }
      }
      assert (outSPW.nrow() == 1);
      // Set nr of channels.
      ScalarColumn<Int> channum(outSPW, "NUM_CHAN");
      channum.fillColumn (itsNrChan);
      // Change the column shapes.
      TableDesc tdesc = inSPW.tableDesc();
      makeArrayColumn (tdesc["CHAN_FREQ"], shape, 0, outSPW);
      makeArrayColumn (tdesc["CHAN_WIDTH"], shape, 0, outSPW);
      makeArrayColumn (tdesc["EFFECTIVE_BW"], shape, 0, outSPW);
      makeArrayColumn (tdesc["RESOLUTION"], shape, 0, outSPW);
      // Create the required column objects.
      ArrayColumn<Double> outFREQ(outSPW, "CHAN_FREQ");
      ArrayColumn<Double> outWIDTH(outSPW, "CHAN_WIDTH");
      ArrayColumn<Double> outBW(outSPW, "EFFECTIVE_BW");
      ArrayColumn<Double> outRESOLUTION(outSPW, "RESOLUTION");
      ScalarColumn<Double> outTOTALBW(outSPW, "TOTAL_BANDWIDTH");
      ScalarColumn<Double> outREFFREQ(outSPW, "REF_FREQUENCY");
      outFREQ.put      (0, info.chanFreqs());
      outWIDTH.put     (0, info.chanWidths());
      outBW.put        (0, info.effectiveBW());
      outRESOLUTION.put(0, info.resolutions());
      outTOTALBW.put   (0, info.totalBW());
      outREFFREQ.put   (0, info.refFreq());
      // Adjust the spwid in the DATA_DESCRIPTION.
      ScalarColumn<Int> spwCol(outDD, "SPECTRAL_WINDOW_ID");
      spwCol.put (0, 0);
    }

    void MSWriter::updateObs (const string& outName)
    {
      Table outObs = Table(outName + "/OBSERVATION", Table::Update);
      // Set nr of channels.
      ArrayColumn<Double> timeRange(outObs, "TIME_RANGE");
      Vector<double> times(2);
      times[0] = itsReader->firstTime() - 0.5 * itsReader->getInfo().timeInterval();
      times[1] = itsReader->lastTime()  + 0.5 * itsReader->getInfo().timeInterval();
      // There should be one row, but loop in case of.
      for (uint i=0; i<outObs.nrow(); ++i) {
        timeRange.put (i, times);
      }
    }

    void MSWriter::updatePhaseCentre (const string& outName, const DPInfo& info)
    {
      Table outField = Table(outName + "/FIELD", Table::Update);
      // Write new phase center.
      ArrayMeasColumn<MDirection> phaseCol (outField, "PHASE_DIR");
      Vector<MDirection> dir(1, info.phaseCenter());
      phaseCol.put (0, dir);
    }
    
    void MSWriter::updateBeam (const std::string& outName, const std::string& outColName, const DPInfo& info)
    {
      const char
        *beamModeFieldName = "LOFAR_APPLIED_BEAM_MODE",
        *beamDirFieldName = "LOFAR_APPLIED_BEAM_DIR";
      
      Table mainTable(outName, Table::Update);
      ArrayColumn<Complex> dataColumn(mainTable, outColName);
      bool fieldsExist = dataColumn.keywordSet().isDefined(beamModeFieldName);
      std::string modeStr;
      switch(info.beamCorrectionMode())
      {
        default:
        case NoBeamCorrection: modeStr = "None"; break;
        case ElementBeamCorrection: modeStr = "Element"; break;
        case ArrayFactorBeamCorrection: modeStr = "ArrayFactor"; break;
        case FullBeamCorrection: modeStr = "Full"; break;
      }
      // If no beam correction has been applied and the LOFAR beam fields don't
      // exist, we have to do nothing (no fields implies no beam correction).
      // If they do exist, we have to make sure they are set to indicate
      // no beam correction.
      if(fieldsExist || info.beamCorrectionMode() != NoBeamCorrection)
      {
        dataColumn.rwKeywordSet().define(beamModeFieldName, modeStr);
        Record record;
        MeasureHolder(info.beamCorrectionDir()).toRecord(record);
        dataColumn.rwKeywordSet().defineRecord(beamDirFieldName, record);
      }
    }

    void MSWriter::writeHistory (Table& ms, const ParameterSet& parset)
    {
      Table histtab(ms.keywordSet().asTable("HISTORY"));
      histtab.reopenRW();
      ScalarColumn<Double> time        (histtab, "TIME");
      ScalarColumn<Int>    obsId       (histtab, "OBSERVATION_ID");
      ScalarColumn<String> message     (histtab, "MESSAGE");
      ScalarColumn<String> application (histtab, "APPLICATION");
      ScalarColumn<String> priority    (histtab, "PRIORITY");
      ScalarColumn<String> origin      (histtab, "ORIGIN");
      ArrayColumn<String>  parms       (histtab, "APP_PARAMS");
      ArrayColumn<String>  cli         (histtab, "CLI_COMMAND");
      // Put all parset entries in a Vector<String>.
      // Some WSRT MSs have a FixedShape APP_PARAMS and CLI_COMMAND column.
      // For them, put all params in a single vector element (with newlines).
      bool fixedShaped = 
        (parms.columnDesc().options() & ColumnDesc::FixedShape) != 0;
      Vector<String> appvec;
      Vector<String> clivec;
      if (fixedShaped) {
        appvec.resize(1);
        clivec.resize(1);
        ostringstream ostr;
        parset.writeStream (ostr);
        appvec[0] = ostr.str();
      } else {
        appvec.resize (parset.size());
        Array<String>::contiter viter = appvec.cbegin();
        for (ParameterSet::const_iterator iter = parset.begin();
             iter != parset.end(); ++iter, ++viter) {
          *viter = iter->first + '=' + iter->second.get();
        }
      }
      uint rownr = histtab.nrow();
      histtab.addRow();
      time.put        (rownr, Time().modifiedJulianDay()*24.*3600.);
      obsId.put       (rownr, 0);
      message.put     (rownr, "parameters");
      application.put (rownr, "NDPPP");
      priority.put    (rownr, "NORMAL");
      origin.put      (rownr, DPPPVersion::AsString());
      parms.put       (rownr, appvec);
      cli.put         (rownr, clivec);
    }

    void MSWriter::writeData (Table& out, const DPBuffer& buf)
    {
      ArrayColumn<Complex> dataCol(out, itsDataColName);
      ArrayColumn<Bool>    flagCol(out, "FLAG");
      ScalarColumn<Bool>   flagRowCol(out, "FLAG_ROW");

      if (buf.getData().empty()) {
        return;
      }

      // Write WEIGHT_SPECTRUM and DATA
      ArrayColumn<Float> weightCol(out, "WEIGHT_SPECTRUM");
      itsBuffer.referenceFilled (buf);
      const Array<Float>& weights = itsReader->fetchWeights (buf, itsBuffer,
                                                             itsTimer);

      // If compressing, flagged values need to be set to NaN, and flagged
      // weights to zero, to decrease the dynamic range
      if (itsStManKeys.stManName == "dysco") {
        Cube<Complex> dataCopy = buf.getData().copy();
        Cube<Complex>::iterator dataIter = dataCopy.begin();
        Cube<float> weightsCopy = weights.copy();
        Cube<float>::iterator weightsIter = weightsCopy.begin();
        for (Cube<bool>::const_iterator flagIter = buf.getFlags().begin();
             flagIter != buf.getFlags().end(); ++flagIter) {
          if (*flagIter) {
            *dataIter = Complex(std::numeric_limits<float>::quiet_NaN(),
                                std::numeric_limits<float>::quiet_NaN());
            *weightsIter = 0.;
          }
          ++dataIter;
          ++weightsIter;
        }
        dataCol.putColumn (dataCopy);
        weightCol.putColumn (weightsCopy);
      }
      else {
        dataCol.putColumn (buf.getData());
        weightCol.putColumn (weights);
      }
      
      flagCol.putColumn (buf.getFlags());
      // A row is flagged if no flags in the row are False.
      Vector<Bool> rowFlags (partialNFalse(buf.getFlags(), IPosition(2,0,1)) == 0u);
      flagRowCol.putColumn (rowFlags);
      if (itsWriteFullResFlags) {
        writeFullResFlags (out, buf);
      }

      // Write UVW
      ArrayColumn<Double> uvwCol(out, "UVW");
      const Array<Double>& uvws = itsReader->fetchUVW (buf, itsBuffer,
                                                       itsTimer);
      uvwCol.putColumn (uvws);
    }

    void MSWriter::writeFullResFlags (Table& out, const DPBuffer& buf)
    {
      // Get the flags.
      const Cube<bool>& flags = itsReader->fetchFullResFlags (buf, itsBuffer,
                                                              itsTimer);
      const IPosition& ofShape = flags.shape();
      if(uint(ofShape[0]) != itsNChanAvg * itsNrChan)
        throw Exception(std::to_string(ofShape[0]) + std::to_string(itsNChanAvg) + '*' + std::to_string(itsNrChan));
      if(uint(ofShape[1]) != itsNTimeAvg)
        throw Exception(std::to_string(ofShape[1]) + std::to_string(itsNTimeAvg));
      // Convert the bools to uChar bits.
      IPosition chShape(ofShape);
      chShape[0] = (ofShape[0] + 7) / 8;
      Cube<uChar> chars(chShape);
      // If their sizes match, do it all in one go.
      // Otherwise we have to iterate.
      if (ofShape[0] == chShape[0]*8) {
        Conversion::boolToBit (chars.data(), flags.data(), flags.size());
      } else {
        assert(ofShape[0] < chShape[0]*8);
        const bool* flagsPtr = flags.data();
        uChar* charsPtr = chars.data();
        for (int i=0; i<ofShape[1]*ofShape[2]; ++i) {
          Conversion::boolToBit (charsPtr, flagsPtr, ofShape[0]);
          flagsPtr += ofShape[0];
          charsPtr += chShape[0];
        }
      }
      ArrayColumn<uChar> fullResCol(out, "LOFAR_FULL_RES_FLAG");
      if (! fullResCol.keywordSet().isDefined ("NCHAN_AVG")) {
        fullResCol.rwKeywordSet().define ("NCHAN_AVG", int(itsNChanAvg));
        fullResCol.rwKeywordSet().define ("NTIME_AVG", int(itsNTimeAvg));
      }
      fullResCol.putColumn (chars);
    } 

    void MSWriter::writeMeta (Table& out, const DPBuffer& buf)
    {
      // Fill ANTENNA1/2.
      ScalarColumn<Int> ant1col(out, "ANTENNA1");
      ScalarColumn<Int> ant2col(out, "ANTENNA2");
      ant1col.putColumn (getInfo().getAnt1());
      ant2col.putColumn (getInfo().getAnt2());
      // Fill all rows that do not change.
      fillSca<Double> (buf.getTime(), out, "TIME");
      fillSca<Double> (buf.getTime(), out, "TIME_CENTROID");
      fillSca<Double> (buf.getExposure(), out, "EXPOSURE");
      fillSca<Double> (itsInterval, out, "INTERVAL");
      fillSca<Int> (0, out, "FEED1");
      fillSca<Int> (0, out, "FEED2");
      fillSca<Int> (0, out, "DATA_DESC_ID");
      fillSca<Int> (0, out, "PROCESSOR_ID");
      fillSca<Int> (0, out, "FIELD_ID");
      fillSca<Int> (0, out, "SCAN_NUMBER");
      fillSca<Int> (0, out, "ARRAY_ID");
      fillSca<Int> (0, out, "OBSERVATION_ID");
      fillSca<Int> (0, out, "STATE_ID");
      Array<Float> arr(IPosition(1, itsNrCorr));
      arr = 1;
      fillArr<Float> (arr, out, "SIGMA");
      fillArr<Float> (arr, out, "WEIGHT");
    }

    void MSWriter::copyMeta (const Table& in, Table& out, bool copyTimeInfo)
    {
      // Copy all rows that do not change.
      copySca<int> (in, out, "ANTENNA1");
      copySca<int> (in, out, "ANTENNA2");
      copySca<int> (in, out, "FEED1");
      copySca<int> (in, out, "FEED2");
      copySca<int> (in, out, "PROCESSOR_ID");
      copySca<int> (in, out, "FIELD_ID");
      copySca<int> (in, out, "SCAN_NUMBER");
      copySca<int> (in, out, "ARRAY_ID");
      copySca<int> (in, out, "OBSERVATION_ID");
      copySca<int> (in, out, "STATE_ID");
      copyArr<float> (in, out, "SIGMA");
      copyArr<float> (in, out, "WEIGHT");
      if (copyTimeInfo) {
        copySca<double>(in, out, "TIME");
        copySca<double>(in, out, "TIME_CENTROID");
        copySca<double>(in, out, "INTERVAL");
        copySca<double>(in, out, "EXPOSURE");
        copyArr<double>(in, out, "UVW");
      }
    }

  } //# end namespace
}
