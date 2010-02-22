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

#include <lofar_config.h>
#include <DPPP/Package__Version.h>
#include <DPPP/MSWriter.h>
#include <DPPP/MSUpdater.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/AverageInfo.h>
#include <Common/ParameterSet.h>
#include <tables/Tables/TableCopy.h>
#include <tables/Tables/DataManInfo.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/StandardStMan.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Containers/Record.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSWriter::MSWriter (MSReader* reader, const std::string& outName,
                        const AverageInfo& avgInfo,
                        const ParameterSet& parset, const string& prefix)
      : itsReader       (reader),
        itsInterval     (avgInfo.timeInterval()),
        itsCopyTimeInfo (avgInfo.ntimeAvg() == 1),
        itsNrCorr       (reader->ncorr()),
        itsNrChan       (avgInfo.nchan()),
        itsNrBl         (reader->nbaselines()),
        itsNrTimes      (avgInfo.ntime()),
        itsOrigNrChan   (avgInfo.origNChan()),
        itsNTimeAvg     (avgInfo.ntimeAvg())
    {
      // Get tilesize (default 1024 KBytes).
      uint tileSize       = parset.getUint (prefix+"tilesize", 1024);
      uint tileNChan      = parset.getUint (prefix+"tilenchan", 8);
      itsCopyCorrData     = parset.getBool (prefix+"copycorrecteddata", false);
      itsCopyModelData    = parset.getBool (prefix+"copymodeldata", false);
      itsWritePreAvgFlags = parset.getBool (prefix+"writefullresflag", true);
      itsDataColName      = parset.getString (prefix+"datacolumn", "DATA");
      // Create the MS.
      createMS (outName, avgInfo, tileSize, tileNChan);
      // Write the parset info into the history.
      writeHistory (itsMS, parset);
      itsMS.flush (true, true);
      cout << "Finished preparing output MS" << endl;
    }

    MSWriter::~MSWriter()
    {}

    bool MSWriter::process (const DPBuffer& buf)
    {
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
      Table in;
      bool copyTimeInfo = itsCopyTimeInfo;
      if (buf.getRowNrs().empty()) {
        in = itsReader->table()(itsReader->getBaseRowNrs());
        copyTimeInfo = false;
      } else {
        in = itsReader->table()(buf.getRowNrs());
      }
      // Copy the input columns that do not change.
      copyMeta (in, out, copyTimeInfo);
      // Write the time related values if not already copied.
      if (!copyTimeInfo) {
        writeTimeInfo (out, buf.getTime());
      }
      // Now write the data and flags.
      writeData (out, buf);
      return true;
    }

    void MSWriter::finish()
    {
      itsMS.flush();
    }

    void MSWriter::show (std::ostream& os)
    {
      os << "MSWriter" << std::endl;
      os << "  output MS:      " << itsMS.tableName() << std::endl;
      os << "  nchan:          " << itsNrChan << std::endl;
      os << "  ncorrelations:  " << itsNrCorr << std::endl;
      os << "  nbaselines:     " << itsNrBl << std::endl;
      os << "  ntimes:         " << itsNrTimes << std::endl;
      os << "  time interval:  " << itsInterval << std::endl;
      os << "  DATA column:    " << itsDataColName << std::endl;
    }

    void MSWriter::makeArrayColumn (ColumnDesc desc, const IPosition& ipos,
                                    TiledColumnStMan* tsm, Table& table)
    {
      desc.setOptions(0);
      desc.setShape(ipos);
      desc.setOptions(ColumnDesc::FixedShape);
      if (table.tableDesc().isColumn(desc.name())) {
        table.removeColumn(desc.name());
      }
      // Use tiled storage manager if given.
      if (tsm == 0) {
        table.addColumn (desc);
      } else {
        table.addColumn (desc, *tsm);
      }
    }

    void MSWriter::createMS (const string& outName, const AverageInfo& avgInfo,
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
      // Setup table creation. Exception is thrown if it exists already.
      SetupNewTable newtab(outName, newdesc, Table::NewNoReplace);
      // First bind all column to SSM.
      // For all columns defined in dminfo the bindings will be overwritten.
      // In this way variable columns like ANTENNA1/2 are bound to SSM.
      StandardStMan ssm("SSMVar", 32768);
      newtab.bindAll (ssm);
      // Bind all columns according to dminfo.
      newtab.bindCreate (dminfo);
      itsMS = Table(newtab);
      {
        // Add DATA column using tsm.
        TiledColumnStMan tsm("TiledData", tileShape);
        makeArrayColumn (tdesc["DATA"], dataShape, &tsm, itsMS);
      }
      {
        // Add FLAG column using tsm.
        // Use larger tile shape because flags are stored as bits.
        IPosition tileShapeF(tileShape);
        tileShapeF[2] *= 8;
        TiledColumnStMan tsmf("TiledFlag", tileShapeF);
        makeArrayColumn(tdesc["FLAG"], dataShape, &tsmf, itsMS);
      }
      if (itsWritePreAvgFlags) {
        // Add LOFAR_FULL_RES_FLAG column using tsm.
        IPosition dataShapeF(2, itsOrigNrChan/8, itsNTimeAvg);
        IPosition tileShapeF(2, itsOrigNrChan/8, 1024);
        TiledColumnStMan tsmf("TiledFullResFlag", tileShapeF);
        ArrayColumnDesc<uChar> padesc("LOFAR_FULL_RES_FLAG",
                                      "flags in original full resolution",
                                      dataShapeF, ColumnDesc::FixedShape);
        makeArrayColumn(padesc, dataShapeF, &tsmf, itsMS);
      }
      {
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
        Matrix<Int> selection(2, 1);
        selection(0, 0) = 0;
        selection(1, 0) = itsNrChan;
        keyset.define("CHANNEL_SELECTION", selection);
        TiledColumnStMan tsmm("ModelData", tileShape);
        makeArrayColumn(mdesc, dataShape, &tsmm, itsMS);
      }
      cout << " copying info and subtables ..." << endl;
      // Copy the info and subtables.
      TableCopy::copyInfo(itsMS, temptable);
      TableCopy::copySubTables(itsMS, temptable);
      updateSpw (outName, avgInfo);
    }

    void MSWriter::updateSpw (const string& outName, const AverageInfo& avgInfo)
    {
      // Fix the SPECTRAL_WINDOW values by updating the values in the subtable.
      IPosition shape(1,itsNrChan);
      Table inSPW  = itsReader->table().keywordSet().asTable("SPECTRAL_WINDOW");
      Table outSPW = Table(outName + "/SPECTRAL_WINDOW", Table::Update);
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
      ROArrayColumn<Double> inFREQ(inSPW, "CHAN_FREQ");
      ROArrayColumn<Double> inWIDTH(inSPW, "CHAN_WIDTH");
      ROArrayColumn<Double> inBW(inSPW, "EFFECTIVE_BW");
      ROArrayColumn<Double> inRESOLUTION(inSPW, "RESOLUTION");
      ArrayColumn<Double> outFREQ(outSPW, "CHAN_FREQ");
      ArrayColumn<Double> outWIDTH(outSPW, "CHAN_WIDTH");
      ArrayColumn<Double> outBW(outSPW, "EFFECTIVE_BW");
      ArrayColumn<Double> outRESOLUTION(outSPW, "RESOLUTION");
      ScalarColumn<Double> outTOTALBW(outSPW, "TOTAL_BANDWIDTH");
      Vector<double> newFreq  (itsNrChan);
      Vector<double> newWidth (itsNrChan);
      Vector<double> newBW    (itsNrChan);
      Vector<double> newRes   (itsNrChan);
      // Loop through all rows.
      for (uint i=0; i<inSPW.nrow(); ++i) {
        Vector<double> oldFreq = inFREQ(i);
        Vector<double> oldWidth = inWIDTH(i);
        Vector<double> oldBW = inBW(i);
        Vector<double> oldRes = inRESOLUTION(i);
        double totalBW = 0;
        uint first = avgInfo.startChan();
        // This loops assumes regularly spaced, adjacent frequency channels.
        for (uint j=0; j<itsNrChan; ++j) { 
          uint last  = first + avgInfo.nchanAvg();
          if (last > avgInfo.startChan() + avgInfo.origNChan()) {
            last = avgInfo.startChan() + avgInfo.origNChan();
          }
          double sf = oldFreq[first]  - 0.5*oldWidth[first];
          double ef = oldFreq[last-1] + 0.5*oldWidth[last-1];
          newFreq[j]  = 0.5 * (sf + ef);
          newWidth[j] = ef - sf;
          double newbw = 0;
          double newres = 0;
          for (uint k=first; k<last; ++k) {
            newbw  += oldBW[k];
            newres += oldRes[k];
          }
          newBW[j]  = newbw;
          newRes[j] = newres;
          totalBW  += newbw;
          first = last;
        }
        outFREQ.put      (i, newFreq);
        outWIDTH.put     (i, newWidth);
        outBW.put        (i, newBW);
        outRESOLUTION.put(i, newRes);
        outTOTALBW.put   (i, totalBW);
      }
    }

    void MSWriter::writeHistory (Table& ms, const ParameterSet& parset)
    {
      Table histtab(ms.keywordSet().asTable("HISTORY"));
      histtab.reopenRW();
      ostringstream allParms;
      parset.writeStream (allParms);
      uInt rownr = histtab.nrow();
      histtab.addRow();
      ScalarColumn<double> time        (histtab, "TIME");
      ScalarColumn<Int>    obsId       (histtab, "OBSERVATION_ID");
      ScalarColumn<String> message     (histtab, "MESSAGE");
      ScalarColumn<String> priority    (histtab, "PRIORITY");
      ScalarColumn<String> origin      (histtab, "ORIGIN");
      ScalarColumn<String> application (histtab, "APPLICATION");
      time.put        (rownr, Time().modifiedJulianDay()*24.*3600.);
      obsId.put       (rownr, 0);
      application.put (rownr, "DPPP");
      message.put     (rownr, allParms.str());
      priority.put    (rownr, "NORMAL");
      origin.put      (rownr, Version::getInfo<DPPPVersion>("DPPP", "full"));
    }

    void MSWriter::writeTimeInfo (Table& out, double time)
    {
      ScalarColumn<double> timeCol (itsMS, "TIME");
      ScalarColumn<double> tcenCol (itsMS, "TIME_CENTROID");
      ScalarColumn<double> intvCol (itsMS, "INTERVAL");
      ScalarColumn<double> expoCol (itsMS, "EXPOSURE");
      ArrayColumn<double>  uvwCol  (itsMS, "UVW");
      for (uint i=0; i<out.nrow(); ++i) {
        timeCol.put (i, time);
        tcenCol.put (i, time);
        intvCol.put (i, itsInterval);
        expoCol.put (i, itsInterval);
      }
    }

    void MSWriter::writeData (Table& out, const DPBuffer& buf)
    {
      ArrayColumn<Complex> dataCol(out, itsDataColName);
      ArrayColumn<bool>    flagCol(out, "FLAG");
      dataCol.putColumn (buf.getData());
      flagCol.putColumn (buf.getFlags());
      if (itsWritePreAvgFlags) {
        writePreAvgFlags (out, buf);
      }
      ArrayColumn<float> weightCol(out, "WEIGHT_SPECTRUM");
      weightCol.putColumn (itsReader->fetchWeights (buf, buf.getRowNrs()));
    }

    void MSWriter::writePreAvgFlags (Table& out, const DPBuffer& buf)
    {
      // Get the flags.
      Cube<bool> flags (itsReader->fetchPreAvgFlags (buf, buf.getRowNrs()));
      const IPosition& ofShape = flags.shape();
      // Convert the bools to uChar bits.
      IPosition chShape(ofShape);
      chShape[0] = (ofShape[0] + 7) / 8;
      Cube<uChar> chars(chShape);
      // If their sizes match, do it all in one go.
      // Otherwise we have to iterate.
      if (ofShape[0] == chShape[0]*8) {
        Conversion::boolToBit (chars.data(), flags.data(), flags.size());
      } else {
        ASSERT (ofShape[0] < chShape[0]*8);
        const bool* flagsPtr = flags.data();
        uChar* charsPtr = chars.data();
        for (int i=0; i<ofShape[2]*ofShape[3]; ++i) {
          Conversion::bitToBool (charsPtr, flagsPtr, ofShape[0]);
          flagsPtr += ofShape[0];
          charsPtr += chShape[0];
        }
      }
      ArrayColumn<uChar> preAvgCol(out, "LOFAR_FULL_RES_FLAG");
      preAvgCol.putColumn (chars);
    } 

    void MSWriter::copyMeta (const Table& in, Table& out, bool copyTimeInfo)
    {
      // Copy all rows that do not change.
      copySca<int> (in, out, "ANTENNA1");
      copySca<int> (in, out, "ANTENNA2");
      copySca<int> (in, out, "FEED1");
      copySca<int> (in, out, "FEED2");
      copySca<int> (in, out, "DATA_DESC_ID");
      copySca<int> (in, out, "PROCESSOR_ID");
      copySca<int> (in, out, "FIELD_ID");
      copySca<int> (in, out, "SCAN_NUMBER");
      copySca<int> (in, out, "ARRAY_ID");
      copySca<int> (in, out, "OBSERVATION_ID");
      copySca<int> (in, out, "STATE_ID");
      copySca<bool>(in, out, "FLAG_ROW");
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
