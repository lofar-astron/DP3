//# MSReader.cc: DPPP step reading from an MS
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
#include <DPPP/MSReader.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/DPLogger.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/RecordGram.h>
#include <measures/Measures/MeasTable.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/TableMeasures/ArrayMeasColumn.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSSelection.h>
#include <casa/Containers/Record.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Quanta/MVTime.h>
#include <casa/OS/Conversion.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSReader::MSReader()
      : itsReadVisData (False),
        itsLastMSTime  (0),
        itsNrRead      (0),
        itsNrInserted  (0)
    {}

    MSReader::MSReader (const string& msName,
                        const ParSet& parset, const string& prefix,
                        bool missingData)
      : itsReadVisData (False),
        itsMissingData (missingData),
        itsLastMSTime  (0),
        itsNrRead      (0),
        itsNrInserted  (0)
    {
      NSTimer::StartStop sstime(itsTimer);
      // Get info from parset.
      itsSpw              = parset.getInt    (prefix+"band", -1);
      itsStartChanStr     = parset.getString (prefix+"startchan", "0");
      itsNrChanStr        = parset.getString (prefix+"nchan", "0");
      string startTimeStr = parset.getString (prefix+"starttime", "");
      string endTimeStr   = parset.getString (prefix+"endtime", "");
      itsUseFlags         = parset.getBool   (prefix+"useflag", true);
      itsDataColName      = parset.getString (prefix+"datacolumn", "DATA");
      itsAutoWeight       = parset.getBool   (prefix+"autoweight", false);
      itsAutoWeightForce  = parset.getBool   (prefix+"forceautoweight", false);
      itsNeedSort         = parset.getBool   (prefix+"sort", false);
      itsSelBL            = parset.getString (prefix+"baseline", string());
      // Try to open the MS and get its full name.
      if (itsMissingData  &&  !Table::isReadable (msName)) {
        DPLOG_WARN_STR ("MeasurementSet " << msName
			<< " not found; dummy data used");
        return;
      }
      itsMS = MeasurementSet (msName, TableLock::AutoNoReadLocking);
      itsMSName = itsMS.tableName();
      // See if a selection on band needs to be done.
      // We assume that DATA_DESC_ID and SPW_ID map 1-1.
      if (itsSpw >= 0) {
        Table subset = itsMS (itsMS.col("DATA_DESC_ID") == itsSpw);
        // If not all is selected, use the selection.
        if (subset.nrow() < itsMS.nrow()) {
          ASSERTSTR (subset.nrow() > 0, "Band " << itsSpw << " not found in "
                     << itsMSName);
          itsMS = subset;
        }
      } else {
        itsSpw = 0;
      }
      // See if a selection on baseline needs to be done.
      if (! itsSelBL.empty()) {
        MSSelection select;
        // Set given selection strings.
        select.setAntennaExpr (itsSelBL);
        // Create a table expression for an MS representing the selection.
        MeasurementSet ms(itsMS);
        TableExprNode node = select.toTableExprNode (&ms);
        Table subset = itsMS(node);
        // If not all is selected, use the selection.
        if (subset.nrow() < itsMS.nrow()) {
          ASSERTSTR (subset.nrow() > 0, "Baselines " << itsSelBL
                     << "not found in " << itsMSName);
          itsMS = subset;
        }
      }
      // Prepare the MS access and get time info.
      double startTime, endTime;
      prepare (startTime, endTime, itsTimeInterval);
      // Start and end time can be given in the parset in case leading
      // or trailing time slots are missing.
      // They can also be used to select part of the MS.
      Quantity qtime;
      itsFirstTime = startTime;
      if (!startTimeStr.empty()) {
        if (!MVTime::read (qtime, startTimeStr)) {
          THROW (LOFAR::Exception, startTimeStr << " is an invalid date/time");
        }
        itsFirstTime = qtime.getValue("s");
        ASSERT (itsFirstTime <= endTime);
      }
      itsLastTime = endTime;
      if (!endTimeStr.empty()) {
        if (!MVTime::read (qtime, endTimeStr)) {
          THROW (LOFAR::Exception, endTimeStr << " is an invalid date/time");
        }
        itsLastTime = qtime.getValue("s");
      }
      ASSERT (itsLastTime >= itsFirstTime);
      // If needed, skip the first times in the MS.
      // It also sets itsFirstTime properly (round to time/interval in MS).
      skipFirstTimes();
      itsNextTime  = itsFirstTime;
      itsStartTime = itsFirstTime - 0.5*itsTimeInterval;
      // Parse the chan expressions.
      // Nr of channels can be used as 'nchan' in the expressions.
      Record rec;
      rec.define ("nchan", itsNrChan);
      TableExprNode node1 (RecordGram::parse(rec, itsStartChanStr));
      TableExprNode node2 (RecordGram::parse(rec, itsNrChanStr));
      // nchan=0 means until the last channel.
      double result;
      node1.get (rec, result);
      itsStartChan = uint(result+0.001);
      node2.get (rec, result);
      uint nrChan = uint(result+0.0001);
      uint nAllChan = itsNrChan;
      ASSERTSTR (itsStartChan < nAllChan,
                 "startchan " << itsStartChan
                 << " exceeds nr of channels in MS (" << nAllChan << ')');
      uint maxNrChan = nAllChan - itsStartChan;
      if (nrChan == 0) {
        itsNrChan = maxNrChan;
      } else {
        itsNrChan = std::min (nrChan, maxNrChan);
      }
      // Are all channels used?
      itsUseAllChan = itsStartChan==0 && itsNrChan==nAllChan;
      // Do the rest of the preparation.
      prepare2();
      // Take subset of channel frequencies if needed.
      // Make sure to copy the subset to get a proper Vector.
      // Form the slicer to get channels and correlations from column.
      itsColSlicer = Slicer(IPosition(2, 0, itsStartChan),
                            IPosition(2, itsNrCorr, itsNrChan));
      // Form the slicer to get channels, corrs, and baselines from array.
      itsArrSlicer = Slicer(IPosition(3, 0, itsStartChan, 0),
                            IPosition(3, itsNrCorr, itsNrChan, itsNrBl));
      // Initialize the flag counters.
      itsFlagCounter.init (getInfo());
    }

    MSReader::~MSReader()
    {}

    void MSReader::updateInfo (const DPInfo&)
    {}

    casa::String MSReader::msName() const
    {
      return itsMSName;
    }

    void MSReader::setReadVisData (bool readVisData)
    {
      itsReadVisData = readVisData;
    }

    bool MSReader::process (const DPBuffer&)
    {
      {
        NSTimer::StartStop sstime(itsTimer);
        itsBuffer.clear();
        // Use time from the current time slot in the MS.
        bool useIter = false;
        while (!itsIter.pastEnd()) {
          // Take time from row 0 in subset.
          double mstime = ROScalarColumn<double>(itsIter.table(), "TIME")(0);
          // Skip time slot and give warning if MS data is not in time order.
          if (mstime < itsLastMSTime) {
            LOG_WARN_STR ("Time at rownr "
                          << itsIter.table().rowNumbers(itsMS)[0]
                          << " of MS " << itsMS.tableName()
                          << " is less than previous time slot");
          } else {
            // Use the time slot if near or < nexttime, but > starttime.
            // In this way we cater for irregular times in some WSRT MSs.
            if (near(mstime, itsNextTime)  ||
                (mstime > itsFirstTime  &&  mstime < itsNextTime)) {
              itsFirstTime -= itsNextTime-mstime;
              itsNextTime = mstime;
              useIter = true;
              break;
            }
            if (mstime > itsNextTime) {
              // A time slot seems to be missing; insert one.
              break;
            }
          }
          // Skip this time slot.
          itsLastMSTime = mstime;
          itsIter.next();
        }
        // Stop if at the end.
        if (itsNextTime > itsLastTime  &&  !near(itsNextTime, itsLastTime)) {
          return false;
        }
        // Fill the buffer.
        itsBuffer.setTime (itsNextTime);
        ///cout << "read time " <<itsBuffer.getTime() - 4472025855.0<<endl;
        if (!useIter) {
          // Need to insert a fully flagged time slot.
          itsBuffer.setRowNrs (Vector<uint>());
          itsBuffer.setExposure (itsTimeInterval);
          itsBuffer.getData().resize  (itsNrCorr, itsNrChan, itsNrBl);
          itsBuffer.getFlags().resize (itsNrCorr, itsNrChan, itsNrBl);
          itsBuffer.getData() = Complex();
          itsBuffer.getFlags() = true;
          // Calculate UVWs for them. Setup UVW object if not done yet.
          calcUVW();
          itsNrInserted++;
        } else {
          itsBuffer.setRowNrs (itsIter.table().rowNumbers(itsMS));
          if (itsMissingData) {
            // Data column not present, so fill a fully flagged time slot.
            itsBuffer.setExposure (itsTimeInterval);
            itsBuffer.getData().resize  (itsNrCorr, itsNrChan, itsNrBl);
            itsBuffer.getFlags().resize (itsNrCorr, itsNrChan, itsNrBl);
            itsBuffer.getData() = Complex();
            itsBuffer.getFlags() = true;
          } else {
            // Set exposure.
            itsBuffer.setExposure (ROScalarColumn<double>
                                   (itsIter.table(), "EXPOSURE")(0));
            // Get data and flags from the MS.
            if (itsReadVisData) {
              ROArrayColumn<Complex> dataCol(itsIter.table(), itsDataColName);
              if (itsUseAllChan) {
                itsBuffer.setData (dataCol.getColumn());
              } else {
                itsBuffer.setData (dataCol.getColumn(itsColSlicer));
              }
            }
            if (itsUseFlags) {
              ROArrayColumn<bool> flagCol(itsIter.table(), "FLAG");
              if (itsUseAllChan) {
                itsBuffer.setFlags (flagCol.getColumn());
              } else {
                itsBuffer.setFlags (flagCol.getColumn(itsColSlicer));
              }
              // Set flags if FLAG_ROW is set.
              ROScalarColumn<bool> flagrowCol(itsIter.table(), "FLAG_ROW");
              for (uint i=0; i<itsIter.table().nrow(); ++i) {
                if (flagrowCol(i)) {
                  itsBuffer.getFlags()
                    (IPosition(3,0,0,i),
                     IPosition(3,itsNrCorr-1,itsNrChan-1,i)) = true;
                }
              }
            } else {
              // Do not use FLAG from the MS.
              itsBuffer.getFlags().resize (itsNrCorr, itsNrChan, itsNrBl);
              itsBuffer.getFlags() = false;
            }
            // Flag invalid data (NaN, infinite).
            const Complex* dataPtr = itsBuffer.getData().data();
            bool* flagPtr = itsBuffer.getFlags().data();
            for (uint i=0; i<itsBuffer.getData().size();) {
              for (uint j=i; j<i+itsNrCorr; ++j) {
                bool flag = (!isFinite(dataPtr[j].real())  ||
                             !isFinite(dataPtr[j].imag()));
                if (flag) {
                  itsFlagCounter.incrCorrelation(j-i);
                }
                if (flag  ||  flagPtr[j]) {
                  // Flag all correlations if a single one is flagged.
                  for (uint k=i; k<i+itsNrCorr; ++k) {
                    flagPtr[k] = true;
                  }
                  break;
                }
              }
              i += itsNrCorr;
            }
          }
          itsLastMSTime = itsNextTime;
          itsNrRead++;
          itsIter.next();
        }
        ASSERTSTR (itsBuffer.getFlags().shape()[2] == int(itsNrBl),
                   "#baselines is not the same for all time slots in the MS");
      }   // end of scope stops the timer.
      // Let the next step in the pipeline process this time slot.
      getNextStep()->process (itsBuffer);
      ///      cout << "Reader: " << itsNextTime-4.75e9<<endl;
      // Do not add to previous time, because it introduces round-off errors.
      itsNextTime = itsFirstTime + (itsNrRead+itsNrInserted) * itsTimeInterval;
      return true;
    }

    void MSReader::finish()
    {
      getNextStep()->finish();
    }

    void MSReader::show (std::ostream& os) const
    {
      os << "MSReader" << std::endl;
      os << "  input MS:       " << itsMSName << std::endl;
      if (itsMS.isNull()) {
        os << "    *** MS does not exist ***" << std::endl;
      } else {
        if (! itsSelBL.empty()) {
          os << "  baseline:       " << itsSelBL << std::endl;
        }
        os << "  band            " << itsSpw << std::endl;
        os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr
           << ')' << std::endl;
        os << "  nchan:          " << getInfo().nchan() << "  (" << itsNrChanStr
           << ')' << std::endl;
        os << "  ncorrelations:  " << getInfo().ncorr() << std::endl;
        uint nrbl = getInfo().nbaselines();
        os << "  nbaselines:     " << nrbl << std::endl;
        os << "  ntimes:         " << itsMS.nrow() / nrbl << std::endl;
        os << "  time interval:  " << getInfo().timeInterval() << std::endl;
        os << "  DATA column:    " << itsDataColName;
        if (itsMissingData) {
          os << "  (not present)";
        }
        os << std::endl;
        os << "  autoweight:     " << itsAutoWeight << std::endl;
      }
    }

    void MSReader::showCounts (std::ostream& os) const
    {
      os << endl << "NaN/infinite data flagged in reader";
      os << endl << "===================================" << endl;
      int64 nrtim = itsNrRead;
      itsFlagCounter.showCorrelation (os, nrtim);
      os << itsNrInserted << " missing time slots were inserted" << endl;
    }

    void MSReader::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " MSReader" << endl;
    }

    void MSReader::prepare (double& firstTime, double& lastTime,
                            double& interval)
    {
      ASSERT (itsMS.nrow() > 0);
      // Test if WEIGHT_SPECTRUM is present.
      TableDesc tdesc = itsMS.tableDesc();
      itsHasWeightSpectrum = false;
      if (tdesc.isColumn("WEIGHT_SPECTRUM")) {
        // The column is there, but it might not contain values. Test row 0.
        itsHasWeightSpectrum =
          ROArrayColumn<float>(itsMS, "WEIGHT_SPECTRUM").isDefined(0);
      }
      // Test if the data column is present.
      if (tdesc.isColumn (itsDataColName)) {
        itsMissingData = false;
      } else {
        if (itsMissingData) {
          // Only give warning if a missing data column is allowed.
          LOG_WARN ("Data column " + itsDataColName +
                    " is missing in " + itsMSName);
        } else {
          THROW (Exception, ("Data column " + itsDataColName +
                             " is missing in " + itsMSName));
        }
      }
      // Test if the full resolution flags are present.
      itsHasFullResFlags = tdesc.isColumn("LOFAR_FULL_RES_FLAG");
      if (itsHasFullResFlags) {
        ROTableColumn fullResFlagCol(itsMS, "LOFAR_FULL_RES_FLAG");
        itsFullResNChanAvg = fullResFlagCol.keywordSet().asInt ("NCHAN_AVG");
        itsFullResNTimeAvg = fullResFlagCol.keywordSet().asInt ("NTIME_AVG");
      } else {
        itsFullResNChanAvg = 1;
        itsFullResNTimeAvg = 1;
      }
      // Get the main table in the correct order.
      // Determine if the data are stored using LofarStMan.
      // If so, we know it is in time order.
      // (sorting on TIME with LofarStMan can be expensive).
      bool needSort = itsNeedSort;
      bool useRaw   = false;
      Record dminfo = itsMS.dataManagerInfo();
      for (unsigned i=0; i<dminfo.nfields(); ++i) {
        Record subrec = dminfo.subRecord(i);
        if (subrec.asString("TYPE") == "LofarStMan") {
          needSort = false;
          useRaw   = true;
          break;
        }
      }
      // Give an error if autoweight is used for a non-raw MS.
      if (itsAutoWeightForce) {
        itsAutoWeight = true;
      } else if (!useRaw && itsAutoWeight) {
        THROW (Exception, "Using autoweight=true cannot be done on DPPP-ed MS");
      }
      // If not in order, sort the main table (also on baseline).
      Table sortms(itsMS);
      Block<String> sortCols(3);
      sortCols[0] = "TIME";
      sortCols[1] = "ANTENNA1";
      sortCols[2] = "ANTENNA2";
      if (needSort) {
        sortms = itsMS.sort(sortCols);
      }
      // Get first and last time and interval from MS.
      firstTime = ROScalarColumn<double>(sortms, "TIME")(0);
      lastTime  = ROScalarColumn<double>(sortms, "TIME")(sortms.nrow()-1);
      interval  = ROScalarColumn<double>(sortms, "INTERVAL")(0);
      // Create iterator over time. Do not sort again.
      itsIter = TableIterator (sortms, Block<String>(1, "TIME"),
                               TableIterator::Ascending,
                               TableIterator::NoSort);
      // Find the nr of corr, chan, and baseline.
      IPosition shp (ROArrayColumn<Complex>(itsMS, "DATA").shape(0));
      itsNrCorr = shp[0];
      itsNrChan = shp[1];
      itsNrBl   = itsIter.table().nrow();
      // Ensure we have only one band by checking the nr of unique baselines.
      Table sortab = itsIter.table().sort(sortCols, Sort::Ascending,
                                          Sort::QuickSort + Sort::NoDuplicates);
      ASSERTSTR (sortab.nrow() == itsNrBl,
                 "The MS appears to have multiple subbands");
      // Get the baseline columns.
      ROScalarColumn<Int> ant1col(itsIter.table(), "ANTENNA1");
      ROScalarColumn<Int> ant2col(itsIter.table(), "ANTENNA2");
      // Keep the row numbers of the first part to be used for the meta info
      // of possibly missing time slots.
      itsBaseRowNrs = itsIter.table().rowNumbers(itsMS);
      // Get the antenna names and positions.
      Table anttab(itsMS.keywordSet().asTable("ANTENNA"));
      ROScalarColumn<String> nameCol (anttab, "NAME");
      ROScalarColumn<Double> diamCol (anttab, "DISH_DIAMETER");
      uint nant = anttab.nrow();
      ROScalarMeasColumn<MPosition> antcol (anttab, "POSITION");
      vector<MPosition> antPos;
      antPos.reserve (nant);
      for (uint i=0; i<nant; ++i) {
        antPos.push_back (antcol(i));
      }
      // Set antenna/baseline info.
      info().set (nameCol.getColumn(), diamCol.getColumn(), antPos,
                  ant1col.getColumn(), ant2col.getColumn());
      // Read the phase reference position from the FIELD subtable.
      // Only use the main value from the PHASE_DIR array.
      // The same for DELAY_DIR and LOFAR_TILE_BEAM_DIR.
      // If LOFAR_TILE_BEAM_DIR does not exist, use DELAY_DIR.
      Table fldtab (itsMS.keywordSet().asTable ("FIELD"));
      AlwaysAssert (fldtab.nrow() == 1, AipsError);
      MDirection phaseCenter, delayCenter, tileBeamDir;
      ROArrayMeasColumn<MDirection> fldcol1 (fldtab, "PHASE_DIR");
      ROArrayMeasColumn<MDirection> fldcol2 (fldtab, "DELAY_DIR");
      phaseCenter = *(fldcol1(0).data());
      delayCenter = *(fldcol2(0).data());
      if (fldtab.tableDesc().isColumn ("LOFAR_TILE_BEAM_DIR")) {
        ROArrayMeasColumn<MDirection> fldcol3 (fldtab, "LOFAR_TILE_BEAM_DIR");
        tileBeamDir = *(fldcol3(0).data());
      } else {
        tileBeamDir = delayCenter;
      }
      // Get the array position using the telescope name from the OBSERVATION
      // subtable. 
      Table obstab (itsMS.keywordSet().asTable ("OBSERVATION"));
      ROScalarColumn<String> telCol(obstab, "TELESCOPE_NAME");
      MPosition arrayPos;
      if (obstab.nrow() == 0  ||
          ! MeasTable::Observatory(arrayPos, telCol(0))) {
        // If not found, use the position of the middle antenna.
        arrayPos = antPos[antPos.size() / 2];
      }
      info().set (arrayPos, phaseCenter, delayCenter, tileBeamDir);
      // Create the UVW calculator.
      itsUVWCalc = UVWCalculator (phaseCenter, arrayPos, antPos);
    }

    void MSReader::prepare2()
    {
      // Set the info.
      uint ntime = uint((itsLastTime - itsFirstTime)/itsTimeInterval + 1.5);
      info().init (itsNrCorr, itsNrChan, ntime, itsStartTime,
                   itsTimeInterval, itsMSName);
      // Read the center frequencies of all channels.
      Table spwtab(itsMS.keywordSet().asTable("SPECTRAL_WINDOW"));
      ROArrayColumn<double> freqCol  (spwtab, "CHAN_FREQ");
      ROArrayColumn<double> widthCol (spwtab, "CHAN_WIDTH");
      ROArrayColumn<double> resolCol (spwtab, "RESOLUTION");
      ROArrayColumn<double> effBWCol (spwtab, "EFFECTIVE_BW");
      ROScalarColumn<Double> refCol  (spwtab, "REF_FREQUENCY");
      Vector<double> chanFreqs   = freqCol(itsSpw);
      Vector<double> chanWidths  = widthCol(itsSpw);
      Vector<double> resolutions = resolCol(itsSpw);
      Vector<double> effectiveBW = effBWCol(itsSpw);
      double refFreq = refCol(itsSpw);
      if (itsUseAllChan) {
        info().set (chanFreqs, chanWidths, resolutions, effectiveBW,
                    sum(effectiveBW), refFreq);
      } else {
        Vector<double> cwSlice(effectiveBW(Slice(itsStartChan, itsNrChan)));
        info().set (chanFreqs  (Slice(itsStartChan, itsNrChan)),
                    chanWidths (Slice(itsStartChan, itsNrChan)),
                    resolutions(Slice(itsStartChan, itsNrChan)),
                    cwSlice,
                    sum(cwSlice), refFreq);
      }
    }

    void MSReader::skipFirstTimes()
    {
      while (!itsIter.pastEnd()) {
        // Take time from row 0 in subset.
        double mstime = ROScalarColumn<double>(itsIter.table(), "TIME")(0);
        // Skip time slot and give warning if MS data is not in time order.
        if (mstime < itsLastMSTime) {
          LOG_WARN_STR ("Time at rownr "
                        << itsIter.table().rowNumbers(itsMS)[0]
                        << " of MS " << itsMS.tableName()
                        << " is less than previous time slot");
        } else {
          // Stop skipping if time equal to itsFirstTime.
          if (near(mstime, itsFirstTime)) {
            itsFirstTime = mstime;
            break;
          }
          // Also stop if time > itsFirstTime.
          // In that case determine the true first time, because itsFirstTime
          // can be a time value that does not coincide with a true time.
          // Note that a time stamp might be missing at this point,
          // so do not simply assume that mstime can be used.
          if (mstime > itsFirstTime) {
            int nrt = int((mstime - itsFirstTime) / itsTimeInterval);
            mstime -= (nrt+1) * itsTimeInterval;   // Add 1 for rounding errors
            if (near(mstime, itsFirstTime)) {
              itsFirstTime = mstime;
            } else {
              itsFirstTime = mstime + itsTimeInterval;
            }
            break;
          }
        }
        // Skip this time slot.
        itsLastMSTime = mstime;
        itsIter.next();
      }
    }

    void MSReader::calcUVW()
    {
      Matrix<double> uvws(3, itsNrBl);
      const Vector<Int>& ant1 = getInfo().getAnt1();
      const Vector<Int>& ant2 = getInfo().getAnt2();
      for (uint i=0; i<itsNrBl; ++i) {
        uvws.column(i) = itsUVWCalc.getUVW (ant1[i], ant2[i], itsNextTime);
      }
      itsBuffer.setUVW (uvws);
    }

    Matrix<double> MSReader::getUVW (const RefRows& rowNrs)
    {
      NSTimer::StartStop sstime(itsTimer);
      // Empty rownrs cannot happen for data, because in that case the buffer
      // should contain UVW for a missing time slot.
      ASSERT (! rowNrs.rowVector().empty());
      ROArrayColumn<double> dataCol(itsMS, "UVW");
      return dataCol.getColumnCells (rowNrs);
    }

    Cube<float> MSReader::getWeights (const RefRows& rowNrs,
                                      const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      if (rowNrs.rowVector().empty()) {
        // rowNrs can be empty if a time slot was inserted.
        Cube<float> weights(itsNrCorr, itsNrChan, itsNrBl);
        weights = 0;
        return weights;
      }
      Cube<float> weights;
      // Get weights for entire spectrum if present.
      if (itsHasWeightSpectrum) {
        ROArrayColumn<float> wsCol(itsMS, "WEIGHT_SPECTRUM");
        // Using getColumnCells(rowNrs,itsColSlicer) fails for LofarStMan.
        // Hence work around it.
        weights.reference (wsCol.getColumnCells (rowNrs));
        if (!itsUseAllChan) {
          // Make a copy, so the weights are consecutive in memory.
          weights.reference (weights(itsArrSlicer).copy());
        }
      } else {
        // No spectrum present; get global weights and assign to each channel.
        ROArrayColumn<float> wCol(itsMS, "WEIGHT");
        Matrix<float> inArr = wCol.getColumnCells (rowNrs);
        Cube<float> outArr(itsNrCorr, itsNrChan, itsNrBl);
        float* inPtr  = inArr.data();
        float* outPtr = outArr.data();
        for (uint i=0; i<itsNrBl; ++i) {
          // If global weights are zero, set them to 1. Some old MSs need that.
          for (uint k=0; k<itsNrCorr; ++k) {
            if (inPtr[k] == 0.) {
              inPtr[k] = 1.;
            }
          }
          for (uint j=0; j<itsNrChan; ++j) {
            for (uint k=0; k<itsNrCorr; ++k) {
              *outPtr++ = inPtr[k];
            }
          }
          inPtr += itsNrCorr;
        }
        weights.reference (outArr);
      }
      if (itsAutoWeight) {
        // Adapt weights using autocorrelations.
        autoWeight (weights, buf);
      }
      return weights;
    }

    void MSReader::autoWeight (Cube<float>& weights, const DPBuffer& buf)
    {
      const double* chanWidths = getInfo().chanWidths().data();
      uint npol  = weights.shape()[0];
      uint nchan = weights.shape()[1];
      uint nbl   = weights.shape()[2];
      // Get the autocorrelations indices.
      const vector<int>& autoInx = getInfo().getAutoCorrIndex();
      // Calculate the weight for each cross-correlation data point.
      const Vector<Int>& ant1 = getInfo().getAnt1();
      const Vector<Int>& ant2 = getInfo().getAnt2();
      const Complex* data = buf.getData().data();
      float* weight = weights.data();
      for (uint i=0; i<nbl; ++i) {
        // Can only be done if both autocorrelations are present.
        if (ant1[i] != ant2[i]  &&
            autoInx[ant1[i]] >= 0  &&  autoInx[ant2[i]] >= 0) {
          // Get offset of both autocorrelations in data array.
          const Complex* auto1 = data + autoInx[ant1[i]]*nchan*npol;
          const Complex* auto2 = data + autoInx[ant2[i]]*nchan*npol;
          for (uint j=0; j<nchan; ++j) {
            if (auto1[0].real() != 0  &&  auto2[0].real() != 0) {
              double w = chanWidths[j] * itsTimeInterval;
              weight[0] *= w / (auto1[0].real() * auto2[0].real());      // XX
              if (npol == 4) {
                if (auto1[3].real() != 0  &&  auto2[3].real() != 0) {
                  weight[1] *= w / (auto1[0].real() * auto2[3].real());  // XY
                  weight[2] *= w / (auto1[3].real() * auto2[0].real());  // YX
                  weight[3] *= w / (auto1[3].real() * auto2[3].real());  // YY
                }
              } else if (npol == 2) {
                if (auto1[1].real() != 0  &&  auto2[1].real() != 0) {
                  weight[1] *= w / (auto1[1].real() * auto2[1].real());  // YY
                }
              }
            }
            // Set pointers to next channel.
            weight += npol;
            auto1  += npol;
            auto2  += npol;
          }
        } else {
          // No autocorrelations for this baseline, so skip it.
          weight += nchan*npol;
        }
      }
    }

    Cube<bool> MSReader::getFullResFlags (const RefRows& rowNrs)
    {
      NSTimer::StartStop sstime(itsTimer);
      // Return empty array if no fullRes flags.
      if (!itsHasFullResFlags  ||  rowNrs.rowVector().empty()) {
        return Cube<bool>();
      }
      ROArrayColumn<uChar> fullResFlagCol(itsMS, "LOFAR_FULL_RES_FLAG");
      int norigchan = itsNrChan * itsFullResNChanAvg;
      int origstart = itsStartChan * itsFullResNChanAvg;
      Array<uChar> chars = fullResFlagCol.getColumnCells (rowNrs);
      // The original flags are kept per channel, not per corr.
      // Per row the flags are stored as uchar[nchar,navgtime].
      // Each char contains a bit per channel, thus nchan/8 chars are needed.
      // ntimeavg is the nr of times used when averaging.
      // Return it as Cube<bool>[norigchan,ntimeavg,nrbl].
      IPosition chShape = chars.shape();
      IPosition ofShape(3, norigchan, chShape[1], chShape[2]);
      Cube<bool> flags(ofShape);
      // Now expand the bits to bools.
      // If all bits to convert are contiguous, do it all in one go.
      // Otherwise we have to iterate.
      if (ofShape[0] == chShape[0]*8) {
        Conversion::bitToBool (flags.data(), chars.data(), flags.size());
      } else {
        ASSERT (ofShape[0] < chShape[0]*8);
        const uChar* charsPtr = chars.data();
        bool* flagsPtr = flags.data();
        for (int i=0; i<ofShape[1]*ofShape[2]; ++i) {
          Conversion::bitToBool (flagsPtr, charsPtr, origstart, ofShape[0]);
          flagsPtr += ofShape[0];
          charsPtr += chShape[0];
        }
      }
      return flags;
    }

    Cube<Complex> MSReader::getData (const String& columnName,
                                     const RefRows& rowNrs)
    {
      NSTimer::StartStop sstime(itsTimer);
      // Empty rownrs cannot happen for data, because in that case the buffer
      // should contain data for a missing time slot.
      ASSERT (! rowNrs.rowVector().empty());
      ROArrayColumn<Complex> dataCol(itsMS, columnName);
      // Also work around LofarStMan/getColumnCells slice problem.
      Cube<Complex> data = dataCol.getColumnCells (rowNrs);
      return (itsUseAllChan ? data : data(itsArrSlicer));
    }

    void MSReader::putFlags (const RefRows& rowNrs,
                             const Cube<bool>& flags)
    {
      if (! rowNrs.rowVector().empty()) {
        itsMS.reopenRW();
        ArrayColumn<bool> flagCol(itsMS, "FLAG");
        ScalarColumn<bool> flagRowCol(itsMS, "FLAG_ROW");
        // Loop over all rows of this subset.
	// (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<bool> flagIter (flags, 2);
        for (uint i=0; i<rows.size(); ++i) {
          flagCol.putSlice (rows[i], itsColSlicer, flagIter.array());
          // If a new flag in a row is clear, the ROW_FLAG should not be set.
          // If all new flags are set, we leave it because we might have a
          // subset of the channels, so other flags might still be clear.
          if (anyEQ (flagIter.array(), False)) {
            flagRowCol.put (i, False);
          }
          flagIter.next();
	}
      }
    }

  } //# end namespace
}
