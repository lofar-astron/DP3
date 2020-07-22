// MSReader.cc: DPPP step reading from an MS
// Copyright (C) 2010
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$
//
// @author Ger van Diepen

#include "MSReader.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"

#include <EveryBeam/load.h>
#include <EveryBeam/LofarMetaDataUtil.h>

#include "../Common/ParameterSet.h"
#include "../Common/BaselineSelect.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>
#include <casacore/casa/version.h>

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/OS/Conversion.h>

#include <iostream>

using namespace casacore;

namespace DP3 {
namespace DPPP {

MSReader::MSReader()
    : itsReadVisData(False), itsLastMSTime(0), itsNrRead(0), itsNrInserted(0) {}

MSReader::MSReader(const string& msName, const ParameterSet& parset,
                   const string& prefix, bool missingData)
    : itsReadVisData(False),
      itsMissingData(missingData),
      itsLastMSTime(0),
      itsNrRead(0),
      itsNrInserted(0) {
  NSTimer::StartStop sstime(itsTimer);
  // Get info from parset.
  itsSpw = parset.getInt(prefix + "band", -1);
  itsStartChanStr = parset.getString(prefix + "startchan", "0");
  itsNrChanStr = parset.getString(prefix + "nchan", "0");
  string startTimeStr = parset.getString(prefix + "starttime", "");
  string endTimeStr = parset.getString(prefix + "endtime", "");
  unsigned int nTimes = parset.getInt(prefix + "ntimes", 0);
  int startTimeSlot = parset.getInt(prefix + "starttimeslot", 0);
  // Can be negative to insert flagged time slots before start
  itsTimeTolerance = parset.getDouble(prefix + "timetolerance", 1e-2);
  itsUseFlags = parset.getBool(prefix + "useflag", true);
  itsDataColName = parset.getString(prefix + "datacolumn", "DATA");
  itsWeightColName =
      parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM");
  itsModelColName = parset.getString(prefix + "modelcolumn", "MODEL_DATA");
  itsAutoWeight = parset.getBool(prefix + "autoweight", false);
  itsAutoWeightForce = parset.getBool(prefix + "forceautoweight", false);
  itsNeedSort = parset.getBool(prefix + "sort", false);
  itsSelBL = parset.getString(prefix + "baseline", string());
  // Try to open the MS and get its full name.
  if (itsMissingData && !Table::isReadable(msName)) {
    DPLOG_WARN_STR("MeasurementSet " << msName
                                     << " not found; dummy data used");
    return;
  }
  itsMS = MeasurementSet(msName, TableLock::AutoNoReadLocking);
  itsSelMS = itsMS;
  itsMSName = itsMS.tableName();
  // See if a selection on band needs to be done.
  // We assume that DATA_DESC_ID and SPW_ID map 1-1.
  if (itsSpw >= 0) {
    DPLOG_INFO_STR(" MSReader selecting spectral window " +
                   std::to_string(itsSpw) + " ...");
    Table subset = itsSelMS(itsSelMS.col("DATA_DESC_ID") == itsSpw);
    // If not all is selected, use the selection.
    if (subset.nrow() < itsSelMS.nrow()) {
      if (subset.nrow() <= 0)
        throw Exception("Band " + std::to_string(itsSpw) + " not found in " +
                        itsMSName);
      itsSelMS = subset;
    }
  } else {
    itsSpw = 0;
  }
  // See if a selection on baseline needs to be done.
  if (!itsSelBL.empty()) {
    DPLOG_INFO_STR(" MSReader selecting baselines ...");
    MSSelection select;

    // Overwrite the error handler to ignore errors for unknown antennas.
    // Borrowed from BaselineSelect.cc
    std::ostringstream os;
    auto curHandler = MSAntennaParse::thisMSAErrorHandler;
#if CASACORE_MAJOR_VERSION < 3 ||    \
    (CASACORE_MAJOR_VERSION == 3 &&  \
     (CASACORE_MINOR_VERSION == 0 || \
      (CASACORE_MINOR_VERSION == 1 && CASACORE_PATCH_VERSION < 2)))
    // In casacore < 3.1.2 thisMSAErrorHandler is a raw pointer,
    // From casacore 3.1.2. it's a CountedPtr
    BaselineSelectErrorHandler errorHandler(os);
    MSAntennaParse::thisMSAErrorHandler = &errorHandler;
#else
    CountedPtr<MSSelectionErrorHandler> errorHandler(
        new BaselineSelectErrorHandler(os));
    MSAntennaParse::thisMSAErrorHandler = errorHandler;
#endif

    // Set given selection strings.
    try {
      select.setAntennaExpr(itsSelBL);
      // Create a table expression for an MS representing the selection.
      MeasurementSet ms(itsSelMS);
      TableExprNode node = select.toTableExprNode(&ms);
      Table subset = itsSelMS(node);
      // If not all is selected, use the selection.
      if (subset.nrow() < itsSelMS.nrow()) {
        if (subset.nrow() <= 0)
          throw Exception("Baselines " + itsSelBL + "not found in " +
                          itsMSName);
        itsSelMS = subset;
      }
      MSAntennaParse::thisMSAErrorHandler = curHandler;
    } catch (const std::exception&) {
      MSAntennaParse::thisMSAErrorHandler = curHandler;
      throw;
    }
  }
  // Prepare the MS access and get time info.
  double startTimeMS = 0., endTimeMS = 0.;
  prepare(startTimeMS, endTimeMS, itsTimeInterval);

  // Start and end time can be given in the parset in case leading
  // or trailing time slots are missing.
  // They can also be used to select part of the MS.

  if (!startTimeStr.empty()) {
    if (startTimeSlot > 0) {
      throw std::runtime_error("Only one of " + prefix + "starttimeslot and " +
                               prefix + "starttime can be specified");
    }
    Quantity qtime;
    if (!MVTime::read(qtime, startTimeStr)) {
      throw std::runtime_error(startTimeStr + " is an invalid date/time");
    }
    double startTimeParset = qtime.getValue("s");
    // the parset specified start time is allowed to be before the msstarttime.
    // In that case, flagged samples are injected.
    if (startTimeParset > endTimeMS)
      throw std::runtime_error("Specified starttime is past end of time axis");

    // Round specified first time to a multiple of itsTimeInterval
    itsFirstTime = startTimeMS + std::ceil((startTimeParset - startTimeMS) /
                                           itsTimeInterval) *
                                     itsTimeInterval;
  } else {
    itsFirstTime = startTimeMS + startTimeSlot * itsTimeInterval;
  }

  if (!endTimeStr.empty()) {
    Quantity qtime;
    if (!MVTime::read(qtime, endTimeStr)) {
      throw std::runtime_error(endTimeStr + " is an invalid date/time");
    }
    double endTimeParset = qtime.getValue("s");
    // Some overlap between the measurement set timerange and the parset range
    // is required :
    if (endTimeParset < startTimeMS + 0.5 * itsTimeInterval) {
      throw std::runtime_error(
          "Specified end time " + endTimeStr +
          " is before the first timestep in the measurement set");
    }
    // Round specified first time to a multiple of itsTimeInterval
    itsLastTime = startTimeMS +
                  std::floor((endTimeParset - startTimeMS) / itsTimeInterval) *
                      itsTimeInterval;
  } else {
    itsLastTime = endTimeMS;
  }

  if (itsLastTime < itsFirstTime)
    throw std::runtime_error("Specified endtime is before specified starttime");
  // If needed, skip the first times in the MS.
  // It also sets itsFirstTime properly (round to time/interval in MS).
  skipFirstTimes();
  if (nTimes > 0) {
    if (!endTimeStr.empty()) {
      throw std::runtime_error("Only one of " + prefix + "ntimes and " +
                               prefix + "endtime can be specified");
    }
    itsLastTime = itsFirstTime + (nTimes - 1) * itsTimeInterval;
  }
  itsNextTime = itsFirstTime;
  itsStartTime = itsFirstTime - 0.5 * itsTimeInterval;
  // Parse the chan expressions.
  // Nr of channels can be used as 'nchan' in the expressions.
  Record rec;
  rec.define("nchan", itsNrChan);
  TableExprNode node1(RecordGram::parse(rec, itsStartChanStr));
  TableExprNode node2(RecordGram::parse(rec, itsNrChanStr));
  // nchan=0 means until the last channel.
  double result;
  node1.get(rec, result);
  itsStartChan = (unsigned int)(result + 0.001);
  node2.get(rec, result);
  unsigned int nrChan = (unsigned int)(result + 0.0001);
  unsigned int nAllChan = itsNrChan;
  if (itsStartChan >= nAllChan)
    throw Exception("startchan " + std::to_string(itsStartChan) +
                    " exceeds nr of channels in MS (" +
                    std::to_string(nAllChan) + ')');
  unsigned int maxNrChan = nAllChan - itsStartChan;
  if (nrChan == 0) {
    itsNrChan = maxNrChan;
  } else {
    itsNrChan = std::min(nrChan, maxNrChan);
  }
  // Are all channels used?
  itsUseAllChan = itsStartChan == 0 && itsNrChan == nAllChan;
  // Do the rest of the preparation.
  prepare2();
  // Take subset of channel frequencies if needed.
  // Make sure to copy the subset to get a proper Vector.
  // Form the slicer to get channels and correlations from column.
  itsColSlicer =
      Slicer(IPosition(2, 0, itsStartChan), IPosition(2, itsNrCorr, itsNrChan));
  // Form the slicer to get channels, corrs, and baselines from array.
  itsArrSlicer = Slicer(IPosition(3, 0, itsStartChan, 0),
                        IPosition(3, itsNrCorr, itsNrChan, itsNrBl));
  // Initialize the flag counters.
  itsFlagCounter.init(getInfo());
}

MSReader::~MSReader() {}

void MSReader::updateInfo(const DPInfo& dpInfo) {
  info().setNThreads(dpInfo.nThreads());
}

std::string MSReader::msName() const { return itsMSName; }

void MSReader::setReadVisData(bool readVisData) {
  itsReadVisData = (itsReadVisData || readVisData);
}

bool MSReader::process(const DPBuffer&) {
  if (itsNrRead == 0) {
    if (itsReadVisData) {
      itsBuffer.getData().resize(itsNrCorr, itsNrChan, itsNrBl);
    }
    if (itsUseFlags) {
      itsBuffer.getFlags().resize(itsNrCorr, itsNrChan, itsNrBl);
    }
    /// cout<<(void*)(itsBuffer.getData().data())<<" upd"<<endl;
  }
  {
    NSTimer::StartStop sstime(itsTimer);
    ///        itsBuffer.clear();
    // Use time from the current time slot in the MS.
    bool useIter = false;
    while (!itsIter.pastEnd()) {
      // Take time from row 0 in subset.
      double mstime = ScalarColumn<double>(itsIter.table(), "TIME")(0);
      // Skip time slot and give warning if MS data is not in time order.
      if (mstime < itsLastMSTime) {
        DPLOG_WARN_STR("Time at rownr " +
                       std::to_string(itsIter.table().rowNumbers(itsMS)[0]) +
                       " of MS " + itsMSName +
                       " is less than previous time slot");
      } else {
        // Use the time slot if near or < nexttime, but > starttime.
        // In this way we cater for irregular times in some WSRT MSs.
        if (nearAbs(mstime, itsNextTime, itsTimeTolerance)) {
          useIter = true;
          break;
        } else if (mstime > itsFirstTime && mstime < itsNextTime) {
          itsFirstTime -= itsNextTime - mstime;
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
    // Stop if at the end, or if there is no data at all
    if ((itsNextTime > itsLastTime && !near(itsNextTime, itsLastTime)) ||
        itsNextTime == 0.) {
      return false;
    }
    // Fill the buffer.
    itsBuffer.setTime(itsNextTime);
    /// cout << "read time " <<itsBuffer.getTime() - 4472025855.0<<endl;
    if (!useIter) {
      // Need to insert a fully flagged time slot.
      itsBuffer.setRowNrs(Vector<rownr_t>());
      itsBuffer.setExposure(itsTimeInterval);
      itsBuffer.getFlags() = true;
      if (itsReadVisData) {
        itsBuffer.getData() = Complex();
      }
      itsNrInserted++;
    } else {
      itsBuffer.setRowNrs(itsIter.table().rowNumbers(itsMS, True));
      if (itsMissingData) {
        // Data column not present, so fill a fully flagged time slot.
        itsBuffer.setExposure(itsTimeInterval);
        itsBuffer.getFlags() = true;
        if (itsReadVisData) {
          itsBuffer.getData() = Complex();
        }
      } else {
        // Set exposure.
        itsBuffer.setExposure(
            ScalarColumn<double>(itsIter.table(), "EXPOSURE")(0));
        // Get data and flags from the MS.
        ///            if (itsNrRead%50 < 4) {
        ///              cout<<(void*)(itsBuffer.getData().data())<<"
        ///              rd1"<<endl;
        ///}
        if (itsReadVisData) {
          ArrayColumn<Complex> dataCol(itsIter.table(), itsDataColName);
          if (itsUseAllChan) {
            dataCol.getColumn(itsBuffer.getData());
          } else {
            dataCol.getColumn(itsColSlicer, itsBuffer.getData());
          }
        }
        /// if (itsNrRead%50 < 4) {
        /// cout<<(void*)(itsBuffer.getData().data())<<" rd2"<<endl;
        ///}
        if (itsUseFlags) {
          ArrayColumn<bool> flagCol(itsIter.table(), "FLAG");
          if (itsUseAllChan) {
            flagCol.getColumn(itsBuffer.getFlags());
          } else {
            flagCol.getColumn(itsColSlicer, itsBuffer.getFlags());
          }
          // Set flags if FLAG_ROW is set.
          ScalarColumn<bool> flagrowCol(itsIter.table(), "FLAG_ROW");
          for (unsigned int i = 0; i < itsIter.table().nrow(); ++i) {
            if (flagrowCol(i)) {
              itsBuffer.getFlags()(
                  IPosition(3, 0, 0, i),
                  IPosition(3, itsNrCorr - 1, itsNrChan - 1, i)) = true;
            }
          }
        } else {
          // Do not use FLAG from the MS.
          itsBuffer.getFlags().resize(itsNrCorr, itsNrChan, itsNrBl);
          itsBuffer.getFlags() = false;
        }
        // Flag invalid data (NaN, infinite).
        flagInfNaN(itsBuffer.getData(), itsBuffer.getFlags(), itsFlagCounter);
      }
      itsLastMSTime = itsNextTime;
      itsNrRead++;
      itsIter.next();
    }
    if (itsBuffer.getFlags().shape()[2] != int(itsNrBl))
      throw Exception(
          "#baselines is not the same for all time slots in the MS");
  }  // end of scope stops the timer.
  // Let the next step in the pipeline process this time slot.
  getNextStep()->process(itsBuffer);
  ///      cout << "Reader: " << itsNextTime-4.75e9<<endl;
  // Do not add to previous time, because it introduces round-off errors.
  itsNextTime = itsFirstTime + (itsNrRead + itsNrInserted) * itsTimeInterval;
  return true;
}

void MSReader::flagInfNaN(const casacore::Cube<casacore::Complex>& dataCube,
                          casacore::Cube<bool>& flagsCube,
                          FlagCounter& flagCounter) {
  int ncorr = dataCube.shape()[0];
  const Complex* dataPtr = dataCube.data();
  bool* flagPtr = flagsCube.data();
  for (unsigned int i = 0; i < dataCube.size();) {
    for (unsigned int j = i; j < i + ncorr; ++j) {
      bool flag =
          (!isFinite(dataPtr[j].real()) || !isFinite(dataPtr[j].imag()));
      if (flag) {
        flagCounter.incrCorrelation(j - i);
      }
      if (flag || flagPtr[j]) {
        // Flag all correlations if a single one is flagged.
        for (unsigned int k = i; k < i + ncorr; ++k) {
          flagPtr[k] = true;
        }
        break;
      }
    }
    i += ncorr;
  }
}

void MSReader::finish() { getNextStep()->finish(); }

void MSReader::show(std::ostream& os) const {
  os << "MSReader\n";
  os << "  input MS:       " << itsMSName << '\n';
  if (itsMS.isNull()) {
    os << "    *** MS does not exist ***\n";
  } else {
    if (!itsSelBL.empty()) {
      os << "  baseline:       " << itsSelBL << '\n';
    }
    os << "  band            " << itsSpw << '\n';
    os << "  startchan:      " << itsStartChan << "  (" << itsStartChanStr
       << ")\n";
    os << "  nchan:          " << getInfo().nchan() << "  (" << itsNrChanStr
       << ")\n";
    os << "  ncorrelations:  " << getInfo().ncorr() << '\n';
    unsigned int nrbl = getInfo().nbaselines();
    os << "  nbaselines:     " << nrbl << '\n';
    os << "  first time:     " << MVTime::Format(MVTime::YMD)
       << MVTime(itsFirstTime / (24 * 3600.)) << '\n';
    os << "  last time:      " << MVTime::Format(MVTime::YMD)
       << MVTime(itsLastTime / (24 * 3600.)) << '\n';
    os << "  ntimes:         " << getInfo().ntime()
       << '\n';  // itsSelMS can contain timeslots that are ignored in process
    os << "  time interval:  " << getInfo().timeInterval() << '\n';
    os << "  DATA column:    " << itsDataColName;
    if (itsMissingData) {
      os << "  (not present)";
    }
    os << '\n';
    os << "  WEIGHT column:  " << itsWeightColName << '\n';
    os << "  autoweight:     " << boolalpha << itsAutoWeight << '\n';
  }
}

void MSReader::showCounts(std::ostream& os) const {
  os << endl << "NaN/infinite data flagged in reader";
  os << endl << "===================================" << endl;
  int64_t nrtim = itsNrRead;
  itsFlagCounter.showCorrelation(os, nrtim);
  os << itsNrInserted << " missing time slots were inserted" << endl;
}

void MSReader::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " MSReader" << endl;
}

void MSReader::prepare(double& firstTime, double& lastTime, double& interval) {
  if (itsSelMS.nrow() == 0) {
    DPLOG_WARN_STR("The selected input does not contain any data.");
  }
  TableDesc tdesc = itsMS.tableDesc();

  itsHasWeightSpectrum = false;
  // if weightcolname is specified to "WEIGHT" then this is used, even
  // if a weight_spectrum is present.
  if (itsWeightColName != "WEIGHT") {
    // Test if specified weight column or WEIGHT_SPECTRUM is present.
    if (tdesc.isColumn(itsWeightColName)) {
      // The column is there, but it might not contain values. Test row 0.
      itsHasWeightSpectrum =
          ArrayColumn<float>(itsSelMS, itsWeightColName).isDefined(0);
      if (!itsHasWeightSpectrum && itsWeightColName != "WEIGHT_SPECTRUM") {
        DPLOG_WARN_STR("Specified weight column " + itsWeightColName +
                       "is not a valid column, using WEIGHT instead");
      }
    }
  }

  // Test if the data column is present.
  if (tdesc.isColumn(itsDataColName)) {
    itsMissingData = false;

    // Read beam keywords of input datacolumn
    ArrayColumn<Complex> dataCol(itsMS, itsDataColName);
    if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
      std::string mode =
          dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE");
      if (mode == "None")
        info().setBeamCorrectionMode(NoBeamCorrection);
      else {
        if (mode == "Element")
          info().setBeamCorrectionMode(ElementBeamCorrection);
        else if (mode == "ArrayFactor")
          info().setBeamCorrectionMode(ArrayFactorBeamCorrection);
        else if (mode == "Full")
          info().setBeamCorrectionMode(FullBeamCorrection);

        String error;
        MeasureHolder mHolder;
        if (!mHolder.fromRecord(
                error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
          throw std::runtime_error(error);
        info().setBeamCorrectionDir(mHolder.asMDirection());
      }
    }
  } else {
    if (itsMissingData) {
      // Only give warning if a missing data column is allowed.
      DPLOG_WARN_STR("Data column " + itsDataColName + " is missing in " +
                     itsMSName);
    } else {
      throw Exception("Data column " + itsDataColName + " is missing in " +
                      itsMSName);
    }
  }

  // Test if the full resolution flags are present.
  itsHasFullResFlags = tdesc.isColumn("LOFAR_FULL_RES_FLAG");
  if (itsHasFullResFlags) {
    ROTableColumn fullResFlagCol(itsMS, "LOFAR_FULL_RES_FLAG");
    itsFullResNChanAvg = fullResFlagCol.keywordSet().asInt("NCHAN_AVG");
    itsFullResNTimeAvg = fullResFlagCol.keywordSet().asInt("NTIME_AVG");
  } else {
    itsFullResNChanAvg = 1;
    itsFullResNTimeAvg = 1;
  }
  // Get the main table in the correct order.
  // Determine if the data are stored using LofarStMan.
  // If so, we know it is in time order.
  // (sorting on TIME with LofarStMan can be expensive).
  bool needSort = itsNeedSort;
  bool useRaw = false;
  Record dminfo = itsMS.dataManagerInfo();
  for (unsigned i = 0; i < dminfo.nfields(); ++i) {
    Record subrec = dminfo.subRecord(i);
    if (subrec.asString("TYPE") == "LofarStMan") {
      needSort = false;
      useRaw = true;
      break;
    }
  }
  // Give an error if autoweight is used for a non-raw MS.
  if (itsAutoWeightForce) {
    itsAutoWeight = true;
  } else if (!useRaw && itsAutoWeight) {
    throw Exception("Using autoweight=true cannot be done on DPPP-ed MS");
  }
  // If not in order, sort the table selection (also on baseline).
  Table sortms(itsSelMS);
  Block<String> sortCols(3);
  sortCols[0] = "TIME";
  sortCols[1] = "ANTENNA1";
  sortCols[2] = "ANTENNA2";
  if (needSort) {
    sortms = itsSelMS.sort(sortCols);
  }
  // Get first and last time and interval from MS.
  if (itsSelMS.nrow() > 0) {
    firstTime = ScalarColumn<double>(sortms, "TIME")(0);
    lastTime = ScalarColumn<double>(sortms, "TIME")(sortms.nrow() - 1);
    interval = ScalarColumn<double>(sortms, "INTERVAL")(0);
  }
  // Create iterator over time. Do not sort again.
  itsIter = TableIterator(sortms, Block<String>(1, "TIME"),
                          TableIterator::Ascending, TableIterator::NoSort);
  // Find the nr of corr, chan, and baseline.
  IPosition shp(ArrayColumn<Complex>(itsSelMS, "DATA").shape(0));
  itsNrCorr = shp[0];
  itsNrChan = shp[1];
  itsNrBl = itsIter.table().nrow();
  // Ensure we have only one band by checking the nr of unique baselines.
  Table sortab = itsIter.table().sort(sortCols, Sort::Ascending,
                                      Sort::QuickSort + Sort::NoDuplicates);
  if (sortab.nrow() != itsNrBl)
    throw Exception("The MS appears to have multiple subbands");
  // Get the baseline columns.
  ScalarColumn<Int> ant1col(itsIter.table(), "ANTENNA1");
  ScalarColumn<Int> ant2col(itsIter.table(), "ANTENNA2");
  // Keep the row numbers of the first part to be used for the meta info
  // of possibly missing time slots.
  itsBaseRowNrs = itsIter.table().rowNumbers(itsMS, True);
  // Get the antenna names and positions.
  Table anttab(itsMS.keywordSet().asTable("ANTENNA"));
  ScalarColumn<String> nameCol(anttab, "NAME");
  ScalarColumn<Double> diamCol(anttab, "DISH_DIAMETER");
  unsigned int nant = anttab.nrow();
  ROScalarMeasColumn<MPosition> antcol(anttab, "POSITION");
  vector<MPosition> antPos;
  antPos.reserve(nant);
  for (unsigned int i = 0; i < nant; ++i) {
    antPos.push_back(antcol(i));
  }
  // Set antenna/baseline info.
  info().set(nameCol.getColumn(), diamCol.getColumn(), antPos,
             ant1col.getColumn(), ant2col.getColumn());

  if (itsAutoWeight) {
    info().setNeedVisData();
    info().setWriteWeights();
  }

  // Read the phase reference position from the FIELD subtable.
  // Only use the main value from the PHASE_DIR array.
  // The same for DELAY_DIR and LOFAR_TILE_BEAM_DIR.
  // If LOFAR_TILE_BEAM_DIR does not exist, use DELAY_DIR.
  Table fldtab(itsMS.keywordSet().asTable("FIELD"));
  if (fldtab.nrow() != 1)
    throw std::runtime_error("Multiple entries in FIELD table");
  MDirection phaseCenter, delayCenter, tileBeamDir;
  ROArrayMeasColumn<MDirection> fldcol1(fldtab, "PHASE_DIR");
  ROArrayMeasColumn<MDirection> fldcol2(fldtab, "DELAY_DIR");
  phaseCenter = *(fldcol1(0).data());
  delayCenter = *(fldcol2(0).data());

  tileBeamDir = everybeam::ReadTileBeamDirection(itsMS);

  // Get the array position using the telescope name from the OBSERVATION
  // subtable.
  Table obstab(itsMS.keywordSet().asTable("OBSERVATION"));
  ScalarColumn<String> telCol(obstab, "TELESCOPE_NAME");
  MPosition arrayPos;
  if (obstab.nrow() == 0 || !MeasTable::Observatory(arrayPos, telCol(0))) {
    // If not found, use the position of the middle antenna.
    arrayPos = antPos[antPos.size() / 2];
  }
  info().set(arrayPos, phaseCenter, delayCenter, tileBeamDir);
  // Create the UVW calculator.
  itsUVWCalc = UVWCalculator(phaseCenter, arrayPos, antPos);
}

void MSReader::prepare2() {
  // Set the info.
  // The 1.5 comes from a) rounding (0.5) + b) the paaltjesprobleem.
  unsigned int ntime =
      (unsigned int)((itsLastTime - itsFirstTime) / itsTimeInterval + 1.5);
  // Read the antenna set.
  Table obstab(itsMS.keywordSet().asTable("OBSERVATION"));
  string antennaSet;
  if (obstab.nrow() > 0 && obstab.tableDesc().isColumn("LOFAR_ANTENNA_SET")) {
    antennaSet = ScalarColumn<String>(obstab, "LOFAR_ANTENNA_SET")(0);
  }
  info().init(itsNrCorr, itsStartChan, itsNrChan, ntime, itsStartTime,
              itsTimeInterval, itsMSName, antennaSet);
  info().setDataColName(itsDataColName);
  info().setWeightColName(itsWeightColName);
  // Read the center frequencies of all channels.
  Table spwtab(itsMS.keywordSet().asTable("SPECTRAL_WINDOW"));
  ArrayColumn<double> freqCol(spwtab, "CHAN_FREQ");
  ArrayColumn<double> widthCol(spwtab, "CHAN_WIDTH");
  ArrayColumn<double> resolCol(spwtab, "RESOLUTION");
  ArrayColumn<double> effBWCol(spwtab, "EFFECTIVE_BW");
  ScalarColumn<Double> refCol(spwtab, "REF_FREQUENCY");
  Vector<double> chanFreqs = freqCol(itsSpw);
  Vector<double> chanWidths = widthCol(itsSpw);
  Vector<double> resolutions = resolCol(itsSpw);
  Vector<double> effectiveBW = effBWCol(itsSpw);
  double refFreq = refCol(itsSpw);
  if (itsUseAllChan) {
    info().set(chanFreqs, chanWidths, resolutions, effectiveBW,
               sum(effectiveBW), refFreq);
  } else {
    Vector<double> cwSlice(effectiveBW(Slice(itsStartChan, itsNrChan)));
    info().set(chanFreqs(Slice(itsStartChan, itsNrChan)),
               chanWidths(Slice(itsStartChan, itsNrChan)),
               resolutions(Slice(itsStartChan, itsNrChan)), cwSlice,
               sum(cwSlice), refFreq);
  }
}

void MSReader::skipFirstTimes() {
  while (!itsIter.pastEnd()) {
    // Take time from row 0 in subset.
    double mstime = ScalarColumn<double>(itsIter.table(), "TIME")(0);
    // Skip time slot and give warning if MS data is not in time order.
    if (mstime < itsLastMSTime) {
      DPLOG_WARN_STR("Time at rownr " +
                     std::to_string(itsIter.table().rowNumbers(itsMS)[0]) +
                     " of MS " + itsMSName +
                     " is less than previous time slot");
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
        mstime -= (nrt + 1) * itsTimeInterval;  // Add 1 for rounding errors
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

void MSReader::calcUVW(double time, DPBuffer& buf) {
  Matrix<double>& uvws = buf.getUVW();
  uvws.resize(3, itsNrBl);
  const Vector<Int>& ant1 = getInfo().getAnt1();
  const Vector<Int>& ant2 = getInfo().getAnt2();
  for (unsigned int i = 0; i < itsNrBl; ++i) {
    uvws.column(i) = itsUVWCalc.getUVW(ant1[i], ant2[i], time);
  }
}

void MSReader::getUVW(const RefRows& rowNrs, double time, DPBuffer& buf) {
  NSTimer::StartStop sstime(itsTimer);
  // Calculate UVWs if empty rownrs (i.e., missing data).
  if (rowNrs.rowVector().empty()) {
    calcUVW(time, buf);
  } else {
    ArrayColumn<double> dataCol(itsMS, "UVW");
    dataCol.getColumnCells(rowNrs, buf.getUVW());
  }
}

void MSReader::getWeights(const RefRows& rowNrs, DPBuffer& buf) {
  NSTimer::StartStop sstime(itsTimer);
  Cube<float>& weights = buf.getWeights();
  // Resize if needed (probably when called for first time).
  if (weights.empty()) {
    weights.resize(itsNrCorr, itsNrChan, itsNrBl);
  }
  if (rowNrs.rowVector().empty()) {
    // rowNrs can be empty if a time slot was inserted.
    weights = 0;
  } else {
    // Get weights for entire spectrum if present.
    if (itsHasWeightSpectrum) {
      ArrayColumn<float> wsCol(itsMS, itsWeightColName);
      // Using getColumnCells(rowNrs,itsColSlicer) fails for LofarStMan.
      // Hence work around it.
      if (itsUseAllChan) {
        wsCol.getColumnCells(rowNrs, weights);
      } else {
        Cube<float> w = wsCol.getColumnCells(rowNrs);
        weights = w(itsArrSlicer);
      }
    } else {
      // No spectrum present; get global weights and assign to each channel.
      ArrayColumn<float> wCol(itsMS, "WEIGHT");
      Matrix<float> inArr = wCol.getColumnCells(rowNrs);
      float* inPtr = inArr.data();
      float* outPtr = weights.data();
      for (unsigned int i = 0; i < itsNrBl; ++i) {
        // Set global weights to 1 if zero. Some old MSs need that.
        for (unsigned int k = 0; k < itsNrCorr; ++k) {
          if (inPtr[k] == 0.) {
            inPtr[k] = 1.;
          }
        }
        for (unsigned int j = 0; j < itsNrChan; ++j) {
          for (unsigned int k = 0; k < itsNrCorr; ++k) {
            *outPtr++ = inPtr[k];
          }
        }
        inPtr += itsNrCorr;
      }
    }
    if (itsAutoWeight) {
      // Adapt weights using autocorrelations.
      autoWeight(weights, buf);
    }
  }
}

void MSReader::autoWeight(Cube<float>& weights, const DPBuffer& buf) {
  const double* chanWidths = getInfo().chanWidths().data();
  unsigned int npol = weights.shape()[0];
  unsigned int nchan = weights.shape()[1];
  unsigned int nbl = weights.shape()[2];
  // Get the autocorrelations indices.
  const vector<int>& autoInx = getInfo().getAutoCorrIndex();
  // Calculate the weight for each cross-correlation data point.
  const Vector<Int>& ant1 = getInfo().getAnt1();
  const Vector<Int>& ant2 = getInfo().getAnt2();
  const Complex* data = buf.getData().data();
  float* weight = weights.data();
  for (unsigned int i = 0; i < nbl; ++i) {
    // Can only be done if both autocorrelations are present.
    if (autoInx[ant1[i]] >= 0 && autoInx[ant2[i]] >= 0) {
      // Get offset of both autocorrelations in data array.
      const Complex* auto1 = data + autoInx[ant1[i]] * nchan * npol;
      const Complex* auto2 = data + autoInx[ant2[i]] * nchan * npol;
      for (unsigned int j = 0; j < nchan; ++j) {
        if (auto1[0].real() != 0 && auto2[0].real() != 0) {
          double w = chanWidths[j] * itsTimeInterval;
          weight[0] *= w / (auto1[0].real() * auto2[0].real());  // XX
          if (npol == 4) {
            if (auto1[3].real() != 0 && auto2[3].real() != 0) {
              weight[1] *= w / (auto1[0].real() * auto2[3].real());  // XY
              weight[2] *= w / (auto1[3].real() * auto2[0].real());  // YX
              weight[3] *= w / (auto1[3].real() * auto2[3].real());  // YY
            }
          } else if (npol == 2) {
            if (auto1[1].real() != 0 && auto2[1].real() != 0) {
              weight[1] *= w / (auto1[1].real() * auto2[1].real());  // YY
            }
          }
        }
        // Set pointers to next channel.
        weight += npol;
        auto1 += npol;
        auto2 += npol;
      }
    } else {
      // No autocorrelations for this baseline, so skip it.
      weight += nchan * npol;
    }
  }
}

bool MSReader::getFullResFlags(const RefRows& rowNrs, DPBuffer& buf) {
  NSTimer::StartStop sstime(itsTimer);
  Cube<bool>& flags = buf.getFullResFlags();
  int norigchan = itsNrChan * itsFullResNChanAvg;
  // Resize if needed (probably when called for first time).
  if (flags.empty()) {
    flags.resize(norigchan, itsFullResNTimeAvg, itsNrBl);
  }
  // Return false if no fullRes flags available.
  if (!itsHasFullResFlags) {
    flags = false;
    return false;
  }
  // Flag everything if data rows are missing.
  if (rowNrs.rowVector().empty()) {
    flags = true;
    return true;
  }
  ArrayColumn<uChar> fullResFlagCol(itsMS, "LOFAR_FULL_RES_FLAG");
  int origstart = itsStartChan * itsFullResNChanAvg;
  Array<uChar> chars = fullResFlagCol.getColumnCells(rowNrs);
  // The original flags are kept per channel, not per corr.
  // Per row the flags are stored as uchar[nchar,navgtime].
  // Each char contains a bit per channel, thus nchan/8 chars are needed.
  // ntimeavg is the nr of times used when averaging.
  // Return it as Cube<bool>[norigchan,ntimeavg,nrbl].
  IPosition chShape = chars.shape();
  if (chShape[1] != itsFullResNTimeAvg || chShape[2] != itsNrBl)
    throw std::runtime_error("Incorrect shape of LOFAR_FULL_RES_FLAG column");
  // Now expand the bits to bools.
  // If all bits to convert are contiguous, do it all in one go.
  // Otherwise we have to iterate.
  if (norigchan == chShape[0] * 8) {
    Conversion::bitToBool(flags.data(), chars.data(), flags.size());
  } else {
    if (norigchan > chShape[0] * 8)
      throw std::runtime_error(
          "Shape of LOFAR_FULL_RES_FLAG column is inconsistent");
    const uChar* charsPtr = chars.data();
    bool* flagsPtr = flags.data();
    for (int i = 0; i < chShape[1] * chShape[2]; ++i) {
      Conversion::bitToBool(flagsPtr, charsPtr, origstart, norigchan);
      flagsPtr += norigchan;
      charsPtr += chShape[0];
    }
  }
  return true;
}

void MSReader::getModelData(const casacore::RefRows& rowNrs,
                            casacore::Cube<casacore::Complex>& arr) {
  NSTimer::StartStop sstime(itsTimer);
  if (rowNrs.rowVector().empty()) {
    arr.resize(itsNrCorr, itsNrChan, itsNrBl);
    arr = Complex();
  } else {
    ArrayColumn<Complex> modelCol(itsMS, itsModelColName);
    if (itsUseAllChan) {
      modelCol.getColumnCells(rowNrs, arr);
    } else {
      modelCol.getColumnCells(rowNrs, itsColSlicer, arr);
    }
  }
}

void MSReader::fillBeamInfo(
    vector<everybeam::Station::Ptr>& vec,
    const casacore::Vector<casacore::String>& antNames) {
  // Get the names of all stations in the MS.
  const Vector<String>& allNames = getInfo().antennaNames();
  // Create a vector holding the beam info of all stations.
  vector<everybeam::Station::Ptr> beams(allNames.size());
  everybeam::ReadStations(itsMS, beams.begin());
  // Copy only the ones for which the station name matches.
  // Note: the order of the station names in both vectors match.
  vec.resize(antNames.size());
  unsigned int ant = 0;
  for (unsigned int i = 0; i < allNames.size(); ++i) {
    if (ant < antNames.size() && allNames[i] == antNames[ant]) {
      vec[ant] = beams[i];
      ant++;
    }
  }
  if (ant != vec.size())
    throw Exception(
        "MSReader::fillBeamInfo -"
        " some stations miss the beam info");
}

}  // namespace DPPP
}  // namespace DP3
