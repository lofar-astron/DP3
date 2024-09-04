// MsReader.cc: DP3 step reading from an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MsReader.h"

#include <iostream>
#include <tuple>

#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>

#include <EveryBeam/load.h>
#include <EveryBeam/msreadutils.h>
#include <EveryBeam/telescope/phasedarray.h>

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

#include <aocommon/logger.h>

#include <dp3/common/Types.h>

#include "../base/BaselineSelection.h"
#include "../base/MS.h"
#include "../common/ParameterSet.h"

using casacore::ArrayColumn;
using casacore::ArrayMeasColumn;
using casacore::Block;
using casacore::Cube;
using casacore::IPosition;
using casacore::Matrix;
using casacore::MDirection;
using casacore::MeasTable;
using casacore::MeasureHolder;
using casacore::MeasurementSet;
using casacore::MPosition;
using casacore::MSAntennaParse;
using casacore::MSSelection;
using casacore::MSSelectionErrorHandler;
using casacore::MVTime;
using casacore::Quantity;
using casacore::Record;
using casacore::RecordGram;
using casacore::RefRows;
using casacore::ScalarColumn;
using casacore::ScalarMeasColumn;
using casacore::Slicer;
using casacore::Table;
using casacore::TableColumn;
using casacore::TableDesc;
using casacore::TableExprNode;
using casacore::TableIterator;
using casacore::TableLock;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

MsReader::MsReader(const casacore::MeasurementSet& ms,
                   const common::ParameterSet& parset,
                   const std::string& prefix, bool allow_missing_data)
    : ms_(ms),
      selection_ms_(ms_),
      data_column_name_(parset.getString(prefix + "datacolumn", "DATA")),
      extra_data_column_names_(parset.getStringVector(
          prefix + "extradatacolumns", std::vector<std::string>())),
      flag_column_name_(parset.getString(prefix + "flagcolumn", "FLAG")),
      weight_column_name_(
          parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM")),
      model_column_name_(
          parset.getString(prefix + "modelcolumn", "MODEL_DATA")),
      start_channel_expression_(parset.getString(prefix + "startchan", "0")),
      n_nchannels_expression_(parset.getString(prefix + "nchan", "0")),
      baseline_selection_(parset.getString(prefix + "baseline", std::string())),
      sort_(parset.getBool(prefix + "sort", false)),
      auto_weight_(parset.getBool(prefix + "autoweight", false)),
      force_auto_weight_(parset.getBool(prefix + "forceautoweight", false)),
      use_flags_(parset.getBool(prefix + "useflag", true)),
      time_tolerance_(parset.getDouble(prefix + "timetolerance", 1e-2)) {
  common::NSTimer::StartStop sstime(timer_);

  if (allow_missing_data && !extra_data_column_names_.empty()) {
    throw std::runtime_error(
        "Settings 'missingdata' and 'extradatacolumns' are mutually "
        "exclusive and cannot be provided both.");
  }

  if (allow_missing_data && ms.isNull()) {
    aocommon::Logger::Warn << "MeasurementSet is empty; dummy data used\n";
    return;
  }

  assert(!HasBda(ms));
  // See if a selection on band needs to be done.
  // We assume that DATA_DESC_ID and SPW_ID map 1-1.
  int spectralWindow = parset.getInt(prefix + "band", -1);
  if (spectralWindow >= 0) {
    aocommon::Logger::Info << " MsReader selecting spectral window "
                           << spectralWindow << " ...\n";
    Table subset =
        selection_ms_(selection_ms_.col("DATA_DESC_ID") == spectralWindow);
    // If not all is selected, use the selection.
    if (subset.nrow() < selection_ms_.nrow()) {
      if (subset.nrow() <= 0)
        throw std::runtime_error("Band " + std::to_string(spectralWindow) +
                                 " not found in " + msName());
      selection_ms_ = subset;
    }
  } else {
    spectralWindow = 0;
  }
  // See if a selection on baseline needs to be done.
  if (!baseline_selection_.empty()) {
    aocommon::Logger::Info << " MsReader selecting baselines ...\n";
    MSSelection select;

    // Overwrite the error handler to ignore errors for unknown antennas.
    // First construct MSSelection, because it resets the error handler.
    dp3::base::LogAntennaParseErrors ignore_unknown_antennas;

    // Set given selection strings.
    select.setAntennaExpr(baseline_selection_);
    // Create a table expression for an MS representing the selection.
    MeasurementSet ms(selection_ms_);
    TableExprNode node = select.toTableExprNode(&ms);
    Table subset = selection_ms_(node);
    // If not all is selected, use the selection.
    if (subset.nrow() < selection_ms_.nrow()) {
      if (subset.nrow() <= 0)
        throw std::runtime_error("Baselines " + baseline_selection_ +
                                 "not found in " + msName());
      selection_ms_ = subset;
    }
  }
  // Prepare the MS access and store MS metadata into infoOut().
  prepare(allow_missing_data);

  ParseTimeSelection(parset, prefix);

  unsigned int start_channel, n_channels;
  std::tie(start_channel, n_channels) = ParseChannelSelection(
      start_channel_expression_, n_nchannels_expression_, getInfoOut().nchan());
  // Are all channels used?
  use_all_channels_ = start_channel == 0 && n_channels == getInfoOut().nchan();

  // Do the rest of the preparation.
  prepare2(spectralWindow, start_channel, n_channels);
  // Take subset of channel frequencies if needed.
  // Make sure to copy the subset to get a proper Vector.
  // Form the slicer to get channels and correlations from column.
  column_slicer_ = Slicer(IPosition(2, 0, start_channel),
                          IPosition(2, getInfoOut().ncorr(), n_channels));
  // Form the slicer to get channels, corrs, and baselines from array.
  array_slicer_ = Slicer(IPosition(3, 0, start_channel, 0),
                         IPosition(3, getInfoOut().ncorr(), n_channels,
                                   getInfoOut().nbaselines()));
  // Initialize the flag counters.
  flag_counter_.init(getInfoOut());
}

std::string MsReader::msName() const { return ms_.tableName(); }

bool MsReader::process(std::unique_ptr<DPBuffer> buffer) {
  const std::array<std::size_t, 3> shape{
      getInfoOut().nbaselines(), getInfoOut().nchan(), getInfoOut().ncorr()};

  // Determine what to read. Depends on what later steps in the chain need.
  if (getFieldsToRead().Data()) {
    buffer->GetData().resize(shape);
    for (std::string columnName : extra_data_column_names_) {
      // The extra data buffer could already exist, for instance because
      // MultiMsReader reuses existing buffers.
      if (!buffer->HasData(columnName)) {
        buffer->AddData(columnName);
      } else {
        buffer->GetData(columnName).resize(shape);
      }
    }
  }
  if (getFieldsToRead().Flags()) {
    buffer->GetFlags().resize(shape);
  }

  double corrected_first_time = getInfoOut().firstTime() + time_correction_;

  {
    common::NSTimer::StartStop sstime(timer_);
    // Use time from the current time slot in the MS.
    bool useIter = false;
    while (!ms_iterator_.pastEnd()) {
      // Take time from row 0 in subset.
      double mstime = ScalarColumn<double>(ms_iterator_.table(), "TIME")(0);
      // Skip time slot and give warning if MS data is not in time order.
      if (mstime < previous_time_) {
        aocommon::Logger::Warn
            << "Time at rownr " << ms_iterator_.table().rowNumbers(ms_)[0]
            << " of MS " << msName() << " is less than previous time slot\n";
      } else {
        // Use the time slot if near nexttime or between [starttime, nexttime].
        // In this way we cater for irregular times in some WSRT MSs.
        if (casacore::nearAbs(mstime, next_time_, time_tolerance_)) {
          useIter = true;
          break;
        } else if (mstime > corrected_first_time && mstime < next_time_) {
          time_correction_ -= next_time_ - mstime;
          corrected_first_time = getInfoOut().firstTime() + time_correction_;

          next_time_ = mstime;
          useIter = true;
          break;
        }
        if (mstime > next_time_) {
          // A time slot seems to be missing; insert one.
          break;
        }
      }
      // Skip this time slot.
      previous_time_ = mstime;
      ms_iterator_.next();
    }

    // Stop if at the end, i.e. the above loop completed without hitting an end
    // condition at all, or if there is no data at all.
    if ((next_time_ > getInfoOut().lastTime() &&
         !casacore::near(next_time_, getInfoOut().lastTime())) ||
        next_time_ == 0.0) {
      return false;
    }

    // Fill the buffer.
    buffer->SetTime(next_time_);
    if (!useIter) {
      // Time slot is missing altogether: insert a fully flagged time slot.
      buffer->SetRowNumbers(casacore::Vector<common::rownr_t>());
      buffer->SetExposure(getInfoOut().timeInterval());
      buffer->GetFlags().fill(true);
      if (getFieldsToRead().Data()) {
        buffer->GetData().fill(std::complex<float>());
        for (std::string columnName : extra_data_column_names_) {
          buffer->GetData(columnName).fill(std::complex<float>());
        }
      }
      n_inserted_++;
    } else {
      // Timeslot exists, but it could still be that the MS does not contain a
      // data column (visibilities) for this timeslot.
      buffer->SetRowNumbers(ms_iterator_.table().rowNumbers(ms_, true));
      if (missing_data_) {
        // Data column not present, so fill a fully flagged time slot.
        buffer->SetExposure(getInfoOut().timeInterval());
        buffer->GetFlags().fill(true);
        if (getFieldsToRead().Data()) {
          buffer->GetData().fill(std::complex<float>());
        }
      } else {
        // Set exposure.
        buffer->SetExposure(
            ScalarColumn<double>(ms_iterator_.table(), "EXPOSURE")(0));
        // Get data (visibilities) and flags from the MS.
        const casacore::IPosition casa_shape(3, shape[2], shape[1], shape[0]);
        if (getFieldsToRead().Data()) {
          ArrayColumn<casacore::Complex> dataCol(ms_iterator_.table(),
                                                 data_column_name_);
          casacore::Cube<casacore::Complex> casa_data(
              casa_shape, buffer->GetData().data(), casacore::SHARE);
          if (use_all_channels_) {
            dataCol.getColumn(casa_data);
          } else {
            dataCol.getColumn(column_slicer_, casa_data);
          }
        }
        if (getFieldsToRead().Flags()) {
          if (use_flags_) {
            ArrayColumn<bool> flagCol(ms_iterator_.table(), flag_column_name_);
            casacore::Cube<bool> casa_flags(
                casa_shape, buffer->GetFlags().data(), casacore::SHARE);

            if (use_all_channels_) {
              flagCol.getColumn(casa_flags);
            } else {
              flagCol.getColumn(column_slicer_, casa_flags);
            }
            // Set flags if FLAG_ROW is set.
            ScalarColumn<bool> flagrowCol(ms_iterator_.table(), "FLAG_ROW");
            for (unsigned int i = 0; i < ms_iterator_.table().nrow(); ++i) {
              if (flagrowCol(i)) {
                casa_flags(IPosition(3, 0, 0, i),
                           IPosition(3, getInfoOut().ncorr() - 1,
                                     getInfoOut().nchan() - 1, i)) = true;
              }
            }

          } else {
            // Do not use FLAG from the MS.
            buffer->GetFlags().fill(false);
          }
          // Flag invalid data (NaN, infinite).
          FlagInfinityNan(*buffer, flag_counter_);
        }
      }

      // Get extra data (e.g. model data) from the MS.
      if (getFieldsToRead().Data()) {
        const casacore::IPosition casa_shape(3, shape[2], shape[1], shape[0]);
        for (std::string columnName : extra_data_column_names_) {
          ArrayColumn<casacore::Complex> extraDataCol(ms_iterator_.table(),
                                                      columnName);
          casacore::Cube<casacore::Complex> casa_data(
              casa_shape, buffer->GetData(columnName).data(), casacore::SHARE);
          if (use_all_channels_) {
            extraDataCol.getColumn(casa_data);
          } else {
            extraDataCol.getColumn(column_slicer_, casa_data);
          }
        }
      }

      previous_time_ = next_time_;
      n_read_++;
      ms_iterator_.next();
    }
    if ((getFieldsToRead().Flags() &&
         (buffer->GetFlags().shape(0) != shape[0])) ||
        (useIter && (buffer->GetRowNumbers().shape() != shape[0]))) {
      throw std::runtime_error(
          "#baselines is not the same for all time slots in the MS");
    }
  }  // end of scope stops the timer.

  if (getFieldsToRead().Uvw())
    getUVW(buffer->GetRowNumbers(), buffer->GetTime(), *buffer);
  if (getFieldsToRead().Weights()) GetWeights(buffer->GetRowNumbers(), *buffer);

  getNextStep()->process(std::move(buffer));
  // Do not add to previous time, because it introduces round-off errors.
  next_time_ = corrected_first_time +
               (n_read_ + n_inserted_) * getInfoOut().timeInterval();
  return true;
}

void MsReader::FlagInfinityNan(DPBuffer& buffer, FlagCounter& flagCounter) {
  const int ncorr = buffer.GetData().shape(2);
  const std::complex<float>* dataPtr = buffer.GetData().data();
  bool* flagPtr = buffer.GetFlags().data();
  for (unsigned int i = 0; i < buffer.GetData().size(); i += ncorr) {
    for (unsigned int j = i; j < i + ncorr; ++j) {
      bool flag = (!std::isfinite(dataPtr[j].real()) ||
                   !std::isfinite(dataPtr[j].imag()));
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
  }
}

void MsReader::finish() { getNextStep()->finish(); }

void MsReader::show(std::ostream& os) const {
  os << "MsReader\n";
  os << "  input MS:           " << msName() << '\n';
  if (ms_.isNull()) {
    os << "    *** MS does not exist ***\n";
  } else {
    if (!baseline_selection_.empty()) {
      os << "  baseline:           " << baseline_selection_ << '\n';
    }
    os << "  band                " << getInfoOut().spectralWindow() << '\n';
    os << "  startchan:          " << getInfoOut().startchan() << "  ("
       << start_channel_expression_ << ")\n";
    os << "  nchan:              " << getInfoOut().nchan() << "  ("
       << n_nchannels_expression_ << ")\n";
    os << "  ncorrelations:      " << getInfoOut().ncorr() << '\n';
    unsigned int nrbl = getInfoOut().nbaselines();
    os << "  nbaselines:         " << nrbl << '\n';
    os << "  first time:         " << MVTime::Format(MVTime::YMD)
       << MVTime(getInfoOut().firstTime() / (24 * 3600.)) << '\n';
    os << "  maximum time:       " << MVTime::Format(MVTime::YMD)
       << MVTime(getInfoOut().lastTime() / (24 * 3600.)) << '\n';
    os << "  ntimes:             " << getInfoOut().ntime()
       << '\n';  // selection_ms_ can contain timeslots that are ignored in
                 // process
    os << "  time interval:      " << getInfoOut().timeInterval() << '\n';
    os << "  DATA column:        " << data_column_name_;
    if (missing_data_) {
      os << "  (not present)";
    }
    os << '\n';
    os << "  extra data columns: ";
    for (std::string columnName : extra_data_column_names_) {
      os << columnName << ", ";
    }
    os << '\n';
    os << "  WEIGHT column:      " << weight_column_name_ << '\n';
    os << "  FLAG column:        " << flag_column_name_ << '\n';
    os << "  autoweight:         " << std::boolalpha << auto_weight_ << '\n';
  }
}

void MsReader::showCounts(std::ostream& os) const {
  os << '\n' << "NaN/infinite data flagged in reader";
  os << '\n' << "===================================" << '\n';
  int64_t nrtim = n_read_;
  flag_counter_.showCorrelation(os, nrtim);
  os << n_inserted_ << " missing time slots were inserted" << '\n';
}

void MsReader::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MsReader" << '\n';
}

void MsReader::ParseTimeSelection(const common::ParameterSet& parset,
                                  const std::string& prefix) {
  // Start and end time can be given in the parset in case leading
  // or trailing time slots are missing.
  // They can also be used to select part of the MS.

  // Get info from parset.
  const std::string startTimeStr = parset.getString(prefix + "starttime", "");
  const std::string endTimeStr = parset.getString(prefix + "endtime", "");
  const unsigned int nTimes = parset.getInt(prefix + "ntimes", 0);
  const int startTimeSlot = parset.getInt(prefix + "starttimeslot", 0);

  // Get time properties from the MS, which prepare() already set.
  double first_time = getInfoOut().firstTime();
  double last_time = getInfoOut().lastTime();
  const double interval = getInfo().timeInterval();

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
    if (startTimeParset > last_time)
      throw std::runtime_error("Specified starttime is past end of time axis");

    // Round specified first time to a multiple of 'interval'.
    first_time +=
        std::ceil((startTimeParset - first_time) / interval) * interval;
  } else {
    first_time += startTimeSlot * interval;
  }

  if (!endTimeStr.empty()) {
    Quantity qtime;
    if (!MVTime::read(qtime, endTimeStr)) {
      throw std::runtime_error(endTimeStr + " is an invalid date/time");
    }
    const double endTimeParset = qtime.getValue("s");

    // Some overlap between the measurement set timerange and the parset range
    // is required :
    if (endTimeParset < getInfoOut().startTime()) {
      throw std::runtime_error(
          "Specified end time " + endTimeStr +
          " is before the first timestep in the measurement set");
    }
    // Round specified last time to a multiple of 'interval'.
    const double ms_first_time = getInfoOut().firstTime();
    last_time =
        ms_first_time +
        std::floor((endTimeParset - ms_first_time) / interval) * interval;
  }

  if (last_time < first_time)
    throw std::runtime_error("Specified endtime is before specified starttime");
  // If needed, skip the first times in the MS.
  // It also sets first_time properly (round to time/interval in MS).
  SkipFirstTimes(first_time, interval);
  if (nTimes > 0) {
    if (!endTimeStr.empty()) {
      throw std::runtime_error("Only one of " + prefix + "ntimes and " +
                               prefix + "endtime can be specified");
    }
    last_time = first_time + (nTimes - 1) * interval;
  }

  infoOut().setTimes(first_time, last_time, interval);
  next_time_ = first_time;
}

std::pair<unsigned int, unsigned int> MsReader::ParseChannelSelection(
    const std::string& start_channel_string,
    const std::string& n_channels_string, unsigned int n_all_channels) {
  Record rec;
  rec.define("nchan", n_all_channels);
  TableExprNode node1(RecordGram::parse(rec, start_channel_string));
  TableExprNode node2(RecordGram::parse(rec, n_channels_string));
  double result;
  node1.get(rec, result);
  const unsigned int start_channel = result + 0.001;
  node2.get(rec, result);
  unsigned int n_channels = result + 0.0001;
  if (start_channel >= n_all_channels)
    throw std::runtime_error("startchan " + std::to_string(start_channel) +
                             " exceeds nr of channels in MS (" +
                             std::to_string(n_all_channels) + ')');
  const unsigned int max_n_channels = n_all_channels - start_channel;
  if (n_channels == 0) {  // nchan=0 means until the last channel.
    n_channels = max_n_channels;
  } else {
    n_channels = std::min(n_channels, max_n_channels);
  }

  return std::make_pair(start_channel, n_channels);
}

void MsReader::prepare(const bool allow_missing_data) {
  // Find the number of correlations and channels.
  IPosition shape(
      ArrayColumn<casacore::Complex>(selection_ms_, "DATA").shape(0));
  const unsigned int n_correlations = shape[0];
  const unsigned int n_channels = shape[1];
  info() = DPInfo(n_correlations, n_channels, base::ReadAntennaSet(ms_));

  if (selection_ms_.nrow() == 0) {
    aocommon::Logger::Warn << "The selected input does not contain any data.\n";
  }
  TableDesc tdesc = ms_.tableDesc();

  has_weight_spectrum_ = false;
  // if weightcolname is specified to "WEIGHT" then this is used, even
  // if a weight_spectrum is present.
  if (weight_column_name_ != "WEIGHT") {
    // Test if specified weight column or WEIGHT_SPECTRUM is present.
    if (tdesc.isColumn(weight_column_name_)) {
      // The column is there, but it might not contain values. Test row 0.
      has_weight_spectrum_ =
          ArrayColumn<float>(selection_ms_, weight_column_name_).isDefined(0);
      if (!has_weight_spectrum_ && weight_column_name_ != "WEIGHT_SPECTRUM") {
        aocommon::Logger::Warn
            << "Specified weight column " << weight_column_name_
            << "is not a valid column, using WEIGHT instead\n";
      }
    }
  }

  // Test if the data column is present.
  missing_data_ = !tdesc.isColumn(data_column_name_);
  if (!missing_data_) {
    // Read beam keywords of input datacolumn
    ArrayColumn<casacore::Complex> dataCol(ms_, data_column_name_);
    if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
      const everybeam::CorrectionMode mode = everybeam::ParseCorrectionMode(
          dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE"));
      info().setBeamCorrectionMode(static_cast<int>(mode));
      if (mode != everybeam::CorrectionMode::kNone) {
        casacore::String error;
        MeasureHolder mHolder;
        if (!mHolder.fromRecord(
                error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
          throw std::runtime_error(error);
        info().setBeamCorrectionDir(mHolder.asMDirection());
      }
    }
  } else if (allow_missing_data) {
    // Only give warning if a missing data column is allowed.
    aocommon::Logger::Warn << "Data column " << data_column_name_
                           << " is missing in " << msName() << '\n';
  } else {
    throw std::runtime_error("Data column " + data_column_name_ +
                             " is missing in " + msName());
  }

  // Test if the extra data columns are present (e.g. model visibilities).
  std::string missing_columns;
  for (std::string columnName : extra_data_column_names_) {
    if (!tdesc.isColumn(columnName)) {
      missing_columns += columnName + ", ";
    }
  }
  if (!missing_columns.empty()) {
    missing_columns.erase(missing_columns.size() - 2);
    throw std::runtime_error("Extra data columns [" + missing_columns +
                             "] are missing in " + msName());
  }

  // Get the main table in the correct order.
  // Determine if the data are stored using LofarStMan.
  // If so, we know it is in time order.
  // (sorting on TIME with LofarStMan can be expensive).
  bool needSort = sort_;
  bool useRaw = false;
  Record dminfo = ms_.dataManagerInfo();
  for (unsigned i = 0; i < dminfo.nfields(); ++i) {
    Record subrec = dminfo.subRecord(i);
    if (subrec.asString("TYPE") == "LofarStMan") {
      needSort = false;
      useRaw = true;
      break;
    }
  }
  // Give an error if autoweight is used for a non-raw MS.
  if (force_auto_weight_) {
    auto_weight_ = true;
  } else if (!useRaw && auto_weight_) {
    throw std::runtime_error(
        "Using autoweight=true cannot be done on DP3-ed MS");
  }
  // If not in order, sort the table selection (also on baseline).
  Table sortms(selection_ms_);
  Block<casacore::String> sortCols(3);
  sortCols[0] = "TIME";
  sortCols[1] = "ANTENNA1";
  sortCols[2] = "ANTENNA2";
  if (needSort) {
    sortms = selection_ms_.sort(sortCols);
  }
  // Get first and last time and interval from MS.
  if (selection_ms_.nrow() > 0) {
    ScalarColumn<double> time_column(sortms, "TIME");
    info().setTimes(time_column(0), time_column(sortms.nrow() - 1),
                    ScalarColumn<double>(sortms, "INTERVAL")(0));
  }
  // Create iterator over time. Do not sort again.
  ms_iterator_ = TableIterator(sortms, Block<casacore::String>(1, "TIME"),
                               TableIterator::Ascending, TableIterator::NoSort);
  const common::rownr_t n_baselines = ms_iterator_.table().nrow();

  // Ensure we have only one band by checking the nr of unique baselines.
  Table sortab = ms_iterator_.table().sort(
      sortCols, casacore::Sort::Ascending,
      casacore::Sort::QuickSort + casacore::Sort::NoDuplicates);
  if (sortab.nrow() != n_baselines)
    throw std::runtime_error("The MS appears to have multiple subbands");
  // Get the baseline columns.
  ScalarColumn<int> ant1col(ms_iterator_.table(), "ANTENNA1");
  ScalarColumn<int> ant2col(ms_iterator_.table(), "ANTENNA2");
  if (ant1col.nrow() != n_baselines || ant2col.nrow() != n_baselines) {
    throw std::runtime_error("Antenna column(s) do not match baseline count");
  }
  // Get the antenna names and positions.
  Table anttab(ms_.keywordSet().asTable("ANTENNA"));
  ScalarColumn<casacore::String> nameCol(anttab, "NAME");
  ScalarColumn<double> diamCol(anttab, "DISH_DIAMETER");
  unsigned int nant = anttab.nrow();
  ScalarMeasColumn<MPosition> antcol(anttab, "POSITION");
  std::vector<MPosition> antPos;
  antPos.reserve(nant);
  for (unsigned int i = 0; i < nant; ++i) {
    antPos.push_back(antcol(i));
  }
  // Set antenna/baseline info.
  casacore::Vector<casacore::String> names = nameCol.getColumn();
  infoOut().setAntennas(std::vector<std::string>(names.begin(), names.end()),
                        diamCol.getColumn().tovector(), antPos,
                        ant1col.getColumn().tovector(),
                        ant2col.getColumn().tovector());

  // Read the phase reference position from the FIELD subtable.
  // Only use the main value from the PHASE_DIR array.
  // The same for DELAY_DIR and LOFAR_TILE_BEAM_DIR.
  // If LOFAR_TILE_BEAM_DIR does not exist, use DELAY_DIR.
  Table fldtab(ms_.keywordSet().asTable("FIELD"));
  if (fldtab.nrow() != 1)
    throw std::runtime_error("Multiple entries in FIELD table");
  ArrayMeasColumn<MDirection> fldcol1(fldtab, "PHASE_DIR");
  ArrayMeasColumn<MDirection> fldcol2(fldtab, "DELAY_DIR");
  const MDirection phaseCenter = *(fldcol1(0).data());
  const MDirection delayCenter = *(fldcol2(0).data());

  MDirection tileBeamDir;
  try {
    tileBeamDir = everybeam::ReadTileBeamDirection(ms_);
  } catch (const std::runtime_error& error) {
    // everybeam throws an exception error if telescope != [LOFAR, AARTFAAC]
    // in that case, default back to "DELAY_DIR"
    tileBeamDir = *(fldcol2(0).data());
  }

  // Get the array position using the telescope name from the OBSERVATION
  // subtable.
  const casacore::Table observation_table(
      ms_.keywordSet().asTable(base::DP3MS::kObservationTable));
  ScalarColumn<casacore::String> telCol(observation_table, "TELESCOPE_NAME");
  MPosition arrayPos;
  if (observation_table.nrow() == 0 ||
      !MeasTable::Observatory(arrayPos, telCol(0))) {
    // If not found, use the position of the middle antenna.
    arrayPos = antPos[antPos.size() / 2];
  }
  info().setArrayInformation(arrayPos, phaseCenter, delayCenter, tileBeamDir);
  // Create the UVW calculator.
  uvw_calculator_ =
      std::make_unique<base::UVWCalculator>(phaseCenter, arrayPos, antPos);
}

void MsReader::prepare2(int spectralWindow, unsigned int start_channel,
                        unsigned int n_channels) {
  infoOut().setMsNames(msName(), data_column_name_, flag_column_name_,
                       weight_column_name_);
  // Read the center frequencies of all channels.
  Table spwtab(ms_.keywordSet().asTable("SPECTRAL_WINDOW"));
  ArrayColumn<double> freqCol(spwtab, "CHAN_FREQ");
  ArrayColumn<double> widthCol(spwtab, "CHAN_WIDTH");
  ArrayColumn<double> resolCol(spwtab, "RESOLUTION");
  ArrayColumn<double> effBWCol(spwtab, "EFFECTIVE_BW");
  ScalarColumn<double> refCol(spwtab, "REF_FREQUENCY");
  std::vector<double> chanFreqs = freqCol(spectralWindow).tovector();
  std::vector<double> chanWidths = widthCol(spectralWindow).tovector();
  std::vector<double> resolutions = resolCol(spectralWindow).tovector();
  std::vector<double> effectiveBW = effBWCol(spectralWindow).tovector();
  const double refFreq = refCol(spectralWindow);
  infoOut().setChannels(std::move(chanFreqs), std::move(chanWidths),
                        std::move(resolutions), std::move(effectiveBW), refFreq,
                        spectralWindow);
  infoOut().SelectChannels(start_channel, n_channels);

  std::set<aocommon::PolarizationEnum> polarizations;
  casacore::MSDataDescription data_description_table = ms_.dataDescription();
  casacore::ScalarColumn<int> polarization_index_column(
      data_description_table,
      casacore::MSDataDescription::columnName(
          casacore::MSDataDescription::POLARIZATION_ID));
  const size_t polarization_index = polarization_index_column(spectralWindow);
  // Populate the polarization
  casacore::MSPolarization pol_table = ms_.polarization();
  casacore::ArrayColumn<int> corr_type_column(
      pol_table, casacore::MSPolarization::columnName(
                     casacore::MSPolarizationEnums::CORR_TYPE));
  casacore::Array<int> corr_type_vec(corr_type_column(polarization_index));
  for (casacore::Array<int>::const_contiter p = corr_type_vec.cbegin();
       p != corr_type_vec.cend(); ++p) {
    polarizations.emplace(aocommon::Polarization::AipsIndexToEnum(*p));
  }
  if (polarizations.size() != 4) {
    throw std::runtime_error(
        "DP3 expects a measurement set with 4 polarizations");
  }
  info().setPolarizations(polarizations);
}

void MsReader::SkipFirstTimes(double& first_time, const double interval) {
  while (!ms_iterator_.pastEnd()) {
    // Take time from row 0 in subset.
    double mstime = ScalarColumn<double>(ms_iterator_.table(), "TIME")(0);
    // Skip time slot and give warning if MS data is not in time order.
    if (mstime < previous_time_) {
      aocommon::Logger::Warn
          << "Time at rownr " << ms_iterator_.table().rowNumbers(ms_)[0]
          << " of MS " << msName() << " is less than previous time slot\n";
    } else {
      // Stop skipping if time equal to itsFirstTime.
      if (casacore::near(mstime, first_time)) {
        first_time = mstime;
        break;
      }
      // Also stop if time > first_time.
      // In that case determine the true first time, because first_time
      // can be a time value that does not coincide with a true time.
      // Note that a time stamp might be missing at this point,
      // so do not simply assume that mstime can be used.
      if (mstime > first_time) {
        int nrt = int((mstime - first_time) / interval);
        mstime -= (nrt + 1) * interval;  // Add 1 for rounding errors
        if (casacore::near(mstime, first_time)) {
          first_time = mstime;
        } else {
          first_time = mstime + interval;
        }
        break;
      }
    }
    // Skip this time slot.
    previous_time_ = mstime;
    ms_iterator_.next();
  }
}

void MsReader::getUVW(const RefRows& rowNrs, double time, DPBuffer& buf) {
  common::NSTimer::StartStop sstime(timer_);
  const unsigned int n_baselines = getInfoOut().nbaselines();
  buf.GetUvw().resize({n_baselines, 3});
  if (rowNrs.rowVector().empty()) {
    // Calculate UVWs if empty rownrs (i.e., missing data).
    const std::vector<int>& ant1 = getInfo().getAnt1();
    const std::vector<int>& ant2 = getInfo().getAnt2();
    for (unsigned int i = 0; i < n_baselines; ++i) {
      xt::view(buf.GetUvw(), i, xt::all()) =
          xt::adapt(uvw_calculator_->getUVW(ant1[i], ant2[i], time));
    }
  } else {  // Load UVW from MS
    ArrayColumn<double> dataCol(ms_, "UVW");
    const casacore::IPosition shape(2, 3, n_baselines);
    casacore::Matrix<double> casa_uvw(shape, buf.GetUvw().data(),
                                      casacore::SHARE);
    dataCol.getColumnCells(rowNrs, casa_uvw);
  }
}

void MsReader::GetWeights(const RefRows& rowNrs, DPBuffer& buf) {
  common::NSTimer::StartStop sstime(timer_);
  const unsigned int n_baselines = getInfoOut().nbaselines();
  const unsigned int n_channels = getInfoOut().nchan();
  const unsigned int n_correlations = getInfoOut().ncorr();
  // Resize if needed (probably when called for first time).
  buf.GetWeights().resize({n_baselines, n_channels, n_correlations});
  DPBuffer::WeightsType& weights = buf.GetWeights();
  const casacore::IPosition shape(3, n_correlations, n_channels, n_baselines);
  casacore::Cube<float> casa_weights(shape, weights.data(), casacore::SHARE);
  if (rowNrs.rowVector().empty()) {
    // rowNrs can be empty if a time slot was inserted (i.e., missing data).
    weights.fill(0.0f);
  } else {
    // Get weights for entire spectrum if present in MS.
    if (has_weight_spectrum_) {
      ArrayColumn<float> wsCol(ms_, weight_column_name_);
      // Using getColumnCells(rowNrs,column_slicer_) fails for LofarStMan.
      // Hence work around it.
      if (use_all_channels_) {
        wsCol.getColumnCells(rowNrs, casa_weights);
      } else {
        Cube<float> w = wsCol.getColumnCells(rowNrs);
        casa_weights = w(array_slicer_);
      }
    } else {
      // No spectrum present; get global weights and assign to each channel.
      ArrayColumn<float> wCol(ms_, "WEIGHT");
      Matrix<float> inArr = wCol.getColumnCells(rowNrs);
      float* inPtr = inArr.data();
      float* outPtr = weights.data();
      for (unsigned int i = 0; i < n_baselines; ++i) {
        // Set global weights to 1 if zero. Some old MSs need that.
        for (unsigned int k = 0; k < n_correlations; ++k) {
          if (inPtr[k] == 0.) {
            inPtr[k] = 1.;
          }
        }
        for (unsigned int j = 0; j < n_channels; ++j) {
          for (unsigned int k = 0; k < n_correlations; ++k) {
            *outPtr++ = inPtr[k];
          }
        }
        inPtr += n_correlations;
      }
    }
    if (auto_weight_) {
      // Adapt weights using autocorrelations.
      AutoWeight(buf);
    }
  }
}

void MsReader::AutoWeight(DPBuffer& buf) {
  const double* chanWidths = getInfo().chanWidths().data();
  DPBuffer::WeightsType& weights = buf.GetWeights();
  const unsigned int nbl = weights.shape(0);
  const unsigned int nchan = weights.shape(1);
  const unsigned int npol = weights.shape(2);
  // Get the autocorrelations indices.
  const std::vector<int>& autoInx = getInfo().getAutoCorrIndex();
  // Calculate the weight for each cross-correlation data point.
  const std::vector<int>& ant1 = getInfo().getAnt1();
  const std::vector<int>& ant2 = getInfo().getAnt2();
  const DPBuffer::DataType& data = buf.GetData();
  for (unsigned int bl = 0; bl < nbl; ++bl) {
    // Can only be done if both autocorrelations are present.
    if (autoInx[ant1[bl]] >= 0 && autoInx[ant2[bl]] >= 0) {
      for (unsigned int chan = 0; chan < nchan; ++chan) {
        // Get offset of both autocorrelations in data array.
        const std::complex<float>* auto1 = &data(autoInx[ant1[bl]], chan, 0);
        const std::complex<float>* auto2 = &data(autoInx[ant2[bl]], chan, 0);
        float* weight = &weights(bl, chan, 0);
        if (auto1[0].real() != 0 && auto2[0].real() != 0) {
          double w = chanWidths[chan] * getInfoOut().timeInterval();
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
      }
    }
  }
}

}  // namespace steps
}  // namespace dp3
