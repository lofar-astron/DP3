// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MSBDAReader.h"

#include <algorithm>
#include <iostream>

#include <dp3/base/BDABuffer.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../base/MS.h"

#include "../common/ParameterSet.h"

#include <casacore/casa/OS/HostInfo.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>
#include <casacore/casa/version.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/OS/Conversion.h>

#include <EveryBeam/load.h>
#include <EveryBeam/msreadutils.h>
#include <EveryBeam/telescope/phasedarray.h>

#include <aocommon/logger.h>

using casacore::ArrayColumn;
using casacore::ArrayMeasColumn;
using casacore::Complex;
using casacore::Cube;
using casacore::IPosition;
using casacore::MDirection;
using casacore::MeasTable;
using casacore::MeasurementSet;
using casacore::MPosition;
using casacore::MS;
using casacore::MVTime;
using casacore::RefRows;
using casacore::ScalarColumn;
using casacore::ScalarMeasColumn;
using casacore::Table;
using casacore::TableColumn;
using casacore::TableDesc;
using casacore::TableLock;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using std::vector;

namespace dp3 {
namespace steps {

MSBDAReader::MSBDAReader(const casacore::MeasurementSet& ms,
                         const common::ParameterSet& parset,
                         const std::string& prefix)
    : ms_(ms),
      data_column_name_(
          parset.getString(prefix + "data_column", MS::columnName(MS::DATA))),
      weight_column_name_(parset.getString(
          prefix + "weightcolumn", MS::columnName(MS::WEIGHT_SPECTRUM))),
      last_ms_time_(0),
      last_ms_interval_(0),
      nread_(0),
      timer_(),
      pool_size_(0),
      desc_id_to_nchan_(),
      bl_to_id_() {
  // InputStep::CreateReader only creates MSBDAReader's for BDA MS's, so it's
  // a programming error, and not a user error, if the MS doesn't have BDA data.
  assert(HasBda(ms));

  const int spectral_window = parset.getInt(prefix + "band", -1);
  if (spectral_window > 0) {
    throw std::invalid_argument(
        "BDA in combination with multiple spectral windows is not implemented");
  }

  const unsigned int nchan = parset.getInt(prefix + "nchan", 0);
  if (nchan > 0) {
    throw std::invalid_argument(
        "BDA in combination with channel filtering is not implemented. Remove "
        "'nchan' from input parset");
  }
  const unsigned int start_chan = parset.getInt(prefix + "startchan", 0);
  if (start_chan > 0) {
    throw std::invalid_argument(
        "BDA in combination with channel filtering is not implemented. Remove "
        "'startchan' from input parset");
  }
  const unsigned int ntimes = parset.getInt(prefix + "ntimes", 0);
  if (ntimes > 0) {
    throw std::invalid_argument(
        "BDA in combination with time filtering is not implemented. Remove "
        "'ntimes' from input parset");
  }
}

MSBDAReader::~MSBDAReader() {}

std::string MSBDAReader::msName() const { return ms_.tableName(); }

bool MSBDAReader::process(std::unique_ptr<base::DPBuffer>) {
  return process(std::unique_ptr<base::BDABuffer>());
}

bool MSBDAReader::process(std::unique_ptr<base::BDABuffer>) {
  common::NSTimer::StartStop sstime(timer_);

  // TODO: Pre-calculate actual required pool size beforehand.
  auto buffer = std::make_unique<base::BDABuffer>(
      info().nbaselines() * info().nchan() * info().ncorr());

  ScalarColumn<int> ant1_col(ms_, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(ms_, MS::columnName(MS::ANTENNA2));
  ArrayColumn<float> weights_col(ms_, weight_column_name_);
  ArrayColumn<casacore::Complex> data_col(ms_, data_column_name_);
  ArrayColumn<double> uvw_col(ms_, MS::columnName(MS::UVW));
  ScalarColumn<double> time_col(ms_, MS::columnName(MS::TIME));
  ScalarColumn<double> interval_col(ms_, MS::columnName(MS::INTERVAL));
  ScalarColumn<double> exposure_col(ms_, MS::columnName(MS::EXPOSURE));
  ScalarColumn<int> data_desc_id_col(ms_, MS::columnName(MS::DATA_DESC_ID));

  // Cache the data that will be add to the buffer
  RefRows cell_range{
      nread_,
      std::min(ms_.nrow() - 1, common::rownr_t(nread_ + info().nbaselines()))};

  auto time = time_col.getColumnCells(cell_range);
  auto interval = interval_col.getColumnCells(cell_range);
  auto exposure = exposure_col.getColumnCells(cell_range);
  auto data_desc_id = data_desc_id_col.getColumnCells(cell_range);

  unsigned i = 0;
  while (nread_ < ms_.nrow() && i < info().nbaselines()) {
    const double ms_time = time[i];
    const double ms_interval = interval[i];

    if (ms_time + ms_interval / 2 <
        last_ms_time_ + last_ms_interval_ / 2 - 0.001) {
      aocommon::Logger::Warn << "Time at rownr " << nread_ << " of MS "
                             << msName()
                             << " is less than previous time slot\n";
      ++i;
      ++nread_;
      continue;
    }

    const auto ant12 = std::make_pair(ant1_col(nread_), ant2_col(nread_));
    casacore::Complex* data_ptr = nullptr;
    std::vector<casacore::Complex> data;
    if (getFieldsToRead().Data()) {
      data = data_col.get(nread_).tovector();
      data_ptr = data.data();
    }
    Cube<float> weights = weights_col.get(nread_);
    Cube<double> uvw = uvw_col.get(nread_);

    const bool success = buffer->AddRow(
        ms_time, interval[i], exposure[i], bl_to_id_[ant12],
        desc_id_to_nchan_[data_desc_id[i]], info().ncorr(), data_ptr, nullptr,
        weights.tovector().data(), nullptr, uvw.tovector().data());
    (void)success;
    assert(success);  // The buffer should always be large enough.

    last_ms_time_ = ms_time;
    last_ms_interval_ = ms_interval;
    ++i;
    ++nread_;
  }

  getNextStep()->process(std::move(buffer));

  // Return true while there are still items remaining
  return nread_ < ms_.nrow();
}

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::updateInfo(const DPInfo& info_in) {
  using MS_Ant = casacore::MSAntenna;
  using MS_SPW = casacore::MSSpectralWindow;
  using MS_Field = casacore::MSField;
  using MS_Obs = casacore::MSObservation;

  Table factors = ms_.keywordSet().asTable(base::DP3MS::kBDAFactorsTable);
  Table spw = ms_.spectralWindow();

  unsigned int nbl = factors.nrow();
  unsigned int ncorr =
      ArrayColumn<Complex>(ms_, MS::columnName(MS::DATA)).shape(0)[0];

  // Required columns to read
  ScalarColumn<int> factor_col(factors, base::DP3MS::kFactor);
  ScalarColumn<double> time_col(ms_, MS::columnName(MS::TIME));
  ScalarColumn<int> ant1_col(factors, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(factors, MS::columnName(MS::ANTENNA2));
  ScalarColumn<int> ids_col(factors, base::DP3MS::kSpectralWindowId);
  ArrayColumn<double> freqs_col(spw, MS_SPW::columnName(MS_SPW::CHAN_FREQ));
  ArrayColumn<double> widths_col(spw, MS_SPW::columnName(MS_SPW::CHAN_WIDTH));

  // Fill info with the data required to repopulate BDA_FACTORS
  std::vector<std::vector<double>> freqs(nbl);
  std::vector<std::vector<double>> widths(nbl);
  std::vector<unsigned int> baseline_factors(nbl);
  std::size_t max_n_channels = 0;
  for (unsigned int i = 0; i < nbl; ++i) {
    unsigned int spw_id = ids_col(i);
    baseline_factors[i] = factor_col(i);
    freqs[i] = freqs_col.get(spw_id).tovector();
    widths[i] = widths_col.get(spw_id).tovector();
    max_n_channels = std::max(max_n_channels, freqs[i].size());

    desc_id_to_nchan_[ids_col(i)] = freqs[i].size();
    bl_to_id_[std::make_pair(ant1_col(i), ant2_col(i))] = i;
  }

  Table anttab(ms_.keywordSet().asTable(base::DP3MS::kAntennaTable));
  ScalarColumn<casacore::String> name_col(anttab,
                                          MS_Ant::columnName(MS_Ant::NAME));
  ScalarColumn<double> diam_col(anttab,
                                MS_Ant::columnName(MS_Ant::DISH_DIAMETER));
  ScalarMeasColumn<MPosition> ant_col(anttab,
                                      MS_Ant::columnName(MS_Ant::POSITION));
  vector<MPosition> antPos;
  antPos.reserve(anttab.nrow());
  for (unsigned int i = 0; i < anttab.nrow(); ++i) {
    antPos.push_back(ant_col(i));
  }

  const unsigned int start_channel = 0;
  const std::string antenna_set = base::ReadAntennaSet(ms_);
  info() = DPInfo(ncorr, max_n_channels, start_channel, antenna_set);

  // Set time info.

  // Get the first time and last time.
  // When all baselines are averaged in time, the first and start time will
  // be incorrect: Together with the (non-averaged) time interval they yield
  // an incorrect start time, for example. AST-1092 should address this issue.
  std::string time_query =
      "select gmin(TIME_CENTROID) as first_time, "
      "gmax(TIME_CENTROID) as last_time "
      "from $1";
  casacore::TaQLResult time_result = casacore::tableCommand(time_query, ms_);
  ScalarColumn<double> first_time_col(time_result.table(), "first_time");
  ScalarColumn<double> last_time_col(time_result.table(), "last_time");

  casacore::Table axis =
      ms_.keywordSet().asTable(base::DP3MS::kBDATimeAxisTable);
  const double interval = axis.col(base::DP3MS::kUnitTimeInterval).getDouble(0);
  is_interval_integer_ = axis.col(base::DP3MS::kIntervalFactors).getBool(0);

  info().setTimes(first_time_col(0), last_time_col(0), interval);

  // Set antenna/baseline info.
  casacore::Vector<casacore::String> names = name_col.getColumn();
  info().setAntennas(std::vector<std::string>(names.begin(), names.end()),
                     diam_col.getColumn().tovector(), antPos,
                     ant1_col.getColumn().tovector(),
                     ant2_col.getColumn().tovector());

  const MS_Field& field = ms_.field();
  if (field.nrow() != 1)
    throw std::runtime_error("Multiple entries in FIELD table");

  ArrayMeasColumn<MDirection> col_phase_dir(
      field, MS_Field::columnName(MS_Field::PHASE_DIR));
  ArrayMeasColumn<MDirection> col_delay_dir(
      field, MS_Field::columnName(MS_Field::DELAY_DIR));

  if (col_phase_dir(0).empty() || col_delay_dir(0).empty())
    throw std::runtime_error("PHASE_DIR, DELAY_DIR are undefined");

  const MDirection phase_center = *(col_phase_dir(0).data());
  const MDirection delay_center = *(col_delay_dir(0).data());

  MDirection tile_beam_dir;
  try {
    tile_beam_dir = everybeam::ReadTileBeamDirection(ms_);
  } catch (const std::runtime_error& error) {
    // everybeam throws an exception error if telescope != [LOFAR, AARTFAAC]
    // in that case, default back to "DELAY_DIR"
    tile_beam_dir = delay_center;
  }

  // Get the array position using the telescope name from the OBSERVATION
  // subtable.
  MS_Obs observation = ms_.observation();

  const ScalarColumn<casacore::String> col_telescope_name(
      observation, MS_Obs::columnName(MS_Obs::TELESCOPE_NAME));
  MPosition array_pos;
  if (observation.nrow() == 0 ||
      !MeasTable::Observatory(array_pos, col_telescope_name(0))) {
    // If not found, use the position of the middle antenna.
    array_pos = antPos[antPos.size() / 2];
  }

  info().setArrayInformation(array_pos, phase_center, delay_center,
                             tile_beam_dir);

  info().update(std::move(baseline_factors));
  info().setChannels(std::move(freqs), std::move(widths));

  const std::string kFlagColumnName = "";  // Reading flags is not supported.
  info().setMsNames(msName(), data_column_name_, kFlagColumnName,
                    weight_column_name_);
  info().setIsBDAIntervalFactorInteger(is_interval_integer_);
}

void MSBDAReader::show(std::ostream& os) const {
  os << "MSBDAReader\n";
  os << "  input MS:       " << msName() << '\n';
  if (ms_.isNull()) {
    os << "    *** MS does not exist ***\n";
  } else {
    os << "  band            " << getInfo().spectralWindow() << '\n';
    os << "  start_chan:     " << 0 << '\n';
    os << "  nchan:          " << getInfo().nchan() << '\n';
    os << "  ncorrelations:  " << getInfo().ncorr() << '\n';
    os << "  nbaselines:     " << getInfo().nbaselines() << '\n';
    os << "  first time:     " << MVTime::Format(MVTime::YMD)
       << MVTime(getInfo().firstTime() / (24 * 3600.)) << '\n';
    os << "  last time:      " << MVTime::Format(MVTime::YMD)
       << MVTime(getInfo().lastTime() / (24 * 3600.)) << '\n';
    os << "  ntimes:         " << getInfo().ntime() << '\n';
    os << "  time interval:  " << getInfo().timeInterval() << '\n';
    os << "  DATA column:    " << data_column_name_;
    os << '\n';
    os << "  WEIGHT column:  " << weight_column_name_ << '\n';
  }
}

void MSBDAReader::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MSBDAReader" << '\n';
}

}  // namespace steps
}  // namespace dp3
