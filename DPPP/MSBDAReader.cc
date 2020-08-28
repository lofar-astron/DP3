// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
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

#include "MSBDAReader.h"
#include "BDABuffer.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"
#include "BDAMS.h"

#include <EveryBeam/load.h>
#include <EveryBeam/lofarreadutils.h>

#include "../Common/ParameterSet.h"
#include "../Common/BaselineSelect.h"

#include <casacore/casa/OS/HostInfo.h>
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
#include <boost/make_unique.hpp>

using casacore::ArrayColumn;
using casacore::Complex;
using casacore::Cube;
using casacore::IPosition;
using casacore::MeasurementSet;
using casacore::MPosition;
using casacore::MS;
using casacore::MVTime;
using casacore::RefRows;
using casacore::ScalarColumn;
using casacore::ScalarMeasColumn;
using casacore::Table;
using casacore::TableLock;

using std::cout;
using std::endl;
using std::vector;

using namespace DP3::DPPP::BDAMS;

namespace DP3 {
namespace DPPP {

MSBDAReader::MSBDAReader()
    : read_vis_data_(false), last_ms_time_(0), nread_(0), max_chan_width_(0) {}

MSBDAReader::MSBDAReader(const string& msName, const ParameterSet& parset,
                         const string& prefix)
    : ms_name_(msName),
      read_vis_data_(false),
      last_ms_time_(0),
      nread_(0),
      max_chan_width_(0) {
  spw_ = parset.getInt(prefix + "band", -1);
  data_col_name_ =
      parset.getString(prefix + "data_column", MS::columnName(MS::DATA));
  weight_col_name_ = parset.getString(prefix + "weightcolumn",
                                      MS::columnName(MS::WEIGHT_SPECTRUM));
}

MSBDAReader::~MSBDAReader() {}

DPStep::MSType MSBDAReader::outputs() const { return BDA; };

void MSBDAReader::updateInfo(const DPInfo& dpInfo) {
  info().setNThreads(dpInfo.nThreads());

  if (!Table::isReadable(ms_name_)) {
    throw std::invalid_argument("No such MS: " + ms_name_);
  }

  ms_ = MeasurementSet(ms_name_, TableLock::AutoNoReadLocking);

  if (!ms_.keywordSet().isDefined(kBDAFactorsTable) ||
      ms_.keywordSet().asTable(kBDAFactorsTable).nrow() == 0) {
    throw std::domain_error(
        "Input MS does not contain BDA data. Table BDA_FACTORS is missing "
        "or not filled");
  }

  // Find the nr of corr, chan, and baseline.

  ncorr_ = IPosition(
      ArrayColumn<Complex>(ms_, MS::columnName(MS::DATA)).shape(0))[0];

  // Read the meta data tables and store required values.
  FillInfoMetaData();

  // Read the antenna set.
  Table obstab(ms_.keywordSet().asTable(kObservationTable));
  string antenna_set;
  if (obstab.nrow() > 0 && obstab.tableDesc().isColumn(kLofarAntennaSet)) {
    antenna_set = ScalarColumn<casacore::String>(obstab, kLofarAntennaSet)(0);
  }

  // Determine the pool size for the bda buffers
  const size_t rows_per_buffer = nbl_ * ncorr_ * max_chan_width_;
  pool_size_ = nbl_ * ncorr_ * max_chan_width_;

  double start_time = 0;
  unsigned int start_chan = 0;
  unsigned int ntime_approx = std::ceil(ms_.nrow() / (float)rows_per_buffer);

  info().init(ncorr_, start_chan, max_chan_width_, ntime_approx, start_time,
              interval_, msName(), antenna_set);
  info().setDataColName(data_col_name_);
  info().setWeightColName(weight_col_name_);
}

string MSBDAReader::msName() const { return ms_.tableName(); }

void MSBDAReader::setReadVisData(bool readVisData) {
  read_vis_data_ = readVisData || read_vis_data_;
}

bool MSBDAReader::process(const DPBuffer&) {
  NSTimer::StartStop sstime(timer_);
  std::unique_ptr<BDABuffer> buffer = boost::make_unique<BDABuffer>(pool_size_);

  ScalarColumn<int> ant1_col(ms_, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(ms_, MS::columnName(MS::ANTENNA2));
  ArrayColumn<Complex> data_col(ms_, MS::columnName(MS::DATA));
  ArrayColumn<float> weights_col(ms_, MS::columnName(MS::WEIGHT_SPECTRUM));
  ArrayColumn<double> uvw_col(ms_, MS::columnName(MS::UVW));
  ScalarColumn<double> interval_col(ms_, MS::columnName(MS::INTERVAL));
  ScalarColumn<double> exposure_col(ms_, MS::columnName(MS::EXPOSURE));
  ScalarColumn<int> data_desc_id_col(ms_, MS::columnName(MS::DATA_DESC_ID));

  // Cache the data that will be add to the buffer
  RefRows cell_range{nread_,
                     std::min(ms_.nrow() - 1, unsigned(nread_ + pool_size_))};
  auto data = data_col.getColumnCells(cell_range);
  auto weights = weights_col.getColumnCells(cell_range);
  auto uvw = uvw_col.getColumnCells(cell_range);
  auto interval = interval_col.getColumnCells(cell_range);
  auto exposure = exposure_col.getColumnCells(cell_range);
  auto data_desc_id = data_desc_id_col.getColumnCells(cell_range);

  unsigned i = 0;
  while (nread_ < ms_.nrow() && buffer->GetRemainingCapacity() > 0) {
    // Take time from row 0 in subset.
    double ms_time = ScalarColumn<double>(ms_, MS::columnName(MS::TIME))(0);

    if (ms_time < last_ms_time_) {
      DPLOG_WARN_STR("Time at rownr " + std::to_string(nread_) + " of MS " +
                     msName() + " is less than previous time slot");
      ++i;
      ++nread_;
      continue;
    }

    size_t blNr = bl_to_id_[std::make_pair(ant1_col(nread_), ant2_col(nread_))];

    buffer->AddRow(ms_time, interval[i], exposure[i], blNr,
                   desc_id_to_nchan_[data_desc_id[i]], info().ncorr(),
                   read_vis_data_ ? data[i].data() : nullptr, nullptr,
                   weights[i].data(), nullptr, uvw[i].data());

    last_ms_time_ = ms_time;
    ++i;
    ++nread_;
  }

  getNextStep()->process(std::move(buffer));

  // Return true while there are still items remaining
  return nread_ < ms_.nrow();
}

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::FillInfoMetaData() {
  Table factors = ms_.keywordSet().asTable(kBDAFactorsTable);
  Table axis = ms_.keywordSet().asTable(kBDATimeAxisTable);
  Table spw = ms_.keywordSet().asTable(kSpectralWindowTable);

  interval_ = axis.col(kUnitTimeInterval).getDouble(0);
  nbl_ = factors.nrow();

  // Required columns to read
  ScalarColumn<int> factor_col(factors, kFactor);
  ScalarColumn<int> ant1_col(factors, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(factors, MS::columnName(MS::ANTENNA2));
  ScalarColumn<int> ids_col(factors, kSpectralWindowId);
  ArrayColumn<double> freqs_col(spw, kChanFreq);
  ArrayColumn<double> widths_col(spw, kChanWidth);

  // Fill info with the data required to repopulate BDA_FACTORS
  std::vector<std::vector<double>> freqs(nbl_);
  std::vector<std::vector<double>> widths(nbl_);
  std::vector<unsigned int> baseline_factors(nbl_);
  for (unsigned int i = 0; i < nbl_; ++i) {
    unsigned int spw_id = ids_col(i);
    baseline_factors[i] = factor_col(i);
    freqs[i] = freqs_col.get(spw_id).tovector();
    widths[i] = widths_col.get(spw_id).tovector();

    max_chan_width_ = std::max(max_chan_width_, freqs[i].size());
    desc_id_to_nchan_[ids_col(i)] = freqs[i].size();
    bl_to_id_[std::make_pair(ant1_col(i), ant2_col(i))] = i;
  }

  Table anttab(ms_.keywordSet().asTable(kAntennaTable));
  ScalarColumn<casacore::String> name_col(anttab, kName);
  ScalarColumn<double> diam_col(anttab, kDishDiameter);
  ROScalarMeasColumn<MPosition> ant_col(anttab, kPosition);
  vector<MPosition> antPos;
  antPos.reserve(anttab.nrow());
  for (unsigned int i = 0; i < anttab.nrow(); ++i) {
    antPos.push_back(ant_col(i));
  }

  // Set antenna/baseline info.
  info().set(name_col.getColumn(), diam_col.getColumn(), antPos,
             ant1_col.getColumn(), ant2_col.getColumn());

  info().update(std::move(baseline_factors));
  info().set(std::move(freqs), std::move(widths));
}

void MSBDAReader::show(std::ostream& os) const {
  os << "MSReader\n";
  os << "  input MS:       " << ms_name_ << '\n';
  if (ms_.isNull()) {
    os << "    *** MS does not exist ***\n";
  } else {
    os << "  band            " << spw_ << '\n';
    os << "  start_chan:      " << 0 << '\n';
    os << "  nchan:          " << getInfo().nchan() << '\n';
    os << "  ncorrelations:  " << getInfo().ncorr() << '\n';
    os << "  nbaselines:     " << nbl_ << '\n';
    os << "  first time:     " << MVTime::Format(MVTime::YMD) << MVTime(0)
       << '\n';
    os << "  last time:      " << MVTime::Format(MVTime::YMD)
       << MVTime(last_ms_time_ / (24 * 3600.)) << '\n';
    os << "  ntimes:         " << getInfo().ntime() << '\n';
    os << "  time interval:  " << getInfo().timeInterval() << '\n';
    os << "  DATA column:    " << data_col_name_;
    os << '\n';
    os << "  WEIGHT column:  " << weight_col_name_ << '\n';
  }
}

void MSBDAReader::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MSBDAReader" << endl;
}

}  // namespace DPPP
}  // namespace DP3
