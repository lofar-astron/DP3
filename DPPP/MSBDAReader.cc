// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MSBDAReader.h"

#include "BDABuffer.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"
#include "MS.h"

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

namespace DP3 {
namespace DPPP {

MSBDAReader::MSBDAReader()
    : ms_(),
      ms_name_(),
      data_col_name_(),
      weight_col_name_(),
      read_vis_data_(false),
      last_ms_time_(0),
      interval_(0),
      spw_(0),
      nread_(0),
      timer_(),
      pool_size_(0),
      desc_id_to_nchan_(),
      bl_to_id_() {}

MSBDAReader::MSBDAReader(const std::string& msName, const ParameterSet& parset,
                         const std::string& prefix)
    : ms_(),
      ms_name_(msName),
      data_col_name_(
          parset.getString(prefix + "data_column", MS::columnName(MS::DATA))),
      weight_col_name_(parset.getString(prefix + "weightcolumn",
                                        MS::columnName(MS::WEIGHT_SPECTRUM))),
      read_vis_data_(false),
      last_ms_time_(0),
      interval_(0),
      spw_(parset.getInt(prefix + "band", -1)),
      nread_(0),
      timer_(),
      pool_size_(0),
      desc_id_to_nchan_(),
      bl_to_id_() {}

MSBDAReader::~MSBDAReader() {}

DPStep::MSType MSBDAReader::outputs() const { return BDA; };

void MSBDAReader::updateInfo(const DPInfo& dpInfo) {
  DPInput::updateInfo(dpInfo);

  if (!Table::isReadable(ms_name_)) {
    throw std::invalid_argument("No such MS: " + ms_name_);
  }

  ms_ = MeasurementSet(ms_name_, TableLock::AutoNoReadLocking);

  if (!ms_.keywordSet().isDefined(DP3MS::kBDAFactorsTable) ||
      ms_.keywordSet().asTable(DP3MS::kBDAFactorsTable).nrow() == 0) {
    throw std::domain_error(
        "Input MS does not contain BDA data. Table BDA_FACTORS is missing "
        "or not filled");
  }

  // Find the nr of correlations, channels, and baselines.

  unsigned int ncorr =
      ArrayColumn<Complex>(ms_, MS::columnName(MS::DATA)).shape(0)[0];

  // Read the meta data tables and store required values.
  FillInfoMetaData();

  // Read the antenna set.
  Table obstab(ms_.keywordSet().asTable(DP3MS::kObservationTable));
  std::string antenna_set;
  if (obstab.nrow() > 0 &&
      obstab.tableDesc().isColumn(DP3MS::kLofarAntennaSet)) {
    antenna_set =
        ScalarColumn<casacore::String>(obstab, DP3MS::kLofarAntennaSet)(0);
  }

  // TODO: Read actual values from the metadata.
  unsigned int start_chan = 0;
  // this ntime must only be used for DP3 progress
  unsigned int ntime = std::ceil(ms_.nrow() / (float)info().nbaselines());
  double start_time = 0;

  // FillInfoMetaData already set the number of channels via DPInfo::set.
  info().init(ncorr, start_chan, info().nchan(), ntime, start_time, interval_,
              msName(), antenna_set);
  info().setDataColName(data_col_name_);
  info().setWeightColName(weight_col_name_);
}

std::string MSBDAReader::msName() const { return ms_.tableName(); }

void MSBDAReader::setReadVisData(bool readVisData) {
  read_vis_data_ = readVisData || read_vis_data_;
}

bool MSBDAReader::process(const DPBuffer&) {
  NSTimer::StartStop sstime(timer_);

  // TODO: Pre-calculate actual required pool size beforehand.
  auto buffer = boost::make_unique<BDABuffer>(info().nbaselines() *
                                              info().nchan() * info().ncorr());

  ScalarColumn<int> ant1_col(ms_, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(ms_, MS::columnName(MS::ANTENNA2));
  ArrayColumn<float> weights_col(ms_, MS::columnName(MS::WEIGHT_SPECTRUM));
  ArrayColumn<casacore::Complex> data_col(ms_, MS::columnName(MS::DATA));
  ArrayColumn<double> uvw_col(ms_, MS::columnName(MS::UVW));
  ScalarColumn<double> time_col(ms_, MS::columnName(MS::TIME));
  ScalarColumn<double> interval_col(ms_, MS::columnName(MS::INTERVAL));
  ScalarColumn<double> exposure_col(ms_, MS::columnName(MS::EXPOSURE));
  ScalarColumn<int> data_desc_id_col(ms_, MS::columnName(MS::DATA_DESC_ID));

  // Cache the data that will be add to the buffer
  RefRows cell_range{
      nread_, std::min(ms_.nrow() - 1, rownr_t(nread_ + info().nbaselines()))};

  auto time = time_col.getColumnCells(cell_range);
  auto interval = interval_col.getColumnCells(cell_range);
  auto exposure = exposure_col.getColumnCells(cell_range);
  auto data_desc_id = data_desc_id_col.getColumnCells(cell_range);

  unsigned i = 0;
  while (nread_ < ms_.nrow() && i < info().nbaselines()) {
    const double ms_time = time[i];

    if (ms_time < last_ms_time_) {
      DPLOG_WARN_STR("Time at rownr " + std::to_string(nread_) + " of MS " +
                     msName() + " is less than previous time slot");
      ++i;
      ++nread_;
      continue;
    }

    const auto ant12 = std::make_pair(ant1_col(nread_), ant2_col(nread_));
    casacore::Complex* data_ptr = nullptr;
    std::vector<casacore::Complex> data;
    if (read_vis_data_) {
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
    ++i;
    ++nread_;
  }

  getNextStep()->process(std::move(buffer));

  // Return true while there are still items remaining
  return nread_ < ms_.nrow();
}

void MSBDAReader::finish() { getNextStep()->finish(); }

void MSBDAReader::FillInfoMetaData() {
  Table factors = ms_.keywordSet().asTable(DP3MS::kBDAFactorsTable);
  Table axis = ms_.keywordSet().asTable(DP3MS::kBDATimeAxisTable);
  Table spw = ms_.keywordSet().asTable(DP3MS::kSpectralWindowTable);

  interval_ = axis.col(DP3MS::kUnitTimeInterval).getDouble(0);
  unsigned int nbl = factors.nrow();

  // Required columns to read
  ScalarColumn<int> factor_col(factors, DP3MS::kFactor);
  ScalarColumn<int> ant1_col(factors, MS::columnName(MS::ANTENNA1));
  ScalarColumn<int> ant2_col(factors, MS::columnName(MS::ANTENNA2));
  ScalarColumn<int> ids_col(factors, DP3MS::kSpectralWindowId);
  ArrayColumn<double> freqs_col(spw, DP3MS::kChanFreq);
  ArrayColumn<double> widths_col(spw, DP3MS::kChanWidth);

  // Fill info with the data required to repopulate BDA_FACTORS
  std::vector<std::vector<double>> freqs(nbl);
  std::vector<std::vector<double>> widths(nbl);
  std::vector<unsigned int> baseline_factors(nbl);
  for (unsigned int i = 0; i < nbl; ++i) {
    unsigned int spw_id = ids_col(i);
    baseline_factors[i] = factor_col(i);
    freqs[i] = freqs_col.get(spw_id).tovector();
    widths[i] = widths_col.get(spw_id).tovector();

    desc_id_to_nchan_[ids_col(i)] = freqs[i].size();
    bl_to_id_[std::make_pair(ant1_col(i), ant2_col(i))] = i;
  }

  Table anttab(ms_.keywordSet().asTable(DP3MS::kAntennaTable));
  ScalarColumn<casacore::String> name_col(anttab, DP3MS::kName);
  ScalarColumn<double> diam_col(anttab, DP3MS::kDishDiameter);
  ROScalarMeasColumn<MPosition> ant_col(anttab, DP3MS::kPosition);
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
    os << "  nbaselines:     " << getInfo().nbaselines() << '\n';
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
