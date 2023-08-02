// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MsColumnReader.h"

#include <string>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableIter.h>

#include "../base/DPLogger.h"
using casacore::ArrayColumn;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

MsColumnReader::MsColumnReader(const common::ParameterSet& parset,
                               const std::string& prefix,
                               const std::string& column)
    : table_(),
      name_(prefix),
      column_name_(parset.getString(prefix + "column", column)) {}

bool MsColumnReader::process(std::unique_ptr<DPBuffer> buffer) {
  buffer->GetData().resize(
      {getInfo().nbaselines(), getInfo().nchan(), getInfo().ncorr()});
  const casacore::IPosition shape(3, getInfo().ncorr(), getInfo().nchan(),
                                  getInfo().nbaselines());
  casacore::Cube<casacore::Complex> data(shape, buffer->GetData().data(),
                                         casacore::SHARE);

  ArrayColumn<casacore::Complex> model_column(table_, column_name_);
  model_column.getColumnCells(buffer->GetRowNumbers(), data);

  getNextStep()->process(std::move(buffer));

  return false;
}

void MsColumnReader::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);

  table_ = casacore::Table(_info.msName());

  // Verify that the column being read has the same shape as the buffer which
  // will be given to this step's process function
  casacore::TableIterator table_iterator_single_time(
      table_, casacore::Block<casacore::String>(1, "TIME"),
      casacore::TableIterator::Ascending, casacore::TableIterator::NoSort);
  size_t n_baselines_column = table_iterator_single_time.table().nrow();

  ArrayColumn<casacore::Complex> model_column(table_, column_name_);

  if ((_info.ncorr() != model_column.shape(0)[0]) ||
      (_info.nchan() != model_column.shape(0)[1]) ||
      _info.nbaselines() != n_baselines_column || _info.metaChanged()) {
    throw std::runtime_error(
        "The column " + column_name_ +
        " has shape NCORR: " + std::to_string(model_column.shape(0)[0]) +
        ", NCHAN: " + std::to_string(model_column.shape(0)[1]) +
        ", NBL: " + std::to_string(n_baselines_column) +
        ", while the input buffer for the current step has shape NCORR: " +
        std::to_string(_info.ncorr()) +
        ", NCHAN: " + std::to_string(_info.nchan()) +
        ", NBL: " + std::to_string(_info.nbaselines()) +
        ".\nAny operation which alters the number of channels, correlations "
        "and baselines can cause a shape mismatch (example: filter, average, "
        "...). Also operations which will change the metadata can cause a "
        "mismatch (example: phase shift, upsample, station adder, ..) \n\n");
  }
}

void MsColumnReader::finish() { getNextStep()->finish(); }

void MsColumnReader::show(std::ostream& os) const {
  os << "MsColumnReader " << name_ << '\n';
  os << "  column:      " << column_name_ << '\n';
  os << '\n';
}

void MsColumnReader::showTimings(std::ostream& os,
                                 [[maybe_unused]] double duration) const {
  os << " MsColumnReader " << name_ << '\n';
}

base::Direction MsColumnReader::GetFirstDirection() const {
  return getInfo().phaseCenterDirection();
}

}  // namespace steps
}  // namespace dp3
