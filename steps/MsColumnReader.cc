// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MsColumnReader.h"

#include <string>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableIter.h>

using casacore::ArrayColumn;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3::steps {

MsColumnReader::MsColumnReader(const common::ParameterSet& parset,
                               const std::string& prefix, MsType input_ms_type,
                               const std::string& column)
    : table_(),
      name_(prefix),
      column_name_(parset.getString(prefix + "column", column)),
      input_ms_type_(input_ms_type) {}

bool MsColumnReader::process(std::unique_ptr<DPBuffer> buffer) {
  buffer->GetData().resize(
      {getInfoOut().nbaselines(), getInfoOut().nchan(), getInfoOut().ncorr()});
  const casacore::IPosition shape(3, getInfoOut().ncorr(), getInfoOut().nchan(),
                                  getInfoOut().nbaselines());
  casacore::Cube<casacore::Complex> data(shape, buffer->GetData().data(),
                                         casacore::SHARE);

  assert(buffer->GetRowNumbers().size() == buffer->GetData().shape(0));

  model_column_.getColumnCells(buffer->GetRowNumbers(), data);

  getNextStep()->process(std::move(buffer));

  return false;
}

bool MsColumnReader::process(std::unique_ptr<base::BdaBuffer> buffer) {
  std::vector<base::BdaBuffer::Row>& rows = buffer->GetRows();
  buffer->AddData();
  std::complex<float>* data = buffer->GetData();
  for (base::BdaBuffer::Row& row : rows) {
    model_column_.get(row.row_nr, row_buffer_, true);
    const size_t n_correlations = row_buffer_.shape()[0];
    if (n_correlations != 4)
      throw std::runtime_error(
          "Column " + column_name_ + " reading error: Row " +
          std::to_string(row.row_nr) + " of the measurement set has " +
          std::to_string(n_correlations) +
          " correlations. Can't process measurement sets with n_correlations "
          "!= 4.");
    const size_t n_channels = row_buffer_.shape()[1];
    if (n_channels != row.n_channels)
      throw std::runtime_error(
          "Row " + std::to_string(row.row_nr) +
          " of the measurement set has an inconsistent number of " +
          std::to_string(n_channels) + " channels in column " + column_name_ +
          ", which doesn't match the input data with " +
          std::to_string(row.n_channels) + " channels.");

    std::copy(row_buffer_.cbegin(), row_buffer_.cend(), &data[row.offset]);
  }
  getNextStep()->process(std::move(buffer));
  return false;
}

void MsColumnReader::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);

  table_ = casacore::Table(_info.msName());
  model_column_ =
      casacore::ArrayColumn<std::complex<float>>(table_, column_name_);

  if (input_ms_type_ == MsType::kRegular) {
    // Verify that the column being read has the same shape as the buffer which
    // will be given to this step's process function
    casacore::TableIterator table_iterator_single_time(
        table_, casacore::Block<casacore::String>(1, "TIME"),
        casacore::TableIterator::Ascending, casacore::TableIterator::NoSort);
    size_t n_baselines_column = table_iterator_single_time.table().nrow();

    if ((_info.ncorr() != model_column_.shape(0)[0]) ||
        (_info.nchan() != model_column_.shape(0)[1]) ||
        _info.nbaselines() > n_baselines_column || _info.metaChanged()) {
      throw std::runtime_error(
          "The column " + column_name_ +
          " has shape NCORR: " + std::to_string(model_column_.shape(0)[0]) +
          ", NCHAN: " + std::to_string(model_column_.shape(0)[1]) +
          ", NBL: " + std::to_string(n_baselines_column) +
          ", while the input buffer for the current step has shape NCORR: " +
          std::to_string(_info.ncorr()) +
          ", NCHAN: " + std::to_string(_info.nchan()) +
          ", NBL: " + std::to_string(_info.nbaselines()) +
          ".\nAny operation which alters the number of channels, correlations "
          "can cause a shape mismatch (example: filter, average, ...). "
          "Also operations which will change the metadata can cause a "
          "mismatch (example: phase shift, upsample, station adder, ..)\n\n");
    }
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
  return getInfoOut().phaseCenterDirection();
}

}  // namespace dp3::steps
