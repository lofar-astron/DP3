// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BDABuffer.h"

#include <algorithm>
#include <limits>

namespace dp3 {
namespace base {

BDABuffer::Row::Row(double _time, double _interval, double _exposure,
                    common::rownr_t _row_nr, std::size_t _baseline_nr,
                    std::size_t _n_channels, std::size_t _n_correlations,
                    std::complex<float>* _data, bool* _flags, float* _weights,
                    bool* _full_res_flags, const double* const _uvw)
    : time(_time),
      interval(_interval),
      exposure(_exposure),
      row_nr(_row_nr),
      baseline_nr(_baseline_nr),
      n_channels(_n_channels),
      n_correlations(_n_correlations),
      data(_data),
      flags(_flags),
      weights(_weights),
      full_res_flags(_full_res_flags),
      uvw{_uvw ? _uvw[0] : std::numeric_limits<double>::quiet_NaN(),
          _uvw ? _uvw[1] : std::numeric_limits<double>::quiet_NaN(),
          _uvw ? _uvw[2] : std::numeric_limits<double>::quiet_NaN()} {}

BDABuffer::BDABuffer(const std::size_t pool_size, const Fields& fields)
    : data_(),
      flags_(),
      weights_(),
      full_res_flags_(),
      rows_(),
      original_capacity_(pool_size),
      remaining_capacity_(pool_size) {
  if (fields.data) {
    data_.reserve(remaining_capacity_);
  }
  if (fields.flags) {
    flags_.reserve(remaining_capacity_);
  }
  if (fields.weights) {
    weights_.reserve(remaining_capacity_);
  }
  if (fields.full_res_flags) {
    full_res_flags_.reserve(remaining_capacity_);
  }
}

// When copying the memory pools in this copy-constructor, the capacity
// of the new memory pools will equal their size. There is therefore
// no remaining capacity in the new copy.
BDABuffer::BDABuffer(const BDABuffer& other, const Fields& fields,
                     const Fields& copy_fields)
    : data_(),
      flags_(),
      weights_(),
      full_res_flags_(),
      rows_(),
      original_capacity_(other.original_capacity_ - other.remaining_capacity_),
      remaining_capacity_(0) {
  if (fields.data) {
    if (copy_fields.data) data_ = other.data_;
    data_.resize(original_capacity_);
  }
  if (fields.flags) {
    if (copy_fields.flags) flags_ = other.flags_;
    flags_.resize(original_capacity_);
  }
  if (fields.weights) {
    if (copy_fields.weights) weights_ = other.weights_;
    weights_.resize(original_capacity_);
  }
  if (fields.full_res_flags) {
    if (copy_fields.full_res_flags) full_res_flags_ = other.full_res_flags_;
    full_res_flags_.resize(original_capacity_);
  }
  CopyRows(other.rows_);
}

void BDABuffer::CopyRows(const std::vector<BDABuffer::Row>& existing_rows) {
  std::complex<float>* row_data = data_.empty() ? nullptr : data_.data();
  bool* row_flags = flags_.empty() ? nullptr : flags_.data();
  float* row_weights = weights_.empty() ? nullptr : weights_.data();
  bool* row_full_res_flags =
      full_res_flags_.empty() ? nullptr : full_res_flags_.data();

  // Note: 'existing_rows' can reference 'rows_' !
  std::vector<Row> new_rows;
  new_rows.reserve(existing_rows.size());
  for (const Row& row : existing_rows) {
    new_rows.emplace_back(row.time, row.interval, row.exposure, row.row_nr,
                          row.baseline_nr, row.n_channels, row.n_correlations,
                          row_data, row_flags, row_weights, row_full_res_flags,
                          row.uvw);
    const std::size_t row_size = row.GetDataSize();
    if (row_data) row_data += row_size;
    if (row_flags) row_flags += row_size;
    if (row_weights) row_weights += row_size;
    if (row_full_res_flags) row_full_res_flags += row_size;
  }
  rows_ = std::move(new_rows);
}

void BDABuffer::SetFields(const Fields& fields) {
  if ((fields.data == !data_.empty()) && (fields.flags == !flags_.empty()) &&
      (fields.weights == !weights_.empty()) &&
      (fields.full_res_flags == !full_res_flags_.empty())) {
    return;
  }

  if (fields.data) {
    data_.resize(original_capacity_);
  } else {
    data_.clear();
    data_.shrink_to_fit();
  }
  if (fields.flags) {
    flags_.resize(original_capacity_);
  } else {
    flags_.clear();
    flags_.shrink_to_fit();
  }
  if (fields.weights) {
    weights_.resize(original_capacity_);
  } else {
    weights_.clear();
    weights_.shrink_to_fit();
  }
  if (fields.full_res_flags) {
    full_res_flags_.resize(original_capacity_);
  } else {
    full_res_flags_.clear();
    full_res_flags_.shrink_to_fit();
  }

  CopyRows(rows_);
}

void BDABuffer::Clear() {
  data_.clear();
  flags_.clear();
  weights_.clear();
  full_res_flags_.clear();
  rows_.clear();
  remaining_capacity_ = original_capacity_;
}

std::size_t BDABuffer::GetNumberOfElements() const {
  return original_capacity_ - remaining_capacity_;
}

bool BDABuffer::AddRow(double time, double interval, double exposure,
                       std::size_t baseline_nr, std::size_t n_channels,
                       std::size_t n_correlations,
                       const std::complex<float>* const data,
                       const bool* const flags, const float* const weights,
                       const bool* const full_res_flags,
                       const double* const uvw) {
  if (!rows_.empty() &&
      TimeIsLessEqual(time + interval / 2,
                      rows_.back().time - rows_.back().interval / 2)) {
    throw std::invalid_argument("Rows are not properly ordered");
  }
  const std::size_t n_elements = n_channels * n_correlations;
  if (remaining_capacity_ < n_elements) {
    return false;
  }
  remaining_capacity_ -= n_elements;
  std::complex<float>* row_data = nullptr;
  bool* row_flags = nullptr;
  float* row_weights = nullptr;
  bool* row_full_res_flags = nullptr;
  if (data_.capacity() > 0) {
    row_data = data_.end();
    if (data) {
      data_.insert(data_.end(), data, data + n_elements);
    } else {
      const float kNaN = std::numeric_limits<float>::quiet_NaN();
      data_.insert(data_.end(), n_elements, {kNaN, kNaN});
    }
  }
  if (flags_.capacity() > 0) {
    row_flags = flags_.end();
    if (flags) {
      flags_.insert(flags_.end(), flags, flags + n_elements);
    } else {
      flags_.insert(flags_.end(), n_elements, false);
    }
  }
  if (weights_.capacity() > 0) {
    row_weights = weights_.end();
    if (weights) {
      weights_.insert(weights_.end(), weights, weights + n_elements);
    } else {
      const float kNaN = std::numeric_limits<float>::quiet_NaN();
      weights_.insert(weights_.end(), n_elements, kNaN);
    }
  }
  if (full_res_flags_.capacity() > 0) {
    row_full_res_flags = full_res_flags_.end();
    if (full_res_flags) {
      full_res_flags_.insert(full_res_flags_.end(), full_res_flags,
                             full_res_flags + n_elements);
    } else {
      full_res_flags_.insert(full_res_flags_.end(), n_elements, false);
    }
  }

  const common::rownr_t row_nr = rows_.empty() ? 0 : rows_.back().row_nr + 1;
  rows_.emplace_back(time, interval, exposure, row_nr, baseline_nr, n_channels,
                     n_correlations, row_data, row_flags, row_weights,
                     row_full_res_flags, uvw);
  return true;
}

void BDABuffer::SetBaseRowNr(common::rownr_t row_nr) {
  for (Row& row : rows_) {
    row.row_nr = row_nr;
    ++row_nr;
  }
}

bool BDABuffer::Row::IsMetadataEqual(const BDABuffer::Row& other) const {
  for (std::size_t i = 0; i < 3; ++i) {
    if (std::isnan(uvw[i])) {
      if (!std::isnan(other.uvw[i])) return false;
    } else {
      if (!BDABuffer::TimeIsEqual(uvw[i], other.uvw[i])) return false;
    }
  }

  return ((BDABuffer::TimeIsEqual(time, other.time)) &&
          (BDABuffer::TimeIsEqual(interval, other.interval)) &&
          (BDABuffer::TimeIsEqual(exposure, other.exposure)) &&
          (row_nr == other.row_nr) && (baseline_nr == other.baseline_nr) &&
          (n_channels == other.n_channels) &&
          (n_correlations == other.n_correlations));
}

bool BDABuffer::IsMetadataEqual(const BDABuffer& other) const {
  if (GetRows().size() == other.GetRows().size()) {
    auto this_row = GetRows().begin();
    auto other_row = other.GetRows().begin();
    while (this_row != GetRows().end()) {
      if (!this_row->IsMetadataEqual(*other_row)) {
        return false;
      }
      ++this_row;
      ++other_row;
    }
    return true;
  } else {
    return false;
  }
}

}  // namespace base
}  // namespace dp3
