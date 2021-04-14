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
  if (fields.data_) {
    data_.reserve(remaining_capacity_);
  }
  if (fields.flags_) {
    flags_.reserve(remaining_capacity_);
  }
  if (fields.weights_) {
    weights_.reserve(remaining_capacity_);
  }
  if (fields.full_res_flags_) {
    full_res_flags_.reserve(remaining_capacity_);
  }
}

// When copying the memory pools in this copy-constructor, the capacity
// of the new memory pools will equal their size. There is therefore
// no remaining capacity in the new copy.
BDABuffer::BDABuffer(const BDABuffer& other)
    : data_(other.data_),
      flags_(other.flags_),
      weights_(other.weights_),
      full_res_flags_(other.full_res_flags_),
      rows_(),
      original_capacity_(other.original_capacity_ - other.remaining_capacity_),
      remaining_capacity_(0) {
  // Copy rows but set their pointers to the new memory pools.
  rows_.reserve(other.rows_.size());
  for (const BDABuffer::Row& row : other.rows_) {
    rows_.emplace_back(row.time, row.interval, row.exposure, row.row_nr,
                       row.baseline_nr, row.n_channels, row.n_correlations,
                       data_.data() + (row.data - other.data_.data()),
                       flags_.data() + (row.flags - other.flags_.data()),
                       weights_.data() + (row.weights - other.weights_.data()),
                       full_res_flags_.data() +
                           (row.full_res_flags - other.full_res_flags_.data()),
                       row.uvw);
  }
}

BDABuffer::BDABuffer(const BDABuffer& other, const Fields& fields)
    : data_(),
      flags_(),
      weights_(),
      full_res_flags_(),
      rows_(),
      original_capacity_(other.original_capacity_ - other.remaining_capacity_),
      remaining_capacity_(0) {
  if (fields.data_) data_ = other.data_;
  if (fields.flags_) flags_ = other.flags_;
  if (fields.weights_) weights_ = other.weights_;
  if (fields.full_res_flags_) full_res_flags_ = other.full_res_flags_;

  // Copy rows but set their pointers to the new memory pools.
  rows_.reserve(other.rows_.size());
  for (const BDABuffer::Row& other_row : other.rows_) {
    std::complex<float>* row_data = nullptr;
    bool* row_flags = nullptr;
    float* row_weights = nullptr;
    bool* row_full_res_flags = nullptr;

    if (fields.data_ && !data_.empty())
      row_data = data_.data() + (other_row.data - other.data_.data());
    if (fields.flags_ && !flags_.empty())
      row_flags = flags_.data() + (other_row.flags - other.flags_.data());
    if (fields.weights_ && !weights_.empty())
      row_weights =
          weights_.data() + (other_row.weights - other.weights_.data());
    if (fields.full_res_flags_ && !full_res_flags_.empty())
      row_full_res_flags =
          full_res_flags_.data() +
          (other_row.full_res_flags - other.full_res_flags_.data());

    rows_.emplace_back(other_row.time, other_row.interval, other_row.exposure,
                       other_row.row_nr, other_row.baseline_nr,
                       other_row.n_channels, other_row.n_correlations, row_data,
                       row_flags, row_weights, row_full_res_flags,
                       other_row.uvw);
  }
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
  if (!rows_.empty() && TimeIsLessEqual(time + interval, rows_.back().time)) {
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

}  // namespace base
}  // namespace dp3
