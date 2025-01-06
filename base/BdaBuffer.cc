// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <dp3/base/BdaBuffer.h>

#include <algorithm>
#include <cassert>
#include <limits>

namespace dp3 {
namespace base {

BdaBuffer::Row::Row(double _time, double _interval, double _exposure,
                    common::rownr_t _row_nr, std::size_t _baseline_nr,
                    std::size_t _n_channels, std::size_t _n_correlations,
                    std::size_t _offset, const double* const _uvw)
    : time(_time),
      interval(_interval),
      exposure(_exposure),
      row_nr(_row_nr),
      baseline_nr(_baseline_nr),
      n_channels(_n_channels),
      n_correlations(_n_correlations),
      offset(_offset),
      uvw{_uvw ? _uvw[0] : std::numeric_limits<double>::quiet_NaN(),
          _uvw ? _uvw[1] : std::numeric_limits<double>::quiet_NaN(),
          _uvw ? _uvw[2] : std::numeric_limits<double>::quiet_NaN()} {}

BdaBuffer::BdaBuffer(const std::size_t pool_size, const Fields& fields)
    : data_(),
      flags_(),
      weights_(),
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
}

// When copying the memory pools in this copy-constructor, the capacity
// of the new memory pools will equal their size. There is therefore
// no remaining capacity in the new copy.
BdaBuffer::BdaBuffer(const BdaBuffer& other, const Fields& fields,
                     const Fields& copy_fields)
    : data_(),
      flags_(),
      weights_(),
      rows_(other.rows_),
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
}

void BdaBuffer::Clear() {
  data_.clear();
  flags_.clear();
  weights_.clear();
  rows_.clear();
  remaining_capacity_ = original_capacity_;
}

bool BdaBuffer::AddRow(double time, double interval, double exposure,
                       std::size_t baseline_nr, std::size_t n_channels,
                       std::size_t n_correlations,
                       const std::complex<float>* const data,
                       const bool* const flags, const float* const weights,
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
  std::size_t offset = 0;
  if (data_.capacity() > 0) {
    offset = data_.size();
    if (data) {
      data_.insert(data_.end(), data, data + n_elements);
    } else {
      const float kNaN = std::numeric_limits<float>::quiet_NaN();
      data_.insert(data_.end(), n_elements, {kNaN, kNaN});
    }
  }
  if (flags_.capacity() > 0) {
    if (offset == 0) {
      offset = flags_.size();
    } else {
      assert(offset == flags_.size());
    }
    if (flags) {
      flags_.insert(flags_.end(), flags, flags + n_elements);
    } else {
      flags_.insert(flags_.end(), n_elements, false);
    }
  }
  if (weights_.capacity() > 0) {
    if (offset == 0) {
      offset = weights_.size();
    } else {
      assert(offset == weights_.size());
    }
    if (weights) {
      weights_.insert(weights_.end(), weights, weights + n_elements);
    } else {
      const float kNaN = std::numeric_limits<float>::quiet_NaN();
      weights_.insert(weights_.end(), n_elements, kNaN);
    }
  }

  const common::rownr_t row_nr = rows_.empty() ? 0 : rows_.back().row_nr + 1;
  rows_.emplace_back(time, interval, exposure, row_nr, baseline_nr, n_channels,
                     n_correlations, offset, uvw);
  return true;
}

void BdaBuffer::SetBaseRowNr(common::rownr_t row_nr) {
  for (Row& row : rows_) {
    row.row_nr = row_nr;
    ++row_nr;
  }
}

bool BdaBuffer::Row::IsMetadataEqual(const BdaBuffer::Row& other) const {
  for (std::size_t i = 0; i < 3; ++i) {
    if (std::isnan(uvw[i])) {
      if (!std::isnan(other.uvw[i])) return false;
    } else {
      if (!BdaBuffer::TimeIsEqual(uvw[i], other.uvw[i])) return false;
    }
  }

  return ((BdaBuffer::TimeIsEqual(time, other.time)) &&
          (BdaBuffer::TimeIsEqual(interval, other.interval)) &&
          (BdaBuffer::TimeIsEqual(exposure, other.exposure)) &&
          (row_nr == other.row_nr) && (baseline_nr == other.baseline_nr) &&
          (n_channels == other.n_channels) &&
          (n_correlations == other.n_correlations));
}

bool BdaBuffer::IsMetadataEqual(const BdaBuffer& other) const {
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
