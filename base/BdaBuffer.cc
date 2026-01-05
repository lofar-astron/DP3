// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>
#include <cassert>
#include <limits>

#include "base/BdaBuffer.h"

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

BdaBuffer::BdaBuffer(const std::size_t pool_size, const common::Fields& fields)
    : data_(),
      flags_(),
      weights_(),
      rows_(),
      original_capacity_(pool_size),
      remaining_capacity_(pool_size) {
  if (fields.Data()) {
    data_[""].reserve(remaining_capacity_);
  }
  if (fields.Flags()) {
    flags_.reserve(remaining_capacity_);
  }
  if (fields.Weights()) {
    weights_.reserve(remaining_capacity_);
  }
}

// When copying the memory pools in this copy-constructor, the capacity
// of the new memory pools will equal their size. There is therefore
// no remaining capacity in the new copy.
BdaBuffer::BdaBuffer(const BdaBuffer& other, const common::Fields& fields)
    : data_(),
      flags_(),
      weights_(),
      rows_(other.rows_),
      original_capacity_(other.original_capacity_ - other.remaining_capacity_),
      remaining_capacity_(0) {
  // aocommon::UVector ensures there is no remaining capacity in a vector copy.
  // The assertions below check that aocommon::UVector still has this behavior.
  if (fields.Data()) {
    data_ = other.data_;  // Copy all data buffers.
    data_[""];            // Add main data buffer if absent.
    // Resize all data buffers.
    for (auto it = data_.begin(); it != data_.end(); ++it) {
      it->second.resize(original_capacity_);
      assert(it->second.capacity() == it->second.size());
    }
  }
  if (fields.Flags()) {
    flags_ = other.flags_;
    flags_.resize(original_capacity_);
    assert(flags_.capacity() == flags_.size());
  }
  if (fields.Weights()) {
    weights_ = other.weights_;
    weights_.resize(original_capacity_);
    assert(weights_.capacity() == weights_.size());
  }
}

void BdaBuffer::MoveData(BdaBuffer& source, const std::string& source_name,
                         const std::string& target_name) {
  if (&source == this && source_name == target_name) return;

  // This code also handles internal moves, when &source == this.
  auto source_iterator = source.data_.find(source_name);
  assert(source_iterator != source.data_.end());
  assert(source_iterator->second.size() == GetNumberOfElements());
  data_.insert_or_assign(target_name, std::move(source_iterator->second));
  source.RemoveData(source_name);
}

void BdaBuffer::Clear() {
  for (auto data_it = data_.begin(); data_it != data_.end(); ++data_it) {
    data_it->second.clear();
  }
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
  for (auto data_it = data_.begin(); data_it != data_.end(); ++data_it) {
    offset = data_it->second.size();
    aocommon::UVector<std::complex<float>>& data_vector = data_it->second;
    if (data_it->first == "" && data) {
      data_vector.insert(data_vector.end(), data, data + n_elements);
    } else {
      data_vector.insert(data_vector.end(), n_elements, {0.0, 0.0});
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
      weights_.insert(weights_.end(), n_elements, 0.0);
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

const std::complex<float>* BdaBuffer::GetData(const std::string& name) const {
  const auto found = data_.find(name);
  if (found != data_.end() && !found->second.empty()) {
    return found->second.data();
  } else {
    return nullptr;
  }
}
std::complex<float>* BdaBuffer::GetData(const std::string& name) {
  const auto found = data_.find(name);
  if (found != data_.end() && !found->second.empty()) {
    return found->second.data();
  } else {
    return nullptr;
  }
}

const std::complex<float>* BdaBuffer::GetData(std::size_t row,
                                              const std::string& name) const {
  const std::complex<float>* data = GetData(name);
  if (data) data += rows_[row].offset;
  return data;
}
std::complex<float>* BdaBuffer::GetData(std::size_t row,
                                        const std::string& name) {
  std::complex<float>* data = GetData(name);
  if (data) data += rows_[row].offset;
  return data;
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
