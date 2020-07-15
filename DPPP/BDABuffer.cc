// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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

#include "BDABuffer.h"

#include <algorithm>
#include <limits>

namespace DP3 {
  namespace DPPP {

    BDABuffer::Row::Row(double time,
                        double interval,
                        rownr_t row_nr,
                        std::size_t baseline_nr,
                        std::size_t n_channels,
                        std::size_t n_correlations,
                        std::complex<float>* data,
                        bool* flags,
                        float* weights,
                        bool* full_res_flags,
                        const double *const uvw)
    : time_(time)
    , interval_(interval)
    , row_nr_(row_nr)
    , baseline_nr_(baseline_nr)
    , n_channels_(n_channels)
    , n_correlations_(n_correlations)
    , data_(data)
    , flags_(flags)
    , weights_(weights)
    , full_res_flags_(full_res_flags)
    , uvw_{ uvw ? uvw[0] : std::numeric_limits<double>::quiet_NaN(),
            uvw ? uvw[1] : std::numeric_limits<double>::quiet_NaN(),
            uvw ? uvw[2] : std::numeric_limits<double>::quiet_NaN() }
    {}

    BDABuffer::BDABuffer(const std::size_t pool_size,
                         const Fields& fields)
    : data_()
    , flags_()
    , weights_()
    , full_res_flags_()
    , rows_()
    , original_capacity_(pool_size)
    , remaining_capacity_(pool_size)
    {
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
    : data_(other.data_)
    , flags_(other.flags_)
    , weights_(other.weights_)
    , full_res_flags_(other.full_res_flags_)
    , rows_()
    , original_capacity_(0)
    , remaining_capacity_(0)
    {
      // Copy rows but set their iterators to the new memory pools.
      auto data_it = data_.begin();
      auto flags_it = flags_.begin();
      auto weights_it = weights_.begin();
      auto full_res_flags_it = full_res_flags_.begin();

      rows_.reserve(other.rows_.size());
      for (const auto& row : other.rows_) {
        rows_.emplace_back(row.time_,
                           row.interval_,
                           row.row_nr_,
                           row.baseline_nr_,
                           row.n_channels_,
                           row.n_correlations_,
                           data_it,
                           flags_it,
                           weights_it,
                           full_res_flags_it,
                           row.uvw_);

        const std::size_t kDataSize = row.n_channels_ * row.n_correlations_;
        data_it += kDataSize;
        flags_it += kDataSize;
        weights_it += kDataSize;
        full_res_flags_it += kDataSize;
        original_capacity_ += kDataSize;
      }
    }

    void BDABuffer::Clear()
    {
      data_.clear();
      flags_.clear();
      weights_.clear();
      full_res_flags_.clear();
      rows_.clear();
      remaining_capacity_ = original_capacity_;
    }

    std::size_t BDABuffer::GetNumberOfElements() const
    {
      return original_capacity_ - remaining_capacity_;
    }

    bool BDABuffer::AddRow(double time,
                           double interval,
                           rownr_t row_nr,
                           std::size_t baseline_nr,
                           std::size_t n_channels,
                           std::size_t n_correlations,
                           const std::complex<float>* const data,
                           const bool* const flags,
                           const float* const weights,
                           const bool* const full_res_flags,
                           const double* const uvw)
    {
      if (!rows_.empty() && TimeIsLessEqual(time + interval, rows_.back().time_)) {
        throw std::invalid_argument("Rows are not properly ordered");
      }

      const std::size_t kDataSize = n_channels * n_correlations;

      if (remaining_capacity_ < kDataSize) {
        return false;
      }
      remaining_capacity_ -= kDataSize;

      std::complex<float>* row_data = nullptr;
      bool* row_flags = nullptr;
      float* row_weights = nullptr;
      bool* row_full_res_flags = nullptr;

      if (data_.capacity() > 0) {
        row_data = data_.end();
        if (data) {
          data_.insert(data_.end(), data, data + kDataSize);
        } else {
          const float kNaN = std::numeric_limits<float>::quiet_NaN();
          data_.insert(data_.end(), kDataSize, { kNaN, kNaN });
        }
      }

      if (flags_.capacity() > 0) {
        row_flags = flags_.end();
        if (flags) {
          flags_.insert(flags_.end(), flags, flags + kDataSize);
        } else {
          flags_.insert(flags_.end(), kDataSize, false);
        }
      }

      if (weights_.capacity() > 0) {
        row_weights = weights_.end();
        if (weights) {
          weights_.insert(weights_.end(), weights, weights + kDataSize);
        } else {
          const float kNaN = std::numeric_limits<float>::quiet_NaN();
          weights_.insert(weights_.end(), kDataSize, kNaN);
        }
      }

      if (full_res_flags_.capacity() > 0) {
        row_full_res_flags = full_res_flags_.end();
        if (full_res_flags) {
          full_res_flags_.insert(full_res_flags_.end(),
                                 full_res_flags,
                                 full_res_flags + kDataSize);
        } else {
          full_res_flags_.insert(full_res_flags_.end(), kDataSize, false);
        }
      }

      rows_.emplace_back(time,
                         interval,
                         row_nr,
                         baseline_nr,
                         n_channels,
                         n_correlations,
                         row_data,
                         row_flags,
                         row_weights,
                         row_full_res_flags,
                         uvw);

      return true;
    }
  }
}
