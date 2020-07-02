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

#include <limits>
  
namespace DP3 {
  namespace DPPP {

    BDABuffer::Row::Row(double time,
                        double interval,
                        rownr_t row_nr,
                        std::size_t baseline_nr,
                        std::size_t n_channels,
                        std::size_t n_correlations,
                        const std::vector<std::complex<float>>::iterator data,
                        const std::vector<bool>::iterator flags,
                        const std::vector<float>::iterator weights,
                        const std::vector<bool>::iterator full_res_flags,
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

    BDABuffer::BDABuffer(const std::size_t pool_size)
    : data_(), flags_(), weights_(), full_res_flags_(), rows_()
    {
      data_.reserve(pool_size);
      flags_.reserve(pool_size);
      weights_.reserve(pool_size);
      full_res_flags_.reserve(pool_size);
    }

    BDABuffer::BDABuffer(const BDABuffer& other)
    : data_(other.data_)
    , flags_(other.flags_)
    , weights_(other.weights_)
    , full_res_flags_(other.full_res_flags_)
    , rows_()
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
      }
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
      if (!rows_.empty() && TimeIsLess(time, rows_.back().time_)) {
        throw std::invalid_argument("Rows are not ordered by start time");
      }
      
      const std::size_t kDataSize = n_channels * n_correlations;

      // Check if there is enough capacity left.
      if ((data_.capacity() - data_.size()) < kDataSize) {
        return false;
      }

      rows_.emplace_back(time,
                         interval,
                         row_nr,
                         baseline_nr,
                         n_channels,
                         n_correlations,
                         data_.end(),
                         flags_.end(),
                         weights_.end(),
                         full_res_flags_.end(),
                         uvw);
      if (data) {
        data_.insert(data_.end(), data, data + kDataSize);
      } else {
        const std::complex<float> nan(std::nanf(""), std::nanf(""));
        data_.insert(data_.end(), kDataSize, nan);
      }

      if (flags) {
        flags_.insert(flags_.end(), flags, flags + kDataSize);
      } else {
        flags_.insert(flags_.end(), kDataSize, false);
      }

      if (weights) {
        weights_.insert(weights_.end(), weights, weights + kDataSize);
      } else {
        const float nan = std::nanf("");
        weights_.insert(weights_.end(), kDataSize, nan);
      }

      if (full_res_flags) {
        full_res_flags_.insert(full_res_flags_.end(),
                               full_res_flags,
                               full_res_flags + kDataSize);
      } else {
        full_res_flags_.insert(full_res_flags_.end(), kDataSize, false);
      }

      return true;
    }

    void BDABuffer::Clear()
    {
      data_.clear();
      flags_.clear();
      weights_.clear();
      full_res_flags_.clear();
      rows_.clear();
    }
  }
}
