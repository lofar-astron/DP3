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

#include "BDAIntervalBuffer.h"

#include <boost/make_unique.hpp>

#include <algorithm>
#include <cassert>

namespace DP3 {
  namespace DPPP {

    BDAIntervalBuffer::BDAIntervalBuffer(const double time,
                                         const double interval)
    : time_(time)
    , interval_(interval)
    , buffers_()
    , rows_()
    {}

    void BDAIntervalBuffer::AddBuffer(const BDABuffer& buffer)
    {
      // Check if 'buffer' is valid.
      if (buffer.GetRows().empty()) {
        return;
      }

      if (!buffers_.empty() &&
          BDABuffer::TimeIsLess(buffer.GetRows().front().time_,
                                buffers_.back()->GetRows().back().time_)) {
        throw std::invalid_argument("New row does not follow existing row.");
      }

      buffers_.push_back(boost::make_unique<BDABuffer>(buffer));

      // Add the rows of the new buffer to valid_rows_.
      for(const auto& row: buffers_.back()->GetRows()) {
        rows_.push_back(&row);
      }
    }

    void BDAIntervalBuffer::Advance(const double interval)
    {
      time_ += interval_;
      interval_ = interval;
      removeOld();
    }

    bool BDAIntervalBuffer::IsComplete() const
    {
      if (buffers_.empty()) {
        return false;
      } else {
        // If the start time of the last input row is later than the end time of
        // the interval, the interval is complete.
        const double kEnd = time_ + interval_;
        const double kLastRowStart = buffers_.back()->GetRows().back().time_;
        return BDABuffer::TimeIsLess(kEnd, kLastRowStart);
      }
    }

    std::unique_ptr<BDABuffer>
    BDAIntervalBuffer::GetBuffer(const BDABuffer::Fields fields) const
    {
      // Count the number of elements in all rows.
      std::size_t pool_size = 0;
      for (const auto& row : rows_) {
        if (BDABuffer::TimeIsGreaterEqual(row->time_, time_ + interval_)) {
          break;
        } else {
          assert(BDABuffer::TimeIsLess(time_, row->time_ + row->interval_));
          pool_size += row->GetDataSize();
        }
      }

      // Create the result buffer and fill it.
      auto result = boost::make_unique<BDABuffer>(pool_size, fields);
      for (const auto& row : rows_) {
        if (BDABuffer::TimeIsGreaterEqual(row->time_, time_ + interval_)) {
          break;
        }

        const double kRowTime = std::max(row->time_, time_);
        const double kRowInterval = std::min(row->time_ + row->interval_,
                                             time_ + interval_) - kRowTime;
        const bool success = result->AddRow(kRowTime,
                                            kRowInterval,
                                            row->row_nr_,
                                            row->baseline_nr_,
                                            row->n_channels_,
                                            row->n_correlations_,
                                            row->data_,
                                            row->flags_,
                                            row->weights_,
                                            row->full_res_flags_,
                                            row->uvw_);
        (void)success;
        assert(success);

        // Adjust weights if the row partially overlaps the interval.
        if (fields[BDABuffer::Field::kWeights] && row->weights_ &&
            BDABuffer::TimeIsLess(kRowInterval, interval_) &&
            BDABuffer::TimeIsLess(kRowInterval, row->interval_)) {

          const double kWeightFactor = kRowInterval / row->interval_;
          float* weights = result->GetWeights(result->GetRows().size() - 1);
          const std::size_t kDataSize = row->GetDataSize();

          for (std::size_t i = 0; i < kDataSize; ++i) {
            *weights *= kWeightFactor;
            ++weights;
          }
        }
      }

      return result;
    }

    void BDAIntervalBuffer::removeOld()
    {
      // Remove old rows from rows_ before removing old buffers from buffers_,
      // since the elements of rows_ point to rows in buffers_.

      // If the start time is after time_, the loop is done, since the rows
      // are ordered by start time.
      auto row_it = rows_.begin();
      while(row_it != rows_.end() &&
            BDABuffer::TimeIsLess((*row_it)->time_, time_)) {
        const double kRowEnd = (*row_it)->time_ + (*row_it)->interval_;
        if (BDABuffer::TimeIsLess(time_, kRowEnd)) {
          ++row_it;
        } else {
          row_it = rows_.erase(row_it);
        }
      }

      // Only buffers with start times before time_ can possibly be removed.
      // If the start time is after time_, the loop is done, since the buffers
      // in buffers_ and the rows in those buffers are ordered by start time.
      auto buffer_it = buffers_.begin();
      while(buffer_it != buffers_.end() &&
            BDABuffer::TimeIsLess((*buffer_it)->GetRows().back().time_, time_ )) {

        // Check if there is a row with an end time after time_.
        bool buffer_is_old = true;
        for (const auto& row : (*buffer_it)->GetRows()) {
          const double kRowEnd = row.time_ + row.interval_;
          if (BDABuffer::TimeIsLess(time_, kRowEnd)) {
            buffer_is_old = false;
            break;
          }
        }

        if (buffer_is_old) {
          buffer_it = buffers_.erase(buffer_it);
        } else {
          ++buffer_it;
        }
      }
    }
  }
}
