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

BDAIntervalBuffer::BDAIntervalBuffer(const double time, const double interval,
                                     const double max_row_interval)
    : completeness_(Completeness::kIncomplete),
      time_(time),
      interval_(interval),
      max_row_interval_(max_row_interval),
      buffers_(),
      current_rows_() {}

void BDAIntervalBuffer::AddBuffer(const BDABuffer& buffer) {
  if (!buffer.GetRows().empty()) {
    buffers_.emplace_back(buffer);

    for (const BDABuffer::Row& row : buffers_.back().GetRows()) {
      if (BDABuffer::TimeIsLess(max_row_interval_, row.interval)) {
        buffers_.pop_back();
        throw std::invalid_argument("Row interval > Max row interval");
      }
      // Add row to current_rows if it overlaps the current interval
      // (it isn't completely before or after the current interval).
      if (!BDABuffer::TimeIsLessEqual(row.time + row.interval, time_) &&
          !BDABuffer::TimeIsLessEqual(time_ + interval_, row.time)) {
        current_rows_.push_back(&row);
      }
    }

    if (completeness_ == Completeness::kIncomplete) {
      // Re-evaluate completeness in next IsComplete call.
      completeness_ = Completeness::kUnknown;
    }
  }
}

void BDAIntervalBuffer::Advance(const double interval) {
  completeness_ = Completeness::kUnknown;
  time_ += interval_;
  interval_ = interval;
  current_rows_.clear();

  // Fill current_rows with rows from the existing buffers.
  // Also detect and remove buffers that only have old rows.
  auto buffer_it = buffers_.begin();
  while (buffer_it != buffers_.end()) {
    bool buffer_is_old = true;
    for (const BDABuffer::Row& row : buffer_it->GetRows()) {
      // Check if the row interval overlaps the bda interval.
      // If row_end <= bda start, it does not overlap and is old -> Ignore.
      // If bda_end <= row start, it does not overlap and is new.
      // In all other cases, it overlaps and is added to current_rows_.
      if (!BDABuffer::TimeIsLessEqual(row.time + row.interval, time_)) {
        buffer_is_old = false;
        if (!BDABuffer::TimeIsLessEqual(time_ + interval_, row.time)) {
          current_rows_.push_back(&row);
        }
      }
    }
    if (buffer_is_old) {
      buffer_it = buffers_.erase(buffer_it);
    } else {
      ++buffer_it;
    }
  }
}

bool BDAIntervalBuffer::IsComplete() const {
  // Try using cached completeness status from earlier calls.
  if (Completeness::kUnknown != completeness_) {
    return Completeness::kComplete == completeness_;
  }
  const double time_complete = time_ + interval_ + max_row_interval_;
  // Evaluate the rows in the BDABuffers backwards.
  //
  // If a row starts at/after time_complete, the interval is complete:
  // Because of max_row_interval_, future rows start at/after kTimeEnd.
  //
  // If a row ends at/before time_complete, the interval is not complete:
  // All possible rows that satisfy the completeness criterion are after
  // that row but we didn't find any while searching backwards.
  //
  // When all rows are checked, the interval is also not complete.
  for (auto buffer_it = buffers_.rbegin();
       (buffer_it != buffers_.rend()) &&
       (Completeness::kUnknown == completeness_);
       ++buffer_it) {
    const std::vector<BDABuffer::Row>& rows = buffer_it->GetRows();
    for (auto row_it = rows.rbegin(); row_it != rows.rend(); ++row_it) {
      if (BDABuffer::TimeIsGreaterEqual(row_it->time, time_complete)) {
        completeness_ = Completeness::kComplete;
        break;
      } else if (BDABuffer::TimeIsLessEqual(row_it->time + row_it->interval,
                                            time_complete)) {
        completeness_ = Completeness::kIncomplete;
        break;
      }
    }
  }

  if (Completeness::kUnknown == completeness_) {
    completeness_ = Completeness::kIncomplete;
  }

  return Completeness::kComplete == completeness_;
}

std::unique_ptr<BDABuffer> BDAIntervalBuffer::GetBuffer(
    const BDABuffer::Fields& fields) const {
  // Count the number of elements in all current rows.
  std::size_t pool_size = 0;
  for (const BDABuffer::Row* row : current_rows_) {
    pool_size += row->n_elements;
  }

  // Create the result buffer and fill it.
  auto result = boost::make_unique<BDABuffer>(pool_size, fields);
  for (const BDABuffer::Row* row : current_rows_) {
    const double row_time = std::max(row->time, time_);
    const double row_interval =
        std::min(row->time + row->interval, time_ + interval_) - row_time;
    const bool success = result->AddRow(
        row_time, row_interval, row->row_nr, row->baseline_nr, row->n_elements,
        row->data, row->flags, row->weights, row->full_res_flags, row->uvw);
    (void)success;
    assert(success);

    // Adjust weights if the row partially overlaps the interval.
    if (fields.weights_ && row->weights &&
        BDABuffer::TimeIsLess(row_interval, interval_) &&
        BDABuffer::TimeIsLess(row_interval, row->interval)) {
      const double weight_factor = row_interval / row->interval;
      float* weights = result->GetRows().back().weights;
      for (std::size_t i = 0; i < row->n_elements; ++i) {
        *weights *= weight_factor;
        ++weights;
      }
    }
  }

  return result;
}

}  // namespace DPPP
}  // namespace DP3
