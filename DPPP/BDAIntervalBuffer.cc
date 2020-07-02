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

    BDAIntervalBuffer::BDAIntervalBuffer(const std::size_t nBaselines,
                                         const double time,
                                         const double interval)
    : time_(time)
    , interval_(interval)
    , buffers_()
    , baselines_(nBaselines)
    {}
    
    void BDAIntervalBuffer::AddBuffer(const BDABuffer& buffer)
    {
      // Check if 'buffer' is valid.
      if (buffer.GetRows().empty()) {
        return;
      }
      
      for(const auto& row : buffer.GetRows()) {
        if (row.baseline_nr_ >= baselines_.size()) {
          throw std::invalid_argument("Invalid baseline");
        }
      }
      
      if (!buffers_.empty() &&
          BDABuffer::TimeIsLess(buffer.GetRows().front().time_,
                                buffers_.back()->GetRows().back().time_)) {
        throw std::invalid_argument("New row does not follow existing row.");
      }
              
      buffers_.push_back(boost::make_unique<BDABuffer>(buffer));
      
      // Add references in baselines_ to the rows of the new buffer.
      for(auto rowIt = buffers_.back()->GetRows().begin();
          rowIt != buffers_.back()->GetRows().end();
          ++rowIt) {
        baselines_[rowIt->baseline_nr_].push_back(rowIt);
      }
    }
    
    void BDAIntervalBuffer::Advance(const double interval)
    {
      time_ += interval_;
      interval_ = interval;
      
      removeOldBaselineRows();
      removeOldBuffers();
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
    
    std::list<BDAIntervalBuffer::BDARowIterator>
    BDAIntervalBuffer::GetBaseline(const std::size_t baseline) const
    {
      if (baseline >= baselines_.size()) {
        throw std::invalid_argument("Invalid baseline");
      }
      
      // Find the first row that has a time greater than or equal to
      // the end of the current interval.
      const auto compare = []( const BDARowIterator& it, const double time)
      {
        return BDABuffer::TimeIsLess(it->time_, time);
      };

      const auto it = std::lower_bound(baselines_[baseline].begin(),
                                       baselines_[baseline].end(),
                                       time_ + interval_,
                                       compare);

      return std::list<BDARowIterator>(baselines_[baseline].begin(), it);
    }
    
    void BDAIntervalBuffer::removeOldBaselineRows()
    {
      for (auto& baseline: baselines_) {
        for (auto rowIterator = baseline.begin();
             rowIterator != baseline.end();
             rowIterator = baseline.erase(rowIterator)) {
          const double kRowEnd = (*rowIterator)->time_ +
                                 (*rowIterator)->interval_;
          if (!BDABuffer::TimeIsLess(kRowEnd, time_)) {
            break;
          }
        }
      }
    }
    
    void BDAIntervalBuffer::removeOldBuffers()
    {
      // If all rows in a buffer have an end time before 'time_',
      // delete the buffer.
      
      auto buffer_iterator = buffers_.begin();
      while(buffer_iterator != buffers_.end()) {
        bool buffer_is_old = true;
        for (const auto& row : (*buffer_iterator)->GetRows()) {
          const double kRowEnd = row.time_ + row.interval_;
          if (BDABuffer::TimeIsLess(time_, kRowEnd)) {
            buffer_is_old = false;
            break;
          }
        }
        if (buffer_is_old) {
          buffer_iterator = buffers_.erase(buffer_iterator);
        } else {
          ++buffer_iterator;
        }
      }
    }
  }
}
