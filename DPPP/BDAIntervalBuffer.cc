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

namespace{
  constexpr double epsilon = 1.0e-8; // For comparing measurement timestamps.
  
  constexpr bool IsLessThan(const double x, const double y) {
    return x < (y - epsilon);
  }
  
  constexpr bool IsGreaterEqual(const double x, const double y) {
    return x > (y - epsilon);
  }
  
  constexpr bool IsEqual(const double x, const double y) {
    return abs(x - y) < epsilon;
  }
}

namespace DP3 {
  namespace DPPP {

    BDAIntervalBuffer::BDAIntervalBuffer(const std::size_t nBaselines,
                                         const double time,
                                         const double interval)
    : time_(time)
    , interval_(interval)
    , buffers_()
    , baselines_(nBaselines)
    , incomplete_()
    {
      initializeIncomplete();
    }
    
    void BDAIntervalBuffer::AddBuffer(const BDABuffer& buffer)
    {
      // Check if 'buffer' is valid.
      for(const auto& row : buffer.GetRows()) {
        if (row.baseline_nr_ >= baselines_.size()) {
          throw std::invalid_argument("Invalid baseline");
        }

        const auto& baseline = baselines_[row.baseline_nr_];
        if (!baseline.empty() &&
            !IsEqual(baseline.back()->time_ + baseline.back()->interval_,
                     row.time_)) {
          throw std::invalid_argument("Start time does not match end time of previous row.");
        }
      }
              
      buffers_.push_back(boost::make_unique<BDABuffer>(buffer));
      
      // Add references in baselines_ to the rows of the new buffer.
      for(auto rowIt = buffers_.back()->GetRows().begin();
          rowIt != buffers_.back()->GetRows().end();
          ++rowIt) {
        baselines_[rowIt->baseline_nr_].push_back(rowIt);
      }
      
      // Update incomplete_
      auto incompleteIterator = incomplete_.begin();
      while(incompleteIterator != incomplete_.end()) { 
        if (baselineIsComplete(*incompleteIterator)) {
          incompleteIterator = incomplete_.erase(incompleteIterator);
        } else {
          ++incompleteIterator;
        }
      }
    }
    
    void BDAIntervalBuffer::Advance(const double interval)
    {
      time_ += interval_;
      interval_ = interval;
      
      removeOldBaselineRows();
      removeOldBuffers();

      initializeIncomplete();
    }
    
    std::list<BDAIntervalBuffer::BDARowIterator>
    BDAIntervalBuffer::GetBaseline(const std::size_t baseline)
    {
      if (baseline >= baselines_.size()) {
        throw std::invalid_argument("Invalid baseline");
      }
      
      // Find the first row that has a time greater than or equal to
      // the end of the current interval.
      const auto compare = []( const BDARowIterator& it, const double time)
      {
        return IsLessThan(it->time_, time);
      };

      const auto it = std::lower_bound(baselines_[baseline].begin(),
                                       baselines_[baseline].end(),
                                       time_ + interval_,
                                       compare);

      return std::list<BDARowIterator>(baselines_[baseline].begin(), it);
    }
    
    bool BDAIntervalBuffer::baselineIsComplete(const std::size_t baselineNr) const
    {
      // AddBuffer() already checks that all rows for a baseline are
      // contiguous: Here, we only have to check the last row.
      assert(baselineNr < baselines_.size());
      const auto& rows = baselines_[baselineNr];
      return !rows.empty() &&
             IsGreaterEqual(rows.back()->time_ + rows.back()->interval_,
                            time_ + interval_);
    }
    
    void BDAIntervalBuffer::initializeIncomplete()
    {
      for (std::size_t i = 0; i < baselines_.size(); ++i) {
        if (baselineIsComplete(i)) {
          incomplete_.erase(i);
        } else {
          incomplete_.insert(i);
        }
      }
    }
    
    void BDAIntervalBuffer::removeOldBaselineRows()
    {
      for (auto& baseline: baselines_) {
        for (auto rowIterator = baseline.begin();
             rowIterator != baseline.end() &&
             IsLessThan((*rowIterator)->time_ + (*rowIterator)->interval_,
                        time_);
             rowIterator = baseline.erase(rowIterator));
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
          if (IsLessThan(time_, row.time_ + row.interval_)) {
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
