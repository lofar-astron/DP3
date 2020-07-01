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

#include <algorithm>
#include <cassert>

namespace{
  constexpr double epsilon = 1.0e-8; // For comparing measurement timestamps.
  
  constexpr bool isLessThan(const double x, const double y) {
    return x < (y - epsilon);
  }
  
  constexpr bool isGreaterEqual(const double x, const double y) {
    return x > (y - epsilon);
  }
  
  constexpr bool isEqual(const double x, const double y) {
    return abs(x - y) < epsilon;
  }
}

namespace DP3 {
  namespace DPPP {

    BDAIntervalBuffer::BDAIntervalBuffer(const std::size_t nBaselines,
                                         const double time,
                                         const double exposure)
    : itsTime(time)
    , itsExposure(exposure)
    , itsBuffers()
    , itsBaselines(nBaselines)
    , itsIncomplete()
    {
      initializeIncomplete();
    }
    
    void BDAIntervalBuffer::addBuffer(std::unique_ptr<BDABuffer> buffer)
    {
      // Check if 'buffer' is valid.
      for(const auto& row : buffer->getRows()) {
        if (row.itsBaselineNr >= itsBaselines.size()) {
          throw std::invalid_argument("Invalid baseline");
        }

        const auto& baseline = itsBaselines[row.itsBaselineNr];
        if (!baseline.empty() &&
            !isEqual(baseline.back()->itsTime + baseline.back()->itsExposure,
                     row.itsTime)) {
          throw std::invalid_argument("Start time does not match end time of previous row.");
        }
      }
              
      // Add references to the rows in itsBaselines.
      for(auto rowIt = buffer->getRows().begin();
          rowIt != buffer->getRows().end();
          ++rowIt) {
        itsBaselines[rowIt->itsBaselineNr].push_back(rowIt);
      }
      
      itsBuffers.push_back(std::move(buffer));
      
      // Update itsIncomplete
      auto incompleteIterator = itsIncomplete.begin();
      while(incompleteIterator != itsIncomplete.end()) { 
        if (baselineIsComplete(*incompleteIterator)) {
          incompleteIterator = itsIncomplete.erase(incompleteIterator);
        } else {
          ++incompleteIterator;
        }
      }
    }
    
    void BDAIntervalBuffer::advance(const double exposure)
    {
      itsTime += itsExposure;
      itsExposure = exposure;
      
      removeOldBaselineRows();
      removeOldBuffers();

      initializeIncomplete();
    }
    
    std::list<BDAIntervalBuffer::BDARowIterator>
    BDAIntervalBuffer::getBaseline(const std::size_t baseline)
    {
      if (baseline >= itsBaselines.size()) {
        throw std::invalid_argument("Invalid baseline");
      }
      
      // Find the first row that has a time greater than or equal to
      // the end of the current interval.
      const auto compare = []( const BDARowIterator& it, const double time)
      {
        return isLessThan(it->itsTime, time);
      };

      const auto it = std::lower_bound(itsBaselines[baseline].begin(),
                                       itsBaselines[baseline].end(),
                                       itsTime + itsExposure,
                                       compare);

      return std::list<BDARowIterator>(itsBaselines[baseline].begin(), it);
    }
    
    bool BDAIntervalBuffer::baselineIsComplete(const std::size_t baselineNr) const
    {
      // addBuffer() already checks that all rows for a baseline are
      // contiguous: Here, we only have to check the last row.
      assert(baselineNr < itsBaselines.size());
      const auto& rows = itsBaselines[baselineNr];
      return !rows.empty() &&
             isGreaterEqual(rows.back()->itsTime + rows.back()->itsExposure,
                            itsTime + itsExposure);
    }
    
    void BDAIntervalBuffer::initializeIncomplete()
    {
      for (std::size_t i = 0; i < itsBaselines.size(); ++i) {
        if (baselineIsComplete(i)) {
          itsIncomplete.erase(i);
        } else {
          itsIncomplete.insert(i);
        }
      }
    }
    
    void BDAIntervalBuffer::removeOldBaselineRows()
    {
      for (auto& baseline: itsBaselines) {
        for (auto rowIterator = baseline.begin();
             rowIterator != baseline.end() &&
             isLessThan((*rowIterator)->itsTime + (*rowIterator)->itsExposure,
                        itsTime);
             rowIterator = baseline.erase(rowIterator));
      }
    }
    
    void BDAIntervalBuffer::removeOldBuffers()
    {
      // If all rows in a buffer have an end time before 'itsTime',
      // delete the buffer.
      
      auto bufferIterator = itsBuffers.begin();
      while(bufferIterator != itsBuffers.end()) {
        bool bufferIsOld = true;
        for (const auto& row : (*bufferIterator)->getRows()) {
          if (isLessThan(itsTime, row.itsTime + row.itsExposure)) {
            bufferIsOld = false;
            break;
          }
        }
        if (bufferIsOld) {
          bufferIterator = itsBuffers.erase(bufferIterator);
        } else {
          ++bufferIterator;
        }
      }
    }
  }
}
