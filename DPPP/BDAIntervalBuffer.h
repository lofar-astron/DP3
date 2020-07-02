// BDAIntervalBuffer.h: Provide BDA data for time intervals.
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

/// @file
/// @brief Provide BDA data for time intervals.
/// @author Maik Nijhuis

#ifndef DPPP_BDAINTERVNALBUFFER_H
#define DPPP_BDAINTERVALBUFFER_H

#include "BDABuffer.h"

#include <memory>
#include <list>
#include <set>
#include <vector>

namespace DP3 {
  namespace DPPP {

    class BDAIntervalBuffer
    {
    public:
      using BDARowIterator = std::vector<BDABuffer::Row>::const_iterator;
    
      /**
       * Constructor.
       * @param n_baselines Number of baselines.
       * @param time Start time for the first interval, in seconds.
       * @param interval Duration for the first interval, in seconds.
       */
      BDAIntervalBuffer(std::size_t n_baselines, double time, double interval);

      /**
       * Add new data to the interval.
       * @param buffer A BDA buffer. The interval buffer copies all necessary content.
       * @throw std::invalid_argument If the buffer is invalid.
       */
      void AddBuffer(const BDABuffer& buffer);

      /**
       * Advance the time interval.
       * 
       * The start time of the new interval becomes the end time of the old
       * interval.
       * 
       * This function deletes data that is no longer needed for the new interval
       * and future time intervals.
       * 
       * @param interval Duration for the new interval, in seconds.
       */
      void Advance(double interval);

      /**
       * Check if all data is present for the current time interval.
       * 
       * @return True if all baselines have data that cover the current
       *         time interval, false otherwise.
       */
      bool IsComplete() const;

      /**
       * Get the BDA rows that match the current time interval for one baseline.
       * @param baseline The index of the requested baseline.
       * @throw std::invalid_argument If baseline is invalid.
       */
      std::list<BDARowIterator> GetBaseline(std::size_t baseline) const;

    private:
      void removeOldBaselineRows();
      void removeOldBuffers();

    private:
      double time_; ///< Start time of current interval.
      double interval_; ///< Duration of current interval.

      std::list<std::unique_ptr<BDABuffer>> buffers_;

      /// Contains an ordered collection of all rows for each baseline.
      std::vector<std::list<BDARowIterator>> baselines_;
    };

  }
}

#endif