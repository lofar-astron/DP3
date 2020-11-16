// BDAIntervalBuffer.h: Provide BDA data for time intervals.
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Provide BDA data for time intervals.
/// @author Maik Nijhuis

#ifndef DPPP_BDAINTERVALBUFFER_H
#define DPPP_BDAINTERVALBUFFER_H

#include "BDABuffer.h"

#include <boost/noncopyable.hpp>

#include <memory>
#include <list>
#include <set>
#include <vector>

namespace DP3 {
namespace DPPP {

class BDAIntervalBuffer : private boost::noncopyable {
 public:
  /**
   * Constructor.
   * @param time Start time for the first interval, in seconds.
   * @param interval Duration for the first interval, in seconds.
   * @param max_row_interval The maximum duration for a BDABuffer row.
   */
  BDAIntervalBuffer(double time, double interval, double max_row_interval);

  /**
   * Add new data to the interval.
   * @param buffer A BDA buffer. The interval buffer copies all necessary
   * content.
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
   * Get a buffer containing the weighted data for the current interval.
   * @param fields Bitset with the requested fields.
   * @return A BDABuffer with the requested data.
   */
  std::unique_ptr<BDABuffer> GetBuffer(
      const BDABuffer::Fields& fields = BDABuffer::Fields()) const;

  /**
   * Access the number of internally stored buffers.
   *
   * Tests use this function for checking that Advance() releases old data.
   * @return The number of stored buffers.
   */
  std::size_t GetPendingBufferCount() const { return buffers_.size(); }

 private:
  enum class Completeness { kUnknown, kComplete, kIncomplete };
  mutable Completeness completeness_;  ///< Cached completeness status.

  double time_;                    ///< Start time of current interval.
  double interval_;                ///< Duration of current interval.
  const double max_row_interval_;  ///< Maximum duration of a BDA row.

  std::list<BDABuffer> buffers_;  ////< Main storage of this class.

  /// Contains all rows for the current interval.
  std::vector<const BDABuffer::Row*> current_rows_;
};

}  // namespace DPPP
}  // namespace DP3

#endif