// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_SOLVER_BUFFER_H
#define DDECAL_BDA_SOLVER_BUFFER_H

#include "../../base/BDABuffer.h"

#include <aocommon/queue.h>

#include <complex>
#include <memory>
#include <vector>

namespace dp3 {
namespace base {
namespace test {
class SolverTester;
}

class BDASolverBuffer {
 public:
  /**
   * Constructor.
   * @param n_directions The number of directions.
   * @param start Start time of the first solution interval.
   * @param interval Length of a solution interval. Should be > 0.0.
   */
  BDASolverBuffer(size_t n_directions, double start, double interval)
      : data_(),
        model_data_(n_directions),
        time_start_(start),
        time_interval_(interval),
        last_complete_interval_(-1),
        data_rows_(),
        model_rows_(n_directions) {
    assert(interval >= 0.0);
    AddInterval();  // Ensure that the current solution interval is valid.
  }

  /**
   * This function takes a buffer with unweighted data and the corresponding
   * weights, and a buffer with model data. It weights these data buffers
   * and stores the result internally.
   * @param data_buffer A buffer with unweighted data and weights.
   * @param model_buffers A vector with model_buffers for each direction.
   * The BDA layout of these buffers should match the layout of the data_buffer.
   * The BDASolverBuffer takes ownership of the model buffers.
   * @throw std::invalid_argument If model_buffers has an invalid size.
   */
  void AppendAndWeight(const BDABuffer& data_buffer,
                       std::vector<std::unique_ptr<BDABuffer>>&& model_buffers);

  /**
   * Clears all internal buffers.
   * Does not affect the solution interval and the number of directions.
   */
  void Clear();

  /**
   * @return True if the current solution interval is complete. An interval
   * is complete if the solver buffer contains a BDA row with a start time
   * greater or equal to the end time of the current solution interval.
   */
  bool IntervalIsComplete() const { return last_complete_interval_ >= 0; }

  /**
   * Advances the current solution interval to the next interval.
   * Releases all internal buffers that only hold data for previous solution
   * intervals.
   */
  void AdvanceInterval();

  /**
   * Get the number of active buffers, which is the number of AppendAndWeight()
   * calls minus the number of buffers that AdvanceInterval() released.
   * This function mainly exists for testing AdvanceInterval().
   * @return The number of active buffers.
   */
  size_t BufferCount() const { return data_.Size(); }

  /**
   * Get the data for the current solution interval.
   * @return Non-modifyable rows with weighted visibilities.
   */
  const std::vector<const BDABuffer::Row*>& GetDataRows() const {
    return data_rows_[0];
  }

  /**
   * Get the model data rows for the current solution interval.
   * @param direction Direction index.
   * @return Non-modifyable rows with weighted model data.
   */
  const std::vector<const BDABuffer::Row*>& GetModelDataRows(
      size_t direction) const {
    return model_rows_[direction][0];
  }

 private:
  void AddInterval();

  /**
   * data_ is a FIFO queue with weighted input data.
   * AppendAndWeight appends items at the end.
   */
  aocommon::Queue<std::unique_ptr<BDABuffer>> data_;
  /**
   * For each direction, a FIFO queue with model data.
   * These queues are managed the same as data_.
   */
  std::vector<aocommon::Queue<std::unique_ptr<BDABuffer>>> model_data_;

  double time_start_;  ///< Start time of current solution interval (seconds).
  double time_interval_;  ///< Solution interval length (seconds).
  /// Index of the last complete solution interval. The value is negative if no
  /// interval is complete.
  int last_complete_interval_;

  /// The data rows for the current and future solution intervals.
  /// The queue always contains one element for the current solution interval.
  aocommon::Queue<std::vector<const BDABuffer::Row*>> data_rows_;

  /// For each direction, the model data rows for each solution interval.
  /// The queues always contain one element for the current solution interval.
  std::vector<aocommon::Queue<std::vector<const BDABuffer::Row*>>> model_rows_;
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_BDA_SOLVER_BUFFER_H
