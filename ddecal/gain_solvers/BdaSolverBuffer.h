// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_SOLVER_BUFFER_H
#define DDECAL_BDA_SOLVER_BUFFER_H

#include "../../base/BDABuffer.h"

#include <aocommon/queue.h>

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

namespace dp3 {
namespace ddecal {

class BdaSolverBuffer {
 public:
  /**
   * Constructor.
   * @param n_directions The number of directions.
   * @param start Start time of the first solution interval.
   * @param interval Length of a solution interval. Should be > 0.0.
   */
  BdaSolverBuffer(size_t n_directions, double start, double interval)
      : data_(),
        time_start_(start),
        time_interval_(interval),
        current_interval_(0),
        last_complete_interval_(-1),
        data_rows_() {
    assert(interval >= 0.0);
    AddInterval(
        n_directions);  // Ensure that the current solution interval is valid.
  }

  /**
   * This function takes a buffer with unweighted data and the corresponding
   * weights, and a buffer with model data. It weights these data buffers
   * and stores the result internally.
   * @param data_buffer A buffer with unweighted data and weights.
   * @param model_buffers A vector with model_buffers for each direction.
   * The BDA layout of these buffers should match the layout of the data_buffer.
   * The BdaSolverBuffer takes ownership of the model buffers.
   * @throw std::invalid_argument If model_buffers has an invalid size.
   */
  void AppendAndWeight(
      std::unique_ptr<base::BDABuffer> data_buffer,
      std::vector<std::unique_ptr<base::BDABuffer>>&& model_buffers);

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
  const std::vector<const base::BDABuffer::Row*>& GetDataRows() const {
    return data_rows_[0].weighted;
  }

  /**
   * Get the model data rows for the current solution interval.
   * @param direction Direction index.
   * @return Rows with weighted model data.
   */
  const std::vector<const base::BDABuffer::Row*>& GetModelDataRows(
      size_t direction) const {
    return data_rows_[0].model[direction];
  }

  void SubtractCorrectedModel(
      const std::vector<std::vector<std::complex<float>>>& solutions,
      size_t n_channel_blocks, bool full_jones, const std::vector<int>& ant1,
      const std::vector<int>& ant2);

  std::vector<std::unique_ptr<base::BDABuffer>> GetDone() {
    std::vector<std::unique_ptr<base::BDABuffer>> result;
    result.swap(done_);
    return result;
  }

 private:
  void AddInterval(size_t n_directions);

  /**
   * @return The relative solution interval index for a given time.
   */
  int RelativeIndex(double time) const {
    return int(std::floor((time - time_start_) / time_interval_)) -
           current_interval_;
  }

  /**
   * The BDASolverBuffer is the owner of both the unweighted and the weighted
   * input data. The data flow is:
   * 1. AppendAndWeight receives unweighted data, weights it, and stores both
   *    the unweighted and the corresponding weighted data in data_.
   * 2. A Solver accesses the weighted data via GetDataRows and
   *    GetModelDataRows, via a SolveData structure.
   * 3. SubtractCorrectedModel() applies the solutions to the model data and
   *    subtracts them from the unweighted data.
   * 4. When a BDA buffer is no longer needed, AdvanceInterval deletes the
   *    weighted data and moves the unweighted data to done_.
   * 5. GetDone() extracts the buffers from done_.
   */
  struct InputData {
    std::unique_ptr<base::BDABuffer> unweighted;
    std::unique_ptr<base::BDABuffer> weighted;
    /// Model data buffer for each direction.
    std::vector<std::unique_ptr<base::BDABuffer>> model;
  };

  /// A FIFO queue with input data. AppendAndWeight appends items.
  aocommon::Queue<InputData> data_;

  /// Fully processed input buffers.
  std::vector<std::unique_ptr<base::BDABuffer>> done_;

  /// Start time of the first solution interval (seconds).
  const double time_start_;
  const double time_interval_;  ///< Solution interval length (seconds).
  int current_interval_;  ///< Absolute index of the current solution interval.

  /// Index of the last complete solution interval, relative to the current
  /// solution interval. The value is negative if no interval is complete.
  int last_complete_interval_;

  struct IntervalRows {
    std::vector<base::BDABuffer::Row*> unweighted;
    std::vector<const base::BDABuffer::Row*> weighted;
    std::vector<std::vector<const base::BDABuffer::Row*>> model;
  };

  /// The data rows for the current and future solution intervals.
  /// The queue always contains one element for the current solution interval.
  aocommon::Queue<IntervalRows> data_rows_;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_BDA_SOLVER_BUFFER_H
