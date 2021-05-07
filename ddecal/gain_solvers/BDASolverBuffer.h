// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_BDA_SOLVER_BUFFER_H
#define DDECAL_BDA_SOLVER_BUFFER_H

#include "../../base/BDABuffer.h"

#include <complex>
#include <memory>
#include <vector>

namespace dp3 {
namespace base {

class BDASolverBuffer {
 public:
  BDASolverBuffer()
      : data_(), model_data_(), time_start_(0.0), time_interval_(0.0) {}

  /**
   * Set the number of directions.
   * @pre Only call this function when the solver buffer is empty.
   * @param directions The number of directions.
   * @throw std::logic_error If the solver buffer is not empty.
   */
  void SetDirections(size_t directions) {
    if (!data_.empty())
      throw std::logic_error(
          "Buffer is not empty while setting number of directions.");

    model_data_.resize(directions);
    current_model_rows_.resize(directions);
    future_model_rows_.resize(directions);
  }

  /**
   * Set the (first) time interval for the solver buffer.
   * @pre Only call this function when the solver buffer is empty.
   * @param start Start time of the first solution interval.
   * @param interval Length of a solution interval.
   * @throw std::invalid_argument If the solution interval is <= 0.0.
   * @throw std::logic_error If the solver buffer is not empty.
   */
  void SetInterval(double start, double interval) {
    if (interval <= 0.0)
      throw std::invalid_argument("Invalid solution interval: " +
                                  std::to_string(interval));
    if (!data_.empty())
      throw std::logic_error(
          "Buffer is not empty while setting solution interval.");

    time_start_ = start;
    time_interval_ = interval;
  }

  /**
   * This function takes a buffer with unweighted data and the corresponding
   * weights, and a buffer with model data. It weights these data buffers
   * and stores the result internally.
   * @pre Only call this function after setting the solution interval.
   * @param data_buffer A buffer with unweighted data and weights.
   * @param model_buffers A vector with model_buffers for each direction.
   * The BDA layout of these buffers should match the layout of the data_buffer.
   * The BDASolverBuffer takes ownership of the model buffers.
   * @throw std::logic_error If the solution interval is not set.
   * @throw std::invalid_argument If model_buffers has an invalid size.
   */
  void AppendAndWeight(const BDABuffer& data_buffer,
                       std::vector<std::unique_ptr<BDABuffer>>&& model_buffers);

  /**
   * Clears all internal buffers.
   */
  void Clear() {
    data_.clear();
    current_data_rows_.clear();
    future_data_rows_.clear();
    for (auto& buffers : model_data_) buffers.clear();
    for (auto& rows : current_model_rows_) rows.clear();
    for (auto& rows : future_model_rows_) rows.clear();
  }

  /**
   * @return The weighted data buffers.
   */
  const std::vector<std::unique_ptr<BDABuffer>>& GetData() const {
    return data_;
  }

  /**
   * Get the data for the current solution interval.
   * @return Non-modifyable rows with weighted visibilities.
   */
  const std::vector<const BDABuffer::Row*>& GetDataRows() const {
    return current_data_rows_;
  }

  /**
   * @param direction Direction index.
   * @return The weighted model buffers for the given direction.
   */
  const std::vector<std::unique_ptr<BDABuffer>>& GetModelData(
      size_t direction) {
    return model_data_[direction];
  }

  /**
   * Get the model data rows for the current solution interval.
   * @param direction Direction index.
   * @return Modifyable rows with weighted model data.
   */
  const std::vector<const BDABuffer::Row*>& GetModelDataRows(size_t direction) {
    return current_model_rows_[direction];
  }

 private:
  /**
   * data_ is a FIFO queue with weighted input data.
   * AppendAndWeight appends items at the end.
   * Advance() may remove items at the front and then update 'first'.
   * If half the vector is unused, it will move all items to the front, resize,
   * the vector and set 'first' to 0.
   */
  std::vector<std::unique_ptr<BDABuffer>> data_;
  /**
   * For each direction, a FIFO queue with model data.
   * These queues are managed the same as data_.
   */
  std::vector<std::vector<std::unique_ptr<BDABuffer>>> model_data_;
  size_t
      first;  ///< Index of the first item in data_ and model_data_[direction].

  double time_start_;  ///< Start time of current solution interval (seconds).
  double time_interval_;  ///< Solution interval length (seconds).

  /// The data rows for the current solution interval.
  std::vector<const BDABuffer::Row*> current_data_rows_;
  /// The data rows for future solution intervals.
  std::vector<const BDABuffer::Row*> future_data_rows_;

  /// For each direction, the model data rows for the current solution interval.
  std::vector<std::vector<const BDABuffer::Row*>> current_model_rows_;
  /// For each direction, the model data rows for future solution intervals.
  std::vector<std::vector<const BDABuffer::Row*>> future_model_rows_;
};

}  // namespace base
}  // namespace dp3

#endif  // DDECAL_BDA_SOLVER_BUFFER_H
