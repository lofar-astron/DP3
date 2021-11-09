// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_BUFFER_H
#define DDECAL_SOLVER_BUFFER_H

#include <complex>
#include <vector>

namespace dp3 {
namespace base {
class DPBuffer;
}

namespace ddecal {

class SolverBuffer {
 public:
  typedef std::complex<float> Complex;

  SolverBuffer();

  /**
   * This function takes (unweighted) data and model data, as well as a
   * weights array, weights these and initializes the buffer from
   * them.
   *
   * @param data_buffers A vector with one buffer for each timestep. Each buffer
   * should contain unweighted data and the corresponding weight values.
   * @param model_buffers Buffers with model data such that
   * model_buffers[time][direction] holds the unweighted model data for the
   * given time step and direction.
   * These buffers have in the same structure as the data.
   * Because the model data is large (e.g. tens of GB in extensive slow gain
   * solves), the data is not copied: The SolverBuffer takes ownership of the
   * model buffers and weights the data in place.
   */
  void AssignAndWeight(
      const std::vector<base::DPBuffer>& data_buffers,
      std::vector<std::vector<base::DPBuffer>>&& model_buffers);

  /**
   * @return The number of time steps in the buffer.
   */
  size_t NTimes() const { return data_.size(); };

  /**
   * @return The number of baselines in the buffer.
   */
  size_t NBaselines() const { return n_baselines_; }

  /**
   * @return The number of channels in the buffer.
   */
  size_t NChannels() const { return n_channels_; }

  /**
   * Get a pointer to the weighted data for a channel.
   * @param time_index The time step index.
   * @param baseline The baseline index in the data.
   * @param channel The channel index in the data.
   * @return A pointer to the correlation value(s) in the requested channel.
   */
  const Complex* DataPointer(size_t time_index, size_t baseline,
                             size_t channel) const;

  /**
   * Copy data from successive channels for all baselines.
   * Use this function instead of using DataPointer() for copying multiple
   * channels, since the internal data format of the SolverBuffer may change.
   * @param time_index The time step index.
   * @param channel_begin Index of the first channel.
   * @param channel_end Index of the next channel. Data up to but not including
   *                    this channel is copied.
   */
  void CopyDataChannels(size_t time_index, size_t channel_begin,
                        size_t channel_end,
                        std::complex<float>* destination) const;

  /**
   * Get a pointer to the weighted model data for a channel.
   * @param time_index The time step index.
   * @param direction The direction index for the model data.
   * @param baseline The baseline index in the model data.
   * @param channel The channel index in the model data.
   * @return A pointer to the correlation value(s) in the requested channel.
   */
  const Complex* ModelDataPointer(size_t time_index, size_t direction,
                                  size_t baseline, size_t channel) const;

 private:
  size_t n_baselines_;
  size_t n_channels_;
  std::vector<std::vector<Complex>> data_;
  std::vector<std::vector<base::DPBuffer>> model_buffers_;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_SOLVER_BUFFER_H
