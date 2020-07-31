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
/// @brief Step for compressing regular data into BDA data.
/// @author Maik Nijhuis

#ifndef DPPP_BDAAVERAGER_H
#define DPPP_BDAAVERAGER_H

#include "DPStep.h"

#include <vector>

namespace DP3 {
namespace Common {
class ParameterSet;
}  // namespace Common
}  // namespace DP3

namespace DP3 {
namespace DPPP {

class BDABuffer;

class BDAAverager : public DPStep {
 public:
  /**
   * Constructor.
   * @param bl_threshold_time The baseline length after which the averager
   * does not perform any time averaging. The time averaging factor for a
   * baseline is bl_threshold_time / baseline length, rounded down.
   * @param bl_threshold_channels The baseline length after which the
   * averager does not perform any channel averaging. The number of channels for
   * a baseline becomes baseline length / bl_threshold_channels * original
   * channel count, rounded up.
   */
  BDAAverager(double bl_threshold_time, double bl_threshold_channels);

  ~BDAAverager() override;

  bool process(const DPBuffer&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void updateInfo(const DPInfo&) override;

 private:
  struct BaselineBuffer {
    BaselineBuffer(std::size_t time_factor, std::size_t n_input_channels,
                   std::size_t n_output_channels, std::size_t n_correlations);
    void Clear();

    std::size_t times_added;        ///< Number of added regular intervals.
    const std::size_t time_factor;  ///< Time averaging factor.
    /// Input channel start and end index for each output channel.
    std::vector<std::size_t> input_channel_indices;
    double time;
    double interval;
    std::vector<std::complex<float>> data;
    std::vector<float> weights;
    double uvw[3];
  };

  void AddBaseline(std::size_t baseline_nr);

  /// Baseline threshold length for time averaging.
  const double bl_threshold_time_;
  /// Baseline threshold length for channel averaging.
  const double bl_threshold_channel_;
  rownr_t next_rownr_;
  std::size_t bda_pool_size_;
  std::unique_ptr<BDABuffer> bda_buffer_;
  std::vector<BaselineBuffer> baseline_buffers_;
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DPPP_BDAAVERAGER_H
