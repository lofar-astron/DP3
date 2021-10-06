// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for compressing regular data into BDA data.
/// @author Maik Nijhuis

#ifndef DPPP_BDAAVERAGER_H
#define DPPP_BDAAVERAGER_H

#include "Step.h"

#include <casacore/casa/Arrays/IPosition.h>

#include <vector>
#include <queue>

namespace dp3 {
namespace common {
class ParameterSet;
}
namespace base {
class BDABuffer;
}
}  // namespace dp3

namespace dp3 {
namespace steps {

class InputStep;

class BDAAverager : public Step {
 public:
  /**
   * Constructor, which uses a parset for configuring the step.
   * @param input InputStep object, for fetching weights, UVW etc.
   * @param parset A ParameterSet that contains the configuration.
   * @param prefix ParameterSet Prefix for obtaining the configuration.
   * @param use_weights_and_flags A flag (true by default) which allows the
   * BdaAverager to ignore the weights and flags. When false, it assumes
   * unflagged data and a weight of 1.0 for all input data.
   */
  BDAAverager(InputStep& input, const common::ParameterSet& parset,
              const std::string& prefix,
              const bool use_weights_and_flags = true);

  ~BDAAverager() override;

  bool process(const base::DPBuffer&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType t) const override { return t == MsType::kRegular; }

  MsType outputs() const override { return MsType::kBda; };

  void set_averaging_params(std::vector<unsigned int> baseline_factors,
                            std::vector<std::vector<double>> freqs,
                            std::vector<std::vector<double>> widths);

  /**
   * Public method, which sets a desired output size (number of rows).
   * @param buffersize Number of rows in the output buffer: these should be
   * given in the order one wishes to see in the output. If the BDAAverager is
   * ready to output a bdabuffer but no size is available, the default value
   * will be used.
   */
  void set_next_desired_buffersize(unsigned int buffersize);

 private:
  struct BaselineBuffer {
    BaselineBuffer(std::size_t time_factor, std::size_t n_input_channels,
                   std::size_t n_output_channels, std::size_t n_correlations);
    void Clear();

    std::size_t times_added;        ///< Number of added regular intervals.
    const std::size_t time_factor;  ///< Time averaging factor.
    /// Input channel start and end index for each output channel.
    std::vector<std::size_t> input_channel_indices;
    double starttime;
    double interval;
    double exposure;
    std::vector<std::complex<float>> data;
    std::vector<float> weights;
    double uvw[3];
  };

  void AddBaseline(std::size_t baseline_nr);

  InputStep& input_;
  common::NSTimer timer_;

  /// Baseline threshold length for time averaging.
  const double bl_threshold_time_;
  /// Baseline threshold length for channel averaging.
  const double bl_threshold_channel_;
  /// Maximum interval / exposure time when doing time averaging.
  const double max_interval_;
  /// Minimum number of channels in the output.
  std::size_t min_channels_;
  /// Name of the step
  std::string name_;

  /// Maximum frequency factor computed from actual baseline lengths
  std::size_t maxfreqfactor_;
  /// Maximum time factor computed from actual baseline lengths
  std::size_t maxtimefactor_;

  common::rownr_t next_rownr_;
  std::size_t bda_pool_size_;
  std::unique_ptr<base::BDABuffer> bda_buffer_;
  std::vector<BaselineBuffer> baseline_buffers_;

  casacore::IPosition expected_input_shape_;

  /// Time averaging factors per baseline (temporarily used to store the
  /// arguments of set_averaging_params())
  std::vector<unsigned int> baseline_factors_;
  /// Center frequency per channel per baseline (temporarily used to store the
  /// arguments of set_averaging_params())
  std::vector<std::vector<double>> freqs_;
  /// Channel width per channel per baseline (temporarily used to store the
  /// arguments of set_averaging_params())
  std::vector<std::vector<double>> widths_;

  const bool use_weights_and_flags_;

  std::queue<std::unique_ptr<base::BDABuffer>> fixed_size_bdabuffers_;
};

}  // namespace steps
}  // namespace dp3

#endif  // DPPP_BDAAVERAGER_H
