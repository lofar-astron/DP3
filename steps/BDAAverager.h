// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for compressing regular data into BDA data.
/// @author Maik Nijhuis

#ifndef DP3_STEPS_BDAAVERAGER_H_
#define DP3_STEPS_BDAAVERAGER_H_

#include "Step.h"
#include "common/Timer.h"

#include <vector>
#include <queue>

namespace dp3 {
namespace common {
class ParameterSet;
}
}  // namespace dp3

namespace dp3 {
namespace steps {

class BdaAverager : public Step {
 public:
  /**
   * Constructor, which uses a parset for configuring the step.
   * @param parset A ParameterSet that contains the configuration.
   * @param prefix ParameterSet Prefix for obtaining the configuration.
   * @param use_weights_and_flags A flag (true by default) which allows the
   * BdaAverager to ignore the weights and flags. When false, it assumes
   * unflagged data and a weight of 1.0 for all input data.
   */
  BdaAverager(const common::ParameterSet& parset, const std::string& prefix,
              const bool use_weights_and_flags = true);

  ~BdaAverager() override;

  common::Fields getRequiredFields() const override {
    common::Fields fields = kDataField | kUvwField;
    if (use_weights_and_flags_) fields |= kFlagsField | kWeightsField;
    return fields;
  }

  common::Fields getProvidedFields() const override {
    // BdaAverager always creates BdaBuffers with all fields.
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void finish() override;

  void show(std::ostream&) const override;

  void updateInfo(const base::DPInfo&) override;

  bool accepts(MsType t) const override { return t == MsType::kRegular; }

  MsType outputs() const override { return MsType::kBda; };

  /**
   * Set the averaging scheme for the BdaAverager.
   * Using this function, a step can internally expand BDA data, process
   * that data using a regular step, and apply the BdaAverager again.
   * SetAveragingParameters can then restore the original averaging scheme so
   * the resulting data is compatible with the original BDA data.
   * Notes:
   * - This function must be called before calling updateInfo().
   * - BdaAverager only supports averaging schemes it would create normally,
   *   without using SetAveragingParameters. Using an incompatible scheme will
   *   make updateInfo() throw an exception.
   * @param info A DPInfo object which contains an existing averaging scheme.
   */
  void SetAveragingParameters(const base::DPInfo& info);

  /**
   * Pushes a request for an output size (number of elements). The sizes are
   * used one by one until the request buffer is empty. When it is empty, a
   * default size is used. This can be used to make output buffers of the
   * same shape and ordering as another averaging step.
   */
  void PushBufferSizeRequest(size_t buffersize);

  /**
   * Computes the total averaging factor (defined as the ratio between
   * non-averaged and averaged visibility counts) for all baselines.
   * Only works after setting the info using updateInfo().
   * @return The averaging factor. A value of 4.2 means 4.2 megabytes of input
   *         visibilities get averaged into 1 megabyte of output visibilities.
   */
  float TotalAveragingFactor() const;

 private:
  struct BaselineBuffer {
    BaselineBuffer(std::size_t time_factor, std::size_t n_input_channels,
                   std::size_t n_output_channels, std::size_t n_correlations);
    void Clear();

    std::size_t times_added;        ///< Number of added regular intervals.
    const std::size_t time_factor;  ///< Time averaging factor.
    /// Input channel start and end index for each output channel.
    /// For example, [0, 3, 5] means there are 5 input channels and 2
    /// averaged output channels.
    /// The first output channel contains the average of input channels 0, 1, 2.
    /// The second output channel contains the average of input channels 3, 4.
    std::vector<std::size_t> input_channel_indices;
    double starttime;
    double interval;
    double exposure;
    std::map<std::string, std::vector<std::complex<float>>> data;
    std::vector<float> weights;
    double uvw[3];
  };

  void AddBaseline(std::size_t baseline_nr);

  void InitializeBdaBuffer(const std::vector<std::string>& data_names);

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
  std::unique_ptr<base::BdaBuffer> bda_buffer_;
  std::vector<BaselineBuffer> baseline_buffers_;

  std::array<std::size_t, 3> expected_input_shape_;

  /// Time averaging factors per baseline (temporarily used to store the
  /// arguments of SetAveragingParameters())
  std::vector<unsigned int> baseline_factors_;
  /// Center frequency per channel per baseline (temporarily used to store the
  /// arguments of SetAveragingParameters())
  std::vector<std::vector<double>> freqs_;
  /// Channel width per channel per baseline (temporarily used to store the
  /// arguments of SetAveragingParameters())
  std::vector<std::vector<double>> widths_;

  const bool use_weights_and_flags_;

  std::queue<size_t> size_requests_;
  std::vector<std::string> data_names_;
};

}  // namespace steps
}  // namespace dp3

#endif  // DP3_STEPS_BDAAVERAGER_H_
