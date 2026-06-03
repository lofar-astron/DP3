// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_CADENCEIMAGER_H_
#define DP3_STEPS_CADENCEIMAGER_H_

#include "steps/Step.h"

#include "base/DPBuffer.h"
#include "common/ParameterSet.h"
#include "common/Timer.h"

#include <memory>
#include <string>

#include <wsclean/inmemoryms.h>

namespace dp3::steps {

/// @brief DP3 step class that creates images using WSClean.
/// This class provides a step that uses the in-memory interface of WSClean
/// to create images at a preset cadence.
///
/// Most WSClean parameters are currently available, though the in-memory
/// interface does not support all command-line options (e.g. beam evaluations).
///
/// The step creates an in-memory MS that contains relevant metadata, along with
/// visibilities and weights. For each timestep the visibilities in the buffer
/// are copied to the in-memory MS until the desired cadence is obtained. After
/// which it calls WSClean with the provided CLI options.
///
/// Note that BDA is currently not supported.
class ImageStep final : public Step {
 public:
  ImageStep(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kUvwField | kWeightsField | kFlagsField;
  }

  common::Fields getProvidedFields() const final { return {}; }

  bool process(std::unique_ptr<base::DPBuffer>) override;

  // TODO: implement the BDA processing.
  bool process(std::unique_ptr<base::BdaBuffer>) override;

  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

 private:
  /// Initialize an single in-memory MeasurementSet and set the required meta
  /// data.
  void InitializeInMemoryMs();

  /// Add rows pertaining to a single time step to the in-memory MS.
  void AddCurrentBufferToInMemoryMs(const base::DPBuffer::DataType& data,
                                    const base::DPBuffer::WeightsType& weights,
                                    const base::DPBuffer::UvwType& uvws,
                                    double time);

  /// Append a counter to the image name prefix such that each of
  /// the FITS files has a unique prefix.
  std::string FormatImageNamePrefix() const;

  std::string name_;
  common::NSTimer timer_;

  /// The number of seconds after which an image will be created.
  double cadence_;

  std::string wsclean_options_;

  /// FITS filename prefix provided to WSClean.
  /// This will overwrite any prefix specified
  /// in the options.
  std::string fits_prefix_;
  size_t image_counter_ = 0;

  // Number of timesteps added to the in-memory MS.
  // Resets when the desired cadence is reached.
  size_t timestep_counter_ = 0;
  size_t n_timesteps_per_image_;

  bool has_frequency_bda_ = false;
  bool has_time_bda_ = false;

  /// The in-memory MS cache that's being filled with
  /// new rows.
  std::unique_ptr<wsclean::InMemoryMs> in_memory_ms_;
};

}  // namespace dp3::steps

#endif
