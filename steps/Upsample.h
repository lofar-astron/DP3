// Upsample.h: DPPP step class to upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to Upsample visibilities
/// @author Tammo Jan Dijkema

#ifndef DPPP_DummyStep_H
#define DPPP_DummyStep_H

#include "InputStep.h"

#include <dp3/base/DPBuffer.h>

#include <utility>

namespace dp3 {
namespace base {
class UVWCalculator;
}
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DPPP step class to Upsample visibilities
class Upsample : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Upsample(const common::ParameterSet&, const std::string& prefix);

  /// Constructor that gets the settings directly.
  Upsample(const std::string& name, unsigned int time_step, bool update_uvw);

  ~Upsample() override;

  common::Fields getRequiredFields() const override {
    common::Fields fields = kFlagsField;
    if (update_uvw_) fields |= kUvwField;
    return fields;
  }

  common::Fields getProvidedFields() const override {
    return update_uvw_ ? kUvwField : common::Fields();
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

 private:
  /// Update the time and exposure of a buffer and if `update_uvw_` is set also
  /// the uvw. Used internally by process function to populate each of
  /// `buffers_` for the appropriate time step.
  void UpdateTimeCentroidExposureAndUvw(std::unique_ptr<base::DPBuffer>& buffer,
                                        double time, double exposure);

  const std::string name_;
  const unsigned int time_step_;
  const bool update_uvw_;

  std::vector<std::unique_ptr<base::DPBuffer>> prev_buffers_;
  std::vector<std::unique_ptr<base::DPBuffer>> buffers_;
  unsigned int first_to_flush_;
  std::unique_ptr<base::UVWCalculator> uvw_calculator_;

  common::NSTimer timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
