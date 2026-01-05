// NullStokes.h: DP3 step class for zeroing out Stokes parameters
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class for zeroing out Sokes parameters
/// @author Matthijs van der Wild

#ifndef DP3_STEPS_NULL_STOKES_H
#define DP3_STEPS_NULL_STOKES_H

#include "InputStep.h"
#include "common/ParameterSet.h"

#include "base/DPBuffer.h"

namespace dp3 {

namespace steps {
/// @brief DP3 step for polarisation modification

/// This class is a Step class that optionally sets
/// Stokes Q or U to zero.

class NullStokes : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  explicit NullStokes(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override { return kDataField; }

  common::Fields getProvidedFields() const override { return kDataField; }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps
  void finish() override;

  /// Update the general info
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters
  void show(std::ostream&) const override;

  /// Show the timings
  void showTimings(std::ostream&, double duration) const override;

 private:
  std::string name_;
  common::NSTimer timer_;
  bool modify_q_;
  bool modify_u_;
};

}  // namespace steps

}  // namespace dp3

#endif
