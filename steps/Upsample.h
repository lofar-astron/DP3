// Upsample.h: DPPP step class to upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to Upsample visibilities
/// @author Tammo Jan Dijkema

#ifndef DPPP_DummyStep_H
#define DPPP_DummyStep_H

#include "InputStep.h"

#include "../base/DPBuffer.h"

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
  Upsample(const common::ParameterSet&, const string& prefix);

  /// Constructor that gets the settings directly.
  Upsample(const std::string& name, unsigned int time_step, bool update_uvw);

  virtual ~Upsample();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

 private:
  const std::string name_;
  const unsigned int time_step_;
  const bool update_uvw_;

  std::vector<base::DPBuffer> prev_buffers_;
  std::vector<base::DPBuffer> buffers_;
  unsigned int first_to_flush_;
  std::unique_ptr<base::UVWCalculator> uvw_calculator_;

  common::NSTimer timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
