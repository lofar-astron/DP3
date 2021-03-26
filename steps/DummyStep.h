// DummyStep.h: DPPP step class to DummyStep visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to DummyStep visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DPPP_DummyStep_H
#define DPPP_DummyStep_H

#include "InputStep.h"

#include "../base/DPBuffer.h"

#include <utility>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DPPP step class to DummyStep visibilities from a source model

/// This class is an empty Step subclass to use as implementation template

class DummyStep : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DummyStep(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~DummyStep();

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

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBuffer;

  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
