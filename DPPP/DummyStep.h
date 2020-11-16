// DummyStep.h: DPPP step class to DummyStep visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to DummyStep visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DPPP_DummyStep_H
#define DPPP_DummyStep_H

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

namespace DP3 {

class ParameterSet;

namespace DPPP {
/// @brief DPPP step class to DummyStep visibilities from a source model

/// This class is an empty DPStep subclass to use as implementation template

class DummyStep : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DummyStep(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~DummyStep();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  DPInput* itsInput;
  string itsName;
  DPBuffer itsBuffer;

  NSTimer itsTimer;
};

}  // namespace DPPP
}  // namespace DP3

#endif
