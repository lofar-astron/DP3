// Upsample.h: DPPP step class to upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to Upsample visibilities
/// @author Tammo Jan Dijkema

#ifndef DPPP_DummyStep_H
#define DPPP_DummyStep_H

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

namespace DP3 {

class ParameterSet;

namespace DPPP {
/// @brief DPPP step class to Upsample visibilities

/// This class is an empty DPStep subclass to use as implementation template

class Upsample : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Upsample(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~Upsample();

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

 private:
  string itsName;
  double itsOldTimeInterval;
  unsigned int itsTimeStep;

  std::vector<DPBuffer> itsPrevBuffers;
  std::vector<DPBuffer> itsBuffers;
  unsigned int itsFirstToFlush;

  NSTimer itsTimer;
};

}  // namespace DPPP
}  // namespace DP3

#endif
