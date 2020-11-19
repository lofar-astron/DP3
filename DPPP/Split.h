// Split.h: DPPP step class to Split visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to Split visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DPPP_Split_H
#define DPPP_Split_H

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

namespace DP3 {

class ParameterSet;

namespace DPPP {
/// @brief DPPP step class to Split visibilities from a source model

/// This class is an empty DPStep subclass to use as implementation template

class Split : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Split(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~Split();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  virtual void addToMS(const string&);

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  string itsName;

  std::vector<std::string> itsReplaceParms;  ///< The names of the parameters
                                             ///< that differ along the substeps
  std::vector<DPStep::ShPtr> itsSubsteps;
  bool itsAddedToMS;  ///< Used in addToMS to prevent recursion
};

}  // namespace DPPP
}  // namespace DP3

#endif
