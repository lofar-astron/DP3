// Counter.h: DPPP step class to count flags
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to count flags
/// @author Ger van Diepen

#ifndef DPPP_COUNTER_H
#define DPPP_COUNTER_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/FlagCounter.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DPPP step class to count flags

/// This class is a Step class counting the number of flags per
/// baseline and channel.
/// It can be used for test purposes to know how many flags have been
/// set by the previous steps.

class Counter : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Counter(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~Counter();

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the flag counts.
  virtual void showCounts(std::ostream&) const;

 private:
  std::string itsName;
  bool itsFlagData;
  unsigned int itsCount;
  bool itsSaveToJson;
  std::string itsJsonFilename;
  base::FlagCounter itsFlagCounter;
};

}  // namespace steps
}  // namespace dp3

#endif
