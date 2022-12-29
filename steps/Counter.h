// Counter.h: DP3 step class to count flags
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to count flags
/// @author Ger van Diepen

#ifndef DP3_STEPS_COUNTER_H_
#define DP3_STEPS_COUNTER_H_

#include <dp3/steps/Step.h>

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"

namespace dp3 {
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
  explicit Counter(const common::ParameterSet&, const string& prefix);

  ~Counter() override;

  common::Fields getRequiredFields() const override {
    if (flag_data_)
      // Visibility data must be read if needed, so NaNs are flagged.
      return kDataField | kFlagsField;
    else
      return kFlagsField;
  }

  common::Fields getProvidedFields() const override { return {}; }

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the flag counts.
  void showCounts(std::ostream&) const override;

 private:
  std::string name_;
  bool flag_data_;
  unsigned int count_;
  bool save_to_json_;
  std::string json_filename_;
  base::FlagCounter flag_counter_;
};

}  // namespace steps
}  // namespace dp3

#endif
