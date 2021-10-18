// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for expanding BDA data to regular data.
/// @author Chiara Salvoni

#ifndef DPPP_BDAEXPANDER_H
#define DPPP_BDAEXPANDER_H

#include "InputStep.h"
#include "Step.h"

#include "../base/BDABuffer.h"
#include "../base/DPBuffer.h"

#include <map>
#include <utility>
#include <vector>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DP3 step that expands BDA data in BDABuffers to regular data in
/// DPBuffers.

/// This class expands bda averaged data to regular data. The averaged data is
/// copied throughout all the timeslots which have been averaged originally.
/// This step doesn't add information, it is just used for backward
/// compatibility to ensure the correct workflow for processing steps which are
/// not yet able to handle BDA data

class BDAExpander : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  BDAExpander(const string &prefix);

  ~BDAExpander() override;

  /// Process the data.
  /// Reads the data from a BDABuffer and fills an internal vector of DPBuffer.
  /// Gives regular buffer to the next step once all the baselines are
  /// available.
  bool process(std::unique_ptr<base::BDABuffer>) override;

  void finish() override;

  void updateInfo(const base::DPInfo &) override;

  void show(std::ostream &) const override;

  void showTimings(std::ostream &, double duration) const override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

 private:
  /// Helper function to copy the data from BDABuffer to DPBuffer
  /// BDA_row: the row in the BDABuffer to copy
  /// bufOut: the regular buffer in which the BDABuffer row is copied
  /// current_bl: the baseline number relative to the BDA_row
  void CopyData(const dp3::base::BDABuffer::Row &bda_row,
                dp3::base::DPBuffer &buf_out, unsigned int current_bl,
                float time_averaging_factor = 1.0);

  class RegularBufferElement {
   public:
    RegularBufferElement() {
      throw std::runtime_error("Default constructor should not be called");
    }
    RegularBufferElement(int n_baseline, unsigned int n_corr,
                         unsigned int n_chan, double current_time,
                         double current_exposure);

    std::vector<bool> baseline_;
    dp3::base::DPBuffer regular_buffer;
  };

  std::map<unsigned int, RegularBufferElement> RB_elements;

  /**
   * This variable has size nbaselines x nchan (without averaging). For each
   * baseline and each non-averaged output channel index, it contains the
   * corresponding BDA input channel index. Examples, when nchan is 4:
   * - No frequency averaging: [0 1 2 3]
   * - Two channels averaged: [0 0 1 1]
   * - All four channels averaged: [0 0 0 0]
   */
  std::vector<std::vector<int>> channels_mapping_;

  unsigned int next_time_slot_to_process_;

  common::NSTimer timer_;
  std::string step_name_;
};

}  // namespace steps
}  // namespace dp3

#endif
