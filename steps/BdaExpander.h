// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for expanding BDA data to regular data.
/// @author Chiara Salvoni

#ifndef DP3_STEPS_BDAEXPANDER_H_
#define DP3_STEPS_BDAEXPANDER_H_

#include <map>
#include <utility>
#include <vector>

#include <dp3/base/BdaBuffer.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/steps/Step.h>

#include "../common/Timer.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step that expands BDA data in BdaBuffers to regular data in
/// DPBuffers.

/// This class expands bda averaged data to regular data. The averaged data is
/// copied throughout all the timeslots which have been averaged originally.
/// This step doesn't add information, it is just used for backward
/// compatibility to ensure the correct workflow for processing steps which are
/// not yet able to handle BDA data

class BdaExpander : public Step {
 public:
  explicit BdaExpander(const std::string& prefix);

  common::Fields getRequiredFields() const override { return kUvwField; }

  common::Fields getProvidedFields() const override {
    // The BdaExpander only provides fields that are already present in the
    // input BdaBuffer -> return an empty Fields object.
    return {};
  }

  /// Process the data.
  /// Reads the data from a BdaBuffer and fills an internal vector of DPBuffer.
  /// Gives regular buffer to the next step once all the baselines are
  /// available.
  bool process(std::unique_ptr<base::BdaBuffer>) override;

  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

  bool accepts(MsType dt) const override { return dt == MsType::kBda; }

 private:
  /// Helper function to copy the data from BdaBuffer to DPBuffer.
  /// @param bda_row The row in the BdaBuffer to copy.
  /// @param buf_out The regular buffer in which the BdaBuffer row is copied.
  /// @param current_bl The baseline number relative to the BDA_row.
  void CopyData(const base::BdaBuffer& bda_buffer,
                const base::BdaBuffer::Row& bda_row,
                std::unique_ptr<dp3::base::DPBuffer>& buf_out,
                int time_averaging_factor = 1);

  struct RegularBufferElement {
    RegularBufferElement(size_t n_baselines, unsigned int n_corr,
                         unsigned int n_chan, double current_time,
                         double current_exposure);

    std::vector<bool> baseline;
    std::unique_ptr<dp3::base::DPBuffer> regular_buffer;
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
