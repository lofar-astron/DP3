// AOFlaggerStep.h: DPPP step class to flag data using rficonsole's
// functionality
//
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
/// @file
/// @author Andre Offringa, Ger van Diepen

#ifndef DP3_STEPS_AOFLAGGERSTEP_H_
#define DP3_STEPS_AOFLAGGERSTEP_H_

#include "../steps/InputStep.h"

#include <dp3/base/DPBuffer.h>
#include "../base/FlagCounter.h"

#include <memory>
#include <mutex>

#include <aoflagger.h>

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DPPP step class to flag using aoflagger's functionality

/// This class is a Step class flagging data points based on the
/// aoflagger library written by Andre Offringa.
/// See the following papers for background information:
/// <ul>
/// <li> Post-correlation radio frequency interference classification
///      methods -- http://arxiv.org/abs/1002.1957
/// <li> A LOFAR RFI detection pipeline and its first results
///      -- http://arxiv.org/abs/1007.2089
/// </ul>
///
/// When a correlation is flagged, all correlations for that data point
/// are flagged. It is possible to specify which correlations have to be
/// taken into account when flagging. Using, say, only XX may boost
/// performance with a factor 4, but miss points to be flagged.
/// It is also possible to specify the order in which the correlations
/// have to be tested.
///
/// It is possible to flag specific baselines only using a selection on
/// baseline length.
/// <br>Furthermore it is possible to only flag the autocorrelations and
/// apply the resulting flags to the crosscorrelations, possibly selected
/// on baseline length.

class AOFlaggerStep : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  AOFlaggerStep(const common::ParameterSet&, const std::string& prefix);

  ~AOFlaggerStep() override;

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField;
  }

  common::Fields getProvidedFields() const override { return kFlagsField; }

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Write the statistics into the MS.
  void addToMS(const std::string& msName) override;

  /// Update the general info.
  /// It is used to adjust the parms if needed.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the flagger counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

 private:
  /// Flag all baselines in the time window (using OpenMP to parallellize).
  /// Process the buffers in the next step.
  void flag(unsigned int rightOverlap);

  /// Format a number as kB, MB, etc.
  static void formatBytes(std::ostream&, double);

  /// Flag a single baseline using the rfistrategy.
  void flagBaseline(unsigned int leftOverlap, unsigned int windowSize,
                    unsigned int rightOverlap, unsigned int bl,
                    base::FlagCounter& counter, aoflagger::Strategy& strategy,
                    aoflagger::QualityStatistics& rfiStats);

  /// Add the flags to the statistics.
  void addStats(aoflagger::QualityStatistics& rfiStats,
                const aoflagger::ImageSet& values,
                const aoflagger::FlagMask& rfiMask,
                const aoflagger::FlagMask& origMask, int bl);

  std::string name_;
  unsigned int buffer_index_;
  unsigned int n_times_;
  std::string strategy_name_;
  unsigned int window_size_;
  unsigned int overlap_;  ///< extra time slots on both sides
  double overlap_percentage_;
  double memory_;  ///< Usable memory in GBytes
  double memory_percentage_;
  double memory_needed_;  ///< Memory needed for data/flags
  bool flag_auto_correlations_;
  bool collect_statistics_;
  std::vector<std::unique_ptr<base::DPBuffer>> buffer_;
  base::FlagCounter flag_counter_;
  common::NSTimer total_timer_;
  common::NSTimer quality_timer_;
  common::NSTimer compute_timer_;
  double move_time_;   ///< data move timer (sum of all threads)
  double flag_time_;   ///< flag timer (sum of all threads)
  double stats_time_;  ///< quality timer (sum of all threads)
  std::vector<double> frequencies_;
  aoflagger::AOFlagger aoflagger_;
  std::mutex mutex_;
  bool has_qstats_;
  aoflagger::QualityStatistics qstats_;
};

}  // namespace steps
}  // namespace dp3

#endif
