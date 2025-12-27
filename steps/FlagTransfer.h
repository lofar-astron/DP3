// FlagTransfer.h: DP3 step class to transfer flags
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to transfer flags from a lower to a MS with higher
/// time/freq resolution
/// @author Mick Veldhuis

#ifndef DP3_STEPS_FLAG_TRANSFER_H_
#define DP3_STEPS_FLAG_TRANSFER_H_

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableIter.h>
#include "base/DPBuffer.h"
#include "steps/Step.h"
#include "Filter.h"
#include "ResultStep.h"
#include "../common/Timer.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

/// \brief DP3 step class to transfer flags from a lower to a MS with higher
/// time/freq resolution

/// This class is a Step class for extrapolating flags of lower time and/or
/// frequency resolution data to a higher resolution MS. This is achieved by
/// buffering the low-resolution flags for a number of time slots and channels
/// depending on the time interval and channel width of the source and target
/// MS, and subsequently writing these flags to DPBuffer.

/// Additionally, a Filter sub-step can be used to transfer flags to a MS that,
/// e.g., contains more baselines than the source MS.

class FlagTransfer final : public Step {
 public:
  /// Construct the object using parameter values from the parset, using the
  /// given prefix.
  explicit FlagTransfer(const common::ParameterSet& parameter_set,
                        const std::string& prefix);

  common::Fields getRequiredFields() const final {
    return filter_step_->getRequiredFields();
  }

  common::Fields getProvidedFields() const final { return kFlagsField; }

  bool process(std::unique_ptr<base::DPBuffer> buffer) final;

  void finish() final;

  void updateInfo(const base::DPInfo&) final;

  void show(std::ostream&) const final;

  void showTimings(std::ostream&, double duration) const final;

 private:
  /// Read FLAG column from the time step pointed to by ms_iterator_
  void ReadSourceMsFlags();

  std::string name_;
  common::NSTimer timer_;
  std::string source_ms_path_;
  casacore::MeasurementSet ms_;
  casacore::TableIterator ms_iterator_;

  /// Flags from the source MeasurementSet (ms_) for a single time step,
  /// with shape ( baseline x channel x correlation )
  base::DPBuffer::FlagsType flags_;

  /// Remainder of time steps for which to hold flags_
  std::size_t timestep_counter_;

  /// Ratios of the time interval between the source and target MeasurementSet
  std::size_t time_averaging_factor_;

  /// Time interval of the data in the source MeasurementSet
  double time_interval_;

  /// Upper edges of source MS frequency bins
  std::vector<double> source_channel_upper_edges_;

  std::shared_ptr<Step> filter_step_;
  std::shared_ptr<ResultStep> result_step_;
};

}  // namespace steps
}  // namespace dp3

#endif
