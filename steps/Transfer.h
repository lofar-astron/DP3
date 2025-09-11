// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @brief DP3 step class to transfer data and flags from a lower to a MS with
/// higher time/freq resolution
/// @author Mick Veldhuis

#ifndef DP3_STEPS_DATA_TRANSFER_H_
#define DP3_STEPS_DATA_TRANSFER_H_

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableIter.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/steps/Step.h>
#include "Filter.h"
#include "ResultStep.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step class to transfer visbility data and flags from a lower to a
/// MS with higher time/freq resolution

/// This class is a Step class for extrapolating visibilities and flags recorded
/// at lower time and/or frequency resolution to a higher resolution MS. This is
/// achieved by buffering the low-resolution data for a number of time slots and
/// channels depending on the time interval and channel width of the source and
/// target MS, and subsequently writing these to DPBuffer.

/// Additionally, a Filter sub-step can be used to transfer data to a MS that,
/// e.g., contains more baselines than the source MS.

class Transfer final : public Step {
 public:
  /// Construct the object using parameter values from the parset, using the
  /// given prefix.
  explicit Transfer(const common::ParameterSet& parameter_set,
                    const std::string& prefix);

  common::Fields getRequiredFields() const final {
    return filter_step_->getRequiredFields();
  }

  common::Fields getProvidedFields() const final {
    common::Fields provided_fields{};
    provided_fields |= (transfer_data_ && output_buffer_name_.empty())
                           ? kDataField
                           : common::Fields();
    provided_fields |= transfer_flags_ ? kFlagsField : common::Fields();
    return provided_fields;
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) final;

  void finish() final;

  void updateInfo(const base::DPInfo&) final;

  void show(std::ostream&) const final;

  void showTimings(std::ostream&, double duration) const final;

 private:
  /// Read the DATA column from the time step pointed to by ms_iterator_
  void ReadSourceMsVisibilities();

  /// Read the FLAG column from the time step pointed to by ms_iterator_
  void ReadSourceMsFlags();

  /// Transfers a single time slot of data from a low to higher-resolution MS
  template <typename T>
  void TransferSingleTimeSlot(
      T source, T& target,
      const std::vector<common::rownr_t>& target_row_numbers,
      const std::vector<double>& target_frequencies) const;

  std::string name_;
  common::NSTimer timer_;
  std::string source_ms_path_;
  std::string source_data_column_;
  casacore::MeasurementSet ms_;
  casacore::TableIterator ms_iterator_;

  /// Visibility and flags from the source MeasurementSet (ms_) for a single
  /// time step, with shape ( baseline x channel x correlation )
  /// @{
  base::DPBuffer::DataType data_;
  base::DPBuffer::FlagsType flags_;
  /// @}

  /// Determines what to transfer from the low high-resolution MS
  /// @{
  bool transfer_data_;
  bool transfer_flags_;
  /// @}

  /// Whether to put the data in a named buffer
  std::string output_buffer_name_;

  /// Remainder of time steps for which to hold data_/flags_
  std::size_t timestep_counter_ = 0;

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
