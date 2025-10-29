// FastPredict.h: DP3 step class to predict visibilities from a source model
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to predict visibilities from a source model

#ifndef DP3_STEPS_FASTPREDICT_H_
#define DP3_STEPS_FASTPREDICT_H_

#ifdef USE_FAST_PREDICT

#include <atomic>

#include <xtensor/xtensor.hpp>

#include "ApplyCal.h"
#include "ResultStep.h"
#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>

#include "../base/ModelComponent.h"
#include "../base/PredictBuffer.h"

#include "../model/Patch.h"

#include <predict/PredictPlan.h>
#include <predict/PredictPlanExecCPU.h>
#include <predict/Predict.h>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief Step class that predicts visibilities, optionally with beam.
/// The Predict class uses one or more instances of this class for predicting
/// data with different regular shapes.
class FastPredict : public ModelDataStep {
 public:
  /**
   * @param sourceList Direction names. If empty, obtain sources from the parset
   */
  FastPredict(const common::ParameterSet&, const std::string& prefix,
              const std::vector<std::string>& source_patterns);

  ~FastPredict() override;

  common::Fields getRequiredFields() const override {
    common::Fields fields = kUvwField;
    if (operation_ == Operation::kAdd || operation_ == Operation::kSubtract) {
      fields |= kDataField;
    }
    if (apply_cal_step_) {
      fields |= base::GetChainRequiredFields(apply_cal_step_);
    }
    return fields;
  }

  common::Fields getProvidedFields() const override {
    // When operation_ == "replace", the output of apply_cal_step_ is passed to
    // the next step. For all other operations, the predicted visibility data
    // is the only provided output.
    common::Fields fields;

    if (output_data_name_.empty()) {
      // The predicted visibilities go to the main data buffer
      fields |= kDataField;
    } else {
      // TODO(AST-1241): Handle these dependencies using Fields.
    }

    if (operation_ == Operation::kReplace && apply_cal_step_) {
      std::shared_ptr<Step> step = apply_cal_step_;
      do {
        fields |= step->getProvidedFields();
        step = step->getNextStep();
      } while (step);
    }
    return fields;
  }

  /// Set the ApplyCal substep and connect it to a ResultStep
  void SetApplyCal(const common::ParameterSet&, const std::string& prefix);

  /// Set the operation type
  void SetOperation(const std::string& type);

  void SetThreadData(std::mutex*) {
    // FastPredict does not require synchronisation, this is done in the
    // RD/Predict module itself. This method is need to satisfy the
    // steps/Predict.h interface.
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Prepare the sources
  void setSources(const std::vector<std::string>& sourcePatterns);

  /// Return the direction of the first patch
  base::Direction GetFirstDirection() const override;

 private:
  enum class Operation { kReplace, kAdd, kSubtract };

  /// The actual constructor
  void Init(const common::ParameterSet&, const std::string& prefix,
            const std::vector<std::string>& sourcePatterns);

  void InitializePlan();

  void RunPlan(base::DPBuffer::DataType& destination, double time);

  void CopyPredictBufferToData(
      base::DPBuffer::DataType& destination,
      const xt::xtensor<double, 4, xt::layout_type::row_major>& buffer);

  std::string name_;
  /// Stores the input data if the operation is add or subtract.
  /// Using a member instead of a local variable avoids allocating memory
  /// for the input data in each process() call.
  xt::xtensor<std::complex<float>, 3> input_data_;
  std::string source_db_name_;
  bool correct_time_smearing_ = false;
  bool correct_freq_smearing_ = false;
  Operation operation_;
  std::string output_data_name_;
  bool apply_beam_ = false;
  std::string coefficients_path_;
  bool use_channel_freq_ = false;
  bool one_beam_per_patch_ = false;
  bool thread_over_baselines_ = false;
  /// If two sources are closer together than given by this setting, they
  /// will be grouped into one patch. Value is in arcsec; zero means don't
  /// group.
  double beam_proximity_limit_ = 0.0;
  double beam_evaluation_interval_ = 0.0;
  double previous_beam_time_ = 0.0;
  bool stokes_i_only_ = false;
  bool any_orientation_is_absolute_ = false;  ///< Any of the Gaussian sources
                                              ///< has absolute orientation
  base::Direction phase_ref_;
  bool moving_phase_ref_ = false;

  std::shared_ptr<ApplyCal> apply_cal_step_;  ///< Optional ApplyCal sub step
  std::shared_ptr<ResultStep> result_step_;   ///< Catches results from ApplyCal

  unsigned int debug_level_ = 0;

  std::vector<double> scaled_ncp_uvw_;

  std::vector<std::pair<size_t, size_t>> baselines_;

  /// Vector containing info on converting baseline uvw to station uvw
  std::vector<int> uvw_split_index_;

  /// UVW coordinates per station (3 coordinates per station)
  xt::xtensor<double, 2> station_uvw_;

  /// The info needed to calculate the station beams.
  everybeam::CorrectionMode beam_mode_ = everybeam::CorrectionMode::kNone;
  everybeam::ElementResponseModel element_response_model_ =
      everybeam::ElementResponseModel::kDefault;
  casacore::MeasFrame meas_frame_;
  casacore::MDirection::Convert meas_converter_;
  std::shared_ptr<everybeam::telescope::Telescope> telescope_;
  predict::PredictPlan predict_plan_;
  predict::Predict predict_;
  std::unique_ptr<predict::PredictPlanExecCPU> predict_plan_exec_;

  std::string direction_str_;  ///< Definition of patches, to pass to applycal
  std::vector<std::shared_ptr<model::Patch>> patch_list_;

  std::vector<std::pair<std::shared_ptr<base::ModelComponent>,
                        std::shared_ptr<model::Patch>>>
      source_list_;

  common::NSTimer timer_;

  /**
   * The total time [µs] of the prediction phase.
   *
   * The total prediction time is the sum of the execution time of all threads
   * in the predict phase. This time can be more than the time measured by
   * @ref timer_.
   *
   * @note Atomic floating point operations require C++20. The timer has a
   * microsecond resolution. So by multiplying by 1e6 the result can be
   * lossless stored in an integral.
   */
  std::atomic<int64_t> predict_time_ = 0;

  std::mutex* measures_mutex_;
  std::mutex mutex_;
};

}  // namespace steps
}  // namespace dp3

#endif  // USE_FAST_PREDICT
#endif  // DP3_BASE_FASTPREDICT_H