// OnePredict.h: DP3 step class to predict visibilities from a source model
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to predict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_ONEPREDICT_H_
#define DP3_STEPS_ONEPREDICT_H_

#include <atomic>

#include <xtensor/xtensor.hpp>

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "ResultStep.h"
#include "base/DP3.h"
#include "base/DPBuffer.h"

#include "../base/ModelComponent.h"
#include "../base/PredictBuffer.h"
#include "../base/PredictModel.h"

#include "../model/Patch.h"
#include "../model/SourceDBUtil.h"

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// @brief Step class that predicts visibilities with optionally beam.
/// The Predict class uses one or more instances of this class for predicting
/// data with different regular shapes.
class OnePredict : public ModelDataStep {
 public:
  /**
   * Constructs the object.
   * @param parset Parameter set with settings for the step.
   * @param prefix Prefix for reading settings from 'parset'.
   * @param sourceList Direction names. If empty, obtain sources from the parset
   */
  OnePredict(const common::ParameterSet&, const std::string& prefix,
             const std::vector<std::string>& source_patterns);

  ~OnePredict() override;

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

  /// When multiple OnePredict steps are running in parallel from multiple
  /// threads, they require synchronisation. This is done with this mutex.
  /// When multiple Predicts steps run serially
  /// (like currently in H5ParmPredict), this function should not be called, as
  /// otherwise they will synchronize needlessly.
  void SetThreadData(std::mutex* measures_mutex) {
    measures_mutex_ = measures_mutex;
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
  void init(const common::ParameterSet&, const std::string& prefix,
            const std::vector<std::string>& sourcePatterns);

  void initializeThreadData();
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);

  void addBeamToData(const model::Patch& patch, size_t buffer_index,
                     aocommon::xt::UTensor<std::complex<double>, 3>& model_data,
                     double time, bool update_beam, size_t thread,
                     aocommon::xt::UTensor<std::complex<double>, 3>& data0,
                     bool stokesIOnly);

  void addBeamToDataRange(
      const model::Patch& patch,
      aocommon::xt::UTensor<std::complex<double>, 3>& model_data, double time,
      size_t thread, aocommon::xt::UTensor<std::complex<double>, 3>& data0,
      const std::pair<size_t, size_t>& baseline_range,
      const std::pair<size_t, size_t>& station_range,
      aocommon::Barrier& barrier, bool stokesIOnly);

  void PredictWithSourceParallelization(base::DPBuffer::DataType& destination,
                                        double time);
  void PredictSourceRange(
      aocommon::xt::UTensor<std::complex<double>, 3>& result, size_t start,
      size_t end, size_t thread_index, std::mutex& mutex, double time,
      bool update_beam);

  /// Assigns @p buffer to @p destination. If @c stokes_i_only_ is set,
  /// only the first and last correlations (e.g. XX and YY) are copied.
  void CopyPredictBufferToData(
      base::DPBuffer::DataType& destination,
      const aocommon::xt::UTensor<std::complex<double>, 3>& buffer);

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
  std::shared_ptr<std::vector<base::PredictBuffer>> predict_buffers_;
  everybeam::BeamMode beam_mode_ = everybeam::BeamMode::kNone;
  everybeam::ElementResponseModel element_response_model_ =
      everybeam::ElementResponseModel::kDefault;
  std::vector<casacore::MeasFrame> meas_frame_;
  std::vector<casacore::MDirection::Convert> meas_convertors_;
  std::shared_ptr<everybeam::telescope::Telescope> telescope_;

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
  /**
   * The total time [µs] of the apply beam phase.
   *
   * Similar to @ref predict_time_.
   */
  std::atomic<int64_t> apply_beam_time_ = 0;

  std::mutex* measures_mutex_;
  std::mutex mutex_;
};

}  // namespace steps
}  // namespace dp3

#endif
