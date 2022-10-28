// OnePredict.h: DP3 step class to predict visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to predict visibilities from a source model
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_ONEPREDICT_H_
#define DP3_STEPS_ONEPREDICT_H_

#include <atomic>

#include "ApplyBeam.h"
#include "ApplyCal.h"
#include "ResultStep.h"
#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>
#include "../base/ModelComponent.h"
#include "../base/Patch.h"
#include "../base/PredictBuffer.h"
#include "../base/SourceDBUtil.h"

namespace aocommon {
class ThreadPool;
}  // namespace aocommon

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

  virtual ~OnePredict();

  common::Fields getRequiredFields() const override {
    common::Fields fields = kUvwField;
    if (operation_ == "add" || operation_ == "subtract") {
      fields |= kDataField;
    }
    if (apply_cal_step_) {
      fields |= base::DP3::GetChainRequiredFields(apply_cal_step_);
    }
    return fields;
  }

  common::Fields getProvidedFields() const override {
    // When operation_ == "replace", the output of apply_cal_step_ is passed to
    // the next step. In all other cases, only kDataField changes.
    common::Fields fields = kDataField;
    if (operation_ == "replace" && apply_cal_step_) {
      std::shared_ptr<Step> step = apply_cal_step_;
      do {
        fields |= step->getProvidedFields();
        step = step->getNextStep();
      } while (step);
    }
    return fields;
  }

  /// Set the applycal substep
  void SetApplyCal(const common::ParameterSet&, const std::string& prefix);

  /// Set the operation type
  void SetOperation(const std::string& type);

  /// When multiple OnePredict steps are running in parallel from multiple
  /// threads, they require synchronisation. This is done with these two
  /// synchronisation structures. When multiple Predicts steps run serially
  /// (like currently in H5ParmPredict), this function should not be called, as
  /// otherwise they will synchronize needlessly.
  ///
  /// It is also possible to make the predict steps share the same threadpool
  /// without further synchronisation, by setting measures_mutex to nullptr.
  void SetThreadData(aocommon::ThreadPool& pool, std::mutex* measures_mutex) {
    thread_pool_ = &pool;
    measures_mutex_ = measures_mutex;
  }

  void SetPredictBuffer(std::shared_ptr<base::PredictBuffer> predict_buffer) {
    predict_buffer_ = std::move(predict_buffer);
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Prepare the sources
  void setSources(const std::vector<string>& sourcePatterns);

  /// Return the direction of the first patch
  base::Direction GetFirstDirection() const override;

 private:
  /// The actual constructor
  void init(const common::ParameterSet&, const std::string& prefix,
            const std::vector<std::string>& sourcePatterns);

  void initializeThreadData();
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);
  void addBeamToData(std::shared_ptr<const base::Patch> patch, double time,
                     size_t thread, size_t nBeamValues,
                     std::complex<double>* data0, bool stokesIOnly);
  void addBeamToData(std::shared_ptr<const base::Patch> patch, double time,
                     size_t thread, size_t nBeamValues,
                     std::complex<double>* data0,
                     const std::pair<size_t, size_t>& baseline_range,
                     const std::pair<size_t, size_t>& station_range,
                     aocommon::Barrier& barrier, bool stokesIOnly);

  std::string name_;
  base::DPBuffer buffer_;
  std::string source_db_name_;
  bool correct_freq_smearing_{false};
  std::string operation_;
  bool apply_beam_{false};
  bool use_channel_freq_{false};
  bool one_beam_per_patch_{false};
  bool thread_over_baselines_{false};
  /// If two sources are closer together than given by this setting, they
  /// will be grouped into one patch. Value is in arcsec; zero means don't
  /// group.
  double beam_proximity_limit_{false};
  bool stokes_i_only_{false};
  bool any_orientation_is_absolute_{false};  ///< Any of the Gaussian sources
                                             ///< has absolute orientation
  base::Direction phase_ref_;
  bool moving_phase_ref_{false};

  std::shared_ptr<ApplyCal> apply_cal_step_;  ///< Optional ApplyCal sub step
  std::shared_ptr<ResultStep> result_step_;   ///< Catches results from ApplyCal

  unsigned int debug_level_{0};

  std::vector<std::pair<size_t, size_t>> baselines_;

  /// Vector containing info on converting baseline uvw to station uvw
  std::vector<int> uvw_split_index_;

  /// UVW coordinates per station (3 coordinates per station)
  casacore::Matrix<double> station_uvw_;

  /// The info needed to calculate the station beams.
  std::shared_ptr<base::PredictBuffer> predict_buffer_;
  everybeam::CorrectionMode beam_mode_{everybeam::CorrectionMode::kNone};
  everybeam::ElementResponseModel element_response_model_{
      everybeam::ElementResponseModel::kDefault};
  std::vector<casacore::MeasFrame> meas_frame_;
  std::vector<casacore::MDirection::Convert> meas_convertors_;
  std::shared_ptr<everybeam::telescope::Telescope> telescope_;

  std::string direction_str_;  ///< Definition of patches, to pass to applycal
  std::vector<std::shared_ptr<base::Patch>> patch_list_;

  std::vector<std::pair<std::shared_ptr<base::ModelComponent>,
                        std::shared_ptr<base::Patch>>>
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
  std::atomic<int64_t> predict_time_{0};
  /**
   * The total time [µs] of the apply beam phase.
   *
   * Similar to @ref predict_time_.
   */
  std::atomic<int64_t> apply_beam_time_{0};

  aocommon::ThreadPool* thread_pool_;
  std::mutex* measures_mutex_;
  std::mutex mutex_;
};

}  // namespace steps
}  // namespace dp3

#endif
