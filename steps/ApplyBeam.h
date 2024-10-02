// ApplyBeam.h: DP3 step class to ApplyBeam visibilities from a source model
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply the beam model (optionally inverted)
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_APPLYBEAM_H_
#define DP3_STEPS_APPLYBEAM_H_

#include "InputStep.h"

#include <dp3/base/DPBuffer.h>

#include <EveryBeam/telescope/telescope.h>

#include <aocommon/matrix2x2.h>
#include <aocommon/barrier.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/measures/Measures/MDirection.h>

#include "../common/ParameterSet.h"

namespace dp3 {
namespace steps {

/// Computes full 2x2 Jones beam matrices using EveryBeam.
size_t ComputeBeam(const base::DPInfo& info, double time,
                   const everybeam::vector3r_t& srcdir,
                   const everybeam::telescope::Telescope* telescope,
                   aocommon::MC2x2* beam_values, bool invert,
                   everybeam::CorrectionMode mode, std::mutex* mutex,
                   const std::vector<size_t>& skip_station_indices);

/// Computes the array factor scalar values.
size_t ComputeArrayFactor(const base::DPInfo& info, double time,
                          const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          everybeam::complex_t* beam_values, bool invert,
                          std::mutex* mutex,
                          const std::vector<size_t>& skip_station_indices);

/**
 * Corrects the values in @p model_data with the precomputed full Jones beam
 * values, and adds the corrected model data to @p data0.
 * @param data0 An array of n_baselines x n_channels x n_correlations
 * (with n_correlations the fastest changing) containing the data.
 * @param model_data Array of same shape as data0; the corrected values are
 * added to these data.
 * @param beam_values Array of n_atenna x n_channels containing the
 * pre-calculated beam matrices.
 */
void ApplyBeamToDataAndAdd(const base::DPInfo& info, size_t n_stations,
                           std::complex<double>* data0,
                           std::complex<double>* model_data, float* weight0,
                           const aocommon::MC2x2* beam_values,
                           bool doUpdateWeights);

/**
 * Corrects the values in @p model_data with the precomputed scalar array
 * factors. Note that unlike @ref ApplyBeamToDataAndAdd(), the data values
 * are assumed to be Stokes I values only. This is used for the optimization
 * when the sky model is unpolarized and only the array factor is applied.
 * @param data0 An array of n_baselines x n_channels
 * (with n_channels the fastest changing) containing the data.
 * @param beam_values Array of n_atenna x n_channels containing the
 * pre-calculated scalar beam values.
 */
void ApplyArrayFactor(const base::DPInfo& info, size_t n_stations,
                      std::complex<double>* data0,
                      const everybeam::complex_t* beam_values);

/// \brief DP3 step class to ApplyBeam visibilities from a source model

/// This class is a Step class to apply the beam model, optionally inverted.
/// The input MeasurementSet it operates on, must have the LOFAR subtables
/// defining the station layout and tiles/dipoles used.

class ApplyBeam final : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ApplyBeam(const common::ParameterSet&, const string& prefix,
            bool substep = false);

  common::Fields getRequiredFields() const override {
    common::Fields fields = kDataField;
    if (itsUpdateWeights) fields |= kWeightsField;
    return fields;
  }

  common::Fields getProvidedFields() const override {
    common::Fields fields = kDataField;
    if (itsUpdateWeights) fields |= kWeightsField;
    return fields;
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override {
    if (use_model_data_) {
      return ProcessModelData(std::move(buffer));
    } else {
      return ProcessData(std::move(buffer));
    }
  }

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  bool invert() { return itsInvert; }

  /**
   * Calculate and apply the beam for processing when
   * parallelizing over baselines. Because the beam is a per-antenna effect,
   * this requires synchronisation, which is performed with the provided
   * barrier.
   */
  static void ApplyBaselineBasedBeam(
      const base::DPInfo& info, double time, std::complex<double>* data0,
      float* weight0, const everybeam::vector3r_t& srcdir,
      const everybeam::telescope::Telescope* telescope,
      aocommon::MC2x2* beam_values,
      const std::pair<size_t, size_t>& baseline_range,
      const std::pair<size_t, size_t>& station_range,
      aocommon::Barrier& barrier, bool invert, everybeam::CorrectionMode mode,
      bool do_update_weights = false, std::mutex* mutex = nullptr,
      const std::vector<size_t>& skip_station_indices = std::vector<size_t>());

  /**
   * Like @ref ApplyBaselineBasedBeam(), but for array factor only.
   */
  static void ApplyBaselineBasedArrayFactor(
      const base::DPInfo& info, double time, std::complex<double>* data0,
      const everybeam::vector3r_t& srcdir,
      const everybeam::telescope::Telescope* telescope,
      everybeam::complex_t* beam_values,
      const std::pair<size_t, size_t>& baseline_range,
      const std::pair<size_t, size_t>& station_range,
      aocommon::Barrier& barrier, bool invert, everybeam::CorrectionMode mode,
      std::mutex* mutex = nullptr,
      const std::vector<size_t>& skip_station_indices = std::vector<size_t>());

 private:
  everybeam::vector3r_t dir2Itrf(const casacore::MDirection& dir,
                                 casacore::MDirection::Convert& measConverter);
  bool ProcessData(std::unique_ptr<base::DPBuffer> buffer);
  bool ProcessModelData(std::unique_ptr<base::DPBuffer> buffer);

  string itsName;
  bool itsInvert;
  bool itsUpdateWeights;
  std::vector<string> itsDirectionStr;
  casacore::MDirection itsDirection;
  bool itsUseChannelFreq;
  std::vector<std::string> itsSkipStationNames;
  std::vector<size_t> itsSkipStationIndices;
  everybeam::CorrectionMode itsMode;
  everybeam::ElementResponseModel itsElementResponseModel;

  /// If a beam has already been applied before running this step, that beam
  /// needs to undone; hence we register that beam info here:
  ///@{
  casacore::MDirection itsDirectionAtStart;
  everybeam::CorrectionMode itsModeAtStart;
  ///@}

  unsigned int itsDebugLevel;

  /// The info needed to calculate the station beams.
  ///@{
  std::unique_ptr<everybeam::telescope::Telescope> telescope_;
  casacore::MeasFrame measure_frame_;
  casacore::MDirection::Convert measure_converter_;
  std::vector<aocommon::MC2x2> beam_values_;
  std::vector<size_t> ant_to_msindex_;
  bool use_model_data_;
  ///@}

  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
