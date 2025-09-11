// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to create dynamic spectra of particular sources
/// @author Mick Veldhuis

#ifndef DP3_STEPS_DYNSPEC_H_
#define DP3_STEPS_DYNSPEC_H_

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/DynSpecFitsWriter.h"

#include "../model/Patch.h"

#include "ResultStep.h"

namespace dp3::steps {

/// @brief DP3 step class that creates visibility averaged dynamic spectra.
/// TODO: (1) Output the weights for each pixel in the spectrum. (2) Select
/// polarizations to output, e.g., only I or instrumental polarizations.
class DynSpec : public Step {
 public:
  /// Type for the dynamic spectra, with axes time x frequency x (4) Stokes
  /// parameters x direction. The FITS format enforces a column-major layout,
  /// hence, we adopt it here as well.
  using DynamicSpectrumTensor =
      xt::xtensor<float, 4, xt::layout_type::column_major>;

  DynSpec(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const final {
    return kDataField | kUvwField | kWeightsField | kFlagsField;
  }

  common::Fields getProvidedFields() const final { return {}; }

  bool process(std::unique_ptr<base::DPBuffer> buffer) final;

  void finish() final;

  void updateInfo(const base::DPInfo&) final;

  void show(std::ostream&) const final;

  void showTimings(std::ostream&, double duration) const final;

 private:
  /// Convert instrumental polarizations to Stokes parameters.
  xt::xtensor<float, 2> ComputeAbsoluteStokesParameters(
      xt::xtensor<std::complex<float>, 2>& baseline_averaged_data) const;

  /// Store spectra in a FITS compatible format.
  void WriteSpectraToDisk();

  size_t time_index_{0};

  /// FITS filename prefix
  std::string fits_prefix_;

  /// File with Name, RA, and Dec of the sources to create spectra for.
  std::string source_file_name_;

  /// Patches read from source_list_.
  std::vector<std::shared_ptr<model::Patch>> source_list_;

  /// Whether to run H5ParmPredict to subtract sources before measuring the
  /// spectrum.
  bool subtract_sources_{false};
  bool subtract_with_h5parmpredict_{false};
  bool subtract_model_column_{false};

  /// Whether calibration solutions are corrected for.
  /// This adds an ApplyCal substep for each source direction.
  bool apply_calibration_solutions_{false};

  /// Apply beam corrections in the direction of each source are applied.
  /// This adds an ApplyBeam substep for each source direction.
  bool apply_beam_correction_{false};
  bool apply_beam_reweighted_{false};

  /// The name of the model column with the foreground, if provided.
  std::string model_column_;

  /// Holds the visibility-averaged spectra for each direction of interest in
  /// memory until written to disk.
  DynamicSpectrumTensor dynamic_spectra_;

  /// Internal steps to predict the sources in the field, shift to the target
  /// direction and capture the results.
  /// {
  std::vector<std::shared_ptr<Step>> first_substeps_;
  std::vector<std::shared_ptr<ResultStep>> results_;
  std::shared_ptr<Step> model_step_;
  std::shared_ptr<ResultStep> model_result_;
  /// }

  std::string name_;
  common::NSTimer total_timer_;
  common::NSTimer substep_timer_;
  common::NSTimer computation_timer_;
  common::NSTimer write_timer_;
};

}  // namespace dp3::steps

#endif
