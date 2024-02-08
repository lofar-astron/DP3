// ApplyCal.h: DP3 step class to ApplyCal visibilities from a source model
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply multiple calibration solutions
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_APPLYCAL_H_
#define DP3_STEPS_APPLYCAL_H_

#include <utility>

#include <dp3/base/DPBuffer.h>
#include "OneApplyCal.h"

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

/// \brief DPPP step class to ApplyCal visibilities from a source model

/// This class is a Step class to apply multiple ParmDB or H5Parm
/// solutions to data.

class ApplyCal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ApplyCal(const common::ParameterSet&, const string& prefix,
           bool substep = false, std::string predictDirection = "");

  ApplyCal() = default;

  ~ApplyCal() override = default;

  common::Fields getRequiredFields() const override {
    // ApplyCal is a dummy step, which is followed by OneApplyCal steps.
    return {};
  }

  common::Fields getProvidedFields() const override {
    // ApplyCal is a dummy step, which is followed by OneApplyCal steps.
    return {};
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Set the next step. It squeezes in the actual OneApplyCal steps
  /// between this ApplyCal step and the next step.
  void setNextStep(Step::ShPtr nextStep) override;

  /// Show the step. When ApplyCal is a step in the main chain, this does
  /// nothing; the nextStep mechanism in DPRun will call show on the actual
  /// OneApplyCals.
  void show(std::ostream&) const override;

  /// Show the timings. When ApplyCal is a step in the main chain, this does
  /// nothing; the nextStep mechanism in DPRun will call show on the actual
  /// OneApplyCals.
  void showTimings(std::ostream&, double duration) const override;

  /// Invert a 2x2 matrix in place
  template <typename NumType>
  static void invert(std::complex<NumType>* v, NumType sigmaMMSE = 0);

  /// Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H
  static void ApplyDiag(const std::complex<float>* gain_a,
                        const std::complex<float>* gain_b,
                        base::DPBuffer& buffer, unsigned int baseline,
                        unsigned int channel, bool update_weights,
                        base::FlagCounter& flag_counter);

  /// Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H,
  /// where the solution is equal for both polarizations
  static void ApplyScalar(const std::complex<float>* gain_a,
                          const std::complex<float>* gain_b,
                          base::DPBuffer& buffer, unsigned int baseline,
                          unsigned int channel, bool update_weights,
                          base::FlagCounter& flag_counter);

  /// Apply a full Jones matrix to a 2x2 visibilities matrix: A.V.B^H
  /// @param buffer ApplyFull updates the visibilities, flags and weights in
  /// this buffer.
  /// @param baseline The baseline index for the visibility.
  /// @param channel The channel index for the visibility.
  static void ApplyFull(const std::complex<float>* gain_a,
                        const std::complex<float>* gain_b,
                        base::DPBuffer& buffer, unsigned int baseline,
                        unsigned int channel, bool update_weights,
                        base::FlagCounter& flag_counter);

  /// Do the same as the combination of BBS + python script
  /// covariance2weight.py (cookbook), except it stores weights per freq.
  /// The diagonal of covariance matrix is transferred to the weights.
  /// Note that the real covariance (mixing of noise terms after which they
  /// are not independent anymore) is not stored.
  /// The input covariance matrix C is assumed to be diagonal with elements
  /// w_i (the weights), the result the diagonal of
  /// (gainA kronecker gainB^H).C.(gainA kronecker gainB^H)^H
  static void ApplyWeights(const std::complex<float>* gain_a,
                           const std::complex<float>* gain_b, float* weight);

 private:
  bool is_sub_step_{false};
  std::vector<std::shared_ptr<OneApplyCal>> apply_cals_;
};

}  // namespace steps
}  // namespace dp3

#endif
