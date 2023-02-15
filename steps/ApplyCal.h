// ApplyCal.h: DP3 step class to ApplyCal visibilities from a source model
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
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
  static void applyDiag(const casacore::Complex* gainA,
                        const casacore::Complex* gainB, casacore::Complex* vis,
                        float* weight, bool* flag, unsigned int bl,
                        unsigned int chan, bool updateWeights,
                        base::FlagCounter& flagCounter);

  /// Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H,
  /// where the solution is equal for both polarizations
  static void applyScalar(const casacore::Complex* gainA,
                          const casacore::Complex* gainB,
                          casacore::Complex* vis, float* weight, bool* flag,
                          unsigned int bl, unsigned int chan,
                          bool updateWeights, base::FlagCounter& flagCounter);

  /// Apply a full Jones matrix to the 2x2 visibilities matrix: A.V.B^H
  static void applyFull(const casacore::Complex* gainA,
                        const casacore::Complex* gainB, casacore::Complex* vis,
                        float* weight, bool* flag, unsigned int bl,
                        unsigned int chan, bool updateWeights,
                        base::FlagCounter& flagCounter);

  /// Do the same as the combination of BBS + python script
  /// covariance2weight.py (cookbook), except it stores weights per freq.
  /// The diagonal of covariance matrix is transferred to the weights.
  /// Note that the real covariance (mixing of noise terms after which they
  /// are not independent anymore) is not stored.
  /// The input covariance matrix C is assumed to be diagonal with elements
  /// w_i (the weights), the result the diagonal of
  /// (gainA kronecker gainB^H).C.(gainA kronecker gainB^H)^H
  static void applyWeights(const casacore::Complex* gainA,
                           const casacore::Complex* gainB, float* weight);

 private:
  bool itsIsSubstep{false};
  string itsName;

  std::vector<OneApplyCal::ShPtr> itsApplyCals;
};

}  // namespace steps
}  // namespace dp3

#endif
