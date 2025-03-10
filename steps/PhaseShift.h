// PhaseShift.h: DPPP step class to shift the data to another phase center
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @author Ger van Diepen

#ifndef DP3_STEPS_PHASESHIFT_H_
#define DP3_STEPS_PHASESHIFT_H_

#include "InputStep.h"

#include <dp3/base/DPBuffer.h>

#include <aocommon/staticfor.h>

#include <casacore/casa/Arrays/Matrix.h>

namespace dp3 {
namespace base {
struct Direction;
}

namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DP3 step class to shift the data to another phase center

/// This class is a Step class to shift the data and UVW coordinates
/// to another phase center. If no phase center is given, a shift is
/// done back to the original phase center.
/// The code is based on the script phaseshift.py by Bas vd Tol.

class PhaseShift : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// This is the standard constructor where the phasecenter must be given.
  PhaseShift(const common::ParameterSet&, const std::string& prefix);

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// This is a constructor for Demixer where the phasecenter has the
  /// given default value.
  PhaseShift(const common::ParameterSet&, const std::string& prefix,
             const std::vector<std::string>& defVal);

  ~PhaseShift() override;

  common::Fields getRequiredFields() const override {
    return kDataField | kUvwField;
  }

  common::Fields getProvidedFields() const override {
    return kDataField | kUvwField;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Fill the Euler rotation matrix for given ra/dec.
  static void fillEulerMatrix(casacore::Matrix<double>& mat,
                              const base::Direction& direction);

  /// Get the phasors resulting from the last process step.
  /// This is used in the Demixer.
  const xt::xtensor<std::complex<double>, 2>& getPhasors() const {
    return itsPhasors;
  }

  /// Get the phase center.
  const std::vector<std::string>& getPhaseCenter() const { return itsCenter; }

 private:
  /// Interpret the phase center specification.
  /// Currently only J2000 RA and DEC can be given.
  casacore::MDirection handleCenter();

  std::string itsName;
  std::vector<std::string> itsCenter;
  std::vector<double> itsFreqC;  ///< freq/C
  casacore::Matrix<double> itsEulerMatrix;
  double itsXYZ[3];  ///< numpy.dot((w-w1).T, T)
  xt::xtensor<std::complex<double>, 2>
      itsPhasors;  ///< phase factor per chan,bl
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
