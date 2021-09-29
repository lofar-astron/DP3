// PhaseShift.h: DPPP step class to shift the data to another phase center
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @author Ger van Diepen

#ifndef DP3_PHASESHIFT_H
#define DP3_PHASESHIFT_H

#include "InputStep.h"

#include "../base/DPBuffer.h"

#include <casacore/casa/Arrays/Matrix.h>

namespace dp3 {
namespace base {
struct Direction;
}

namespace common {
class ParameterSet;
}

namespace steps {

/// @brief DPPP step class to shift the data to another phase center

/// This class is a Step class to shift the data and UVW coordinates
/// to another phase center. If no phase center is given, a shift is
/// done back to the original phase center.
/// The code is based on the script phaseshift.py by Bas vd Tol.

class PhaseShift : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// This is the standard constructor where the phasecenter must be given.
  PhaseShift(InputStep*, const common::ParameterSet&, const string& prefix);

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// This is a constructor for Demixer where the phasecenter has the
  /// given default value.
  PhaseShift(InputStep*, const common::ParameterSet&, const string& prefix,
             const std::vector<string>& defVal);

  virtual ~PhaseShift();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Fill the transformation matrix for given ra/dec.
  static void fillTransMatrix(casacore::Matrix<double>& mat,
                              const base::Direction& direction);

  /// Get the phasors resulting from the last process step.
  /// This is used in the Demixer.
  const casacore::Matrix<casacore::DComplex>& getPhasors() const {
    return itsPhasors;
  }

  /// Get the phase center.
  const std::vector<string>& getPhaseCenter() const { return itsCenter; }

 private:
  /// Interpret the phase center specification.
  /// Currently only J2000 RA and DEC can be given.
  casacore::MDirection handleCenter();

  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBuf;
  std::vector<string> itsCenter;
  std::vector<double> itsFreqC;      ///< freq/C
  casacore::Matrix<double> itsMat1;  ///< TT in phasehift.py
  double itsXYZ[3];                  ///< numpy.dot((w-w1).T, T)
  casacore::Matrix<casacore::DComplex>
      itsPhasors;  ///< phase factor per chan,bl
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
