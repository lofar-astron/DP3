//# ApplyCal.h: DPPP step class to ApplyCal visibilities from a source model
//# Copyright (C) 2013
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id:
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_ApplyCal_H
#define DPPP_ApplyCal_H

// @file
// @brief DPPP step class to apply multiple calibration solutions

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>

#include <DPPP/OneApplyCal.h>

#include <utility>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to apply multiple ParmDB or H5Parm
    // solutions to data.

    class ApplyCal: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      ApplyCal (DPInput*, const ParameterSet&, const string& prefix,
                bool substep=false, std::string predictDirection="");

      // Empty constructor
      ApplyCal ();

      virtual ~ApplyCal();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Set the next step. It squeezes in the actual OneApplyCal steps
      // between this ApplyCal step and the next step.
      virtual void setNextStep (DPStep::ShPtr nextStep);

      // Show the step. When ApplyCal is a step in the main chain, this does
      // nothing; the nextStep mechanism in DPRun will call show on the actual
      // OneApplyCals.
      virtual void show(std::ostream&) const;

      // Show the timings. When ApplyCal is a step in the main chain, this does
      // nothing; the nextStep mechanism in DPRun will call show on the actual
      // OneApplyCals.
      virtual void showTimings (std::ostream&, double duration) const;

      // Invert a 2x2 matrix in place
      static void invert (casacore::DComplex* v, double sigmaMMSE=0);

      // Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      static void applyDiag (const casacore::DComplex* gainA,
                             const casacore::DComplex* gainB,
                             casacore::Complex* vis, float* weight, bool* flag,
                             uint bl, uint chan, bool updateWeights,
                             FlagCounter& flagCounter);

      // Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H,
      // where the solution is equal for both polarizations
      static void applyScalar(const casacore::DComplex* gainA,
                              const casacore::DComplex* gainB,
                              casacore::Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool updateWeights,
                              FlagCounter& flagCounter);

      // Apply a full Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      static void applyFull (const casacore::DComplex* gainA,
                             const casacore::DComplex* gainB,
                             casacore::Complex* vis, float* weight, bool* flag,
                             uint bl, uint chan, bool updateWeights,
                             FlagCounter& flagCounter);

      // Do the same as the combination of BBS + python script
      // covariance2weight.py (cookbook), except it stores weights per freq.
      // The diagonal of covariance matrix is transferred to the weights.
      // Note that the real covariance (mixing of noise terms after which they
      // are not independent anymore) is not stored.
      // The input covariance matrix C is assumed to be diagonal with elements
      // w_i (the weights), the result the diagonal of
      // (gainA kronecker gainB^H).C.(gainA kronecker gainB^H)^H
      static void applyWeights (const casacore::DComplex* gainA,
                                const casacore::DComplex* gainB,
                                float* weight);

    private:
      //# Data members.
      bool             itsIsSubstep;
      string           itsName;

      std::vector<OneApplyCal::ShPtr> itsApplyCals;
    };

  } //# end namespace
}

#endif
