//# PhaseShift.h: DPPP step class to shift the data to another phase center
//# Copyright (C) 2010
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
//# $Id$
//#
//# @author Ger van Diepen

#ifndef DPPP_PHASESHIFT_H
#define DPPP_PHASESHIFT_H

// @file
// @brief DPPP step class to shift the data to another phase center

#include "DPInput.h"
#include "DPBuffer.h"

#include <casacore/casa/Arrays/Matrix.h>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to shift the data and UVW coordinates
    // to another phase center. If no phase center is given, a shift is
    // done back to the original phase center.
    //
    // The code is based on the script phaseshift.py by Bas vd Tol.

    class PhaseShift: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      // This is the standard constructor where the phasecenter must be given.
      PhaseShift (DPInput*, const ParameterSet&, const string& prefix);

      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      // This is a constructor for Demixer where the phasecenter has the
      // given default value.
      PhaseShift (DPInput*, const ParameterSet&, const string& prefix,
                  const std::vector<string>& defVal);

      virtual ~PhaseShift();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      // Fill the transformation matrix for given ra/dec.
      static void fillTransMatrix (casacore::Matrix<double>& mat,
                                   double ra, double dec);

      // Get the phasors resulting from the last process step.
      // This is used in the Demixer.
      const casacore::Matrix<casacore::DComplex>& getPhasors() const
        { return itsPhasors; }

      // Get the phase center.
      const std::vector<string>& getPhaseCenter() const
        { return itsCenter; }

    private:
      // Interpret the phase center specification.
      // Currently only J2000 RA and DEC can be given.
      casacore::MDirection handleCenter();
      
      //# Data members.
      DPInput*             itsInput;
      string               itsName;
      DPBuffer             itsBuf;
      std::vector<string>       itsCenter;
      std::vector<double>       itsFreqC;      //# freq/C
      casacore::Matrix<double> itsMat1;       //# TT in phasehift.py
      double               itsXYZ[3];     //# numpy.dot((w-w1).T, T)
      casacore::Matrix<casacore::DComplex> itsPhasors; //# phase factor per chan,bl
      NSTimer              itsTimer;
    };

  } //# end namespace
}

#endif
