//# PhaseShift.h: DPPP step class to average in time and/or freq
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
// @brief DPPP step class to average in time and/or freq

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <measures/Measures/UVWMachine.h>

namespace LOFAR {

  namespace DPPP {
    class ParSet;

    // @ingroup NDPPP

    // This class is a DPStep class calculating the weighted average of
    // data in time and/or frequency.
    // <br>
    // Only unflagged data points are used. The average is calculated as
    // <tt>sum(data*weight) / sum(weight)</tt> and the sum of the weights
    // is the weight of the new data point. If all data point to use are
    // flagged, the resulting data point and weight are set to zero and flagged.
    //
    // It keeps track of the FullResFlags. It sets them if the corresponding
    // data point is flagged. Note that multiple FullResFlags elements map to
    // a single data point if some averaging was done before.

    class PhaseShift: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      PhaseShift (DPInput*, const ParSet&, const string& prefix);

      virtual ~PhaseShift();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Interpret the phase center specification.
      // Currently only J2000 RA and DEC can be given.
      casa::MDirection handleCenter();
      
      //# Data members.
      DPInput*          itsInput;
      string            itsName;
      vector<string>    itsCenter;
      vector<double>    itsFreqC;      //# freq/C
      casa::UVWMachine* itsMachine;
      NSTimer           itsTimer;
    };

  } //# end namespace
}

#endif
