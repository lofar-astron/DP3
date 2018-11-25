//# ScaleData.h: DPPP step class for freq-dependent scaling of the data
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
//# $Id: ScaleData.h 23223 2012-12-07 14:09:42Z schoenmakers $
//#
//# @author Ger van Diepen

#ifndef DPPP_SCALEDATA_H
#define DPPP_SCALEDATA_H

// @file
// @brief DPPP step class for freq-dependent scaling of the data

#include "DPInput.h"
#include "DPBuffer.h"

#include <casacore/casa/Arrays/Cube.h>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class scaling the data using a polynomial
    // in frequency (in MHz) for LBA and HBA. The coefficients can be given
    // as ParSet parameters having a default determined by Adam Deller.
    //
    // The polynomial coefficients can depend on station by giving them
    // per station name regular expression. The default coefficients are
    // used for the station not matching any regular expression.
    //
    // The data are multiplied with a factor sqrt(scale[ant1] * scale[ant2]).
    // An extra scale factor can be applied to compensate for the different
    // number of dipoles or tiles or for missing ones. By default that
    // extra scale factor is only applied to stations using the default
    // coefficients, because it is assumed that coefficients are scaled well
    // when specifying them explicitly for stations.

    class ScaleData: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      ScaleData (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~ScaleData();

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

    private:
      // Fill the scale factors for stations having different nr of tiles.
      void fillSizeScaleFactors (uint nNominal, std::vector<double>& fact);

      //# Data members.
      string             itsName;
      bool               itsScaleSizeGiven;
      bool               itsScaleSize;
      std::vector<string>     itsStationExp;  // station regex strings
      std::vector<string>     itsCoeffStr;    // coeff per station regex
      std::vector<std::vector<double> > itsStationFactors; // scale factor per station,freq
      casacore::Cube<double> itsFactors;     // scale factor per baseline,freq,pol
      NSTimer            itsTimer;
    };

  } //# end namespace
}

#endif
