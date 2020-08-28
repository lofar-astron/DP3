// AddNoiseLBA.h: DPPP step class to add LBA random noise to data
// Copyright (C) 2013
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief DPPP step class to add LBA random noise to data
/// @author Tammo Jan Dijkema
//

#ifndef DPPP_AddNoiseLBA_H
#define DPPP_AddNoiseLBA_H

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

#define POL_DEGREE 5

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    /// @brief DPPP step class to AddNoiseLBA visibilities from a source model

    /// This class is an empty DPStep subclass to use as implementation template

    class AddNoiseLBA: public DPStep
    {
    public:
      /// Construct the object.
      /// Parameters are obtained from the parset using the given prefix.
      AddNoiseLBA (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~AddNoiseLBA();

      /// Process the data.
      /// It keeps the data.
      /// When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      /// Finish the processing of this step and subsequent steps.
      virtual void finish();

      /// Update the general info.
      virtual void updateInfo (const DPInfo&);

      /// Show the step parameters.
      virtual void show (std::ostream&) const;

      /// Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      DPInput*         itsInput;
      string           itsName;
      DPBuffer         itsBuffer;

      NSTimer          itsTimer;
      // SET mode = 0 data are modified with the noise > data = data + noise 
      // ADD mode = 1 a new array is created > newdata = data + noise
      int              mode;
      // coefficients for polynomial interpolation (from constant -first- to highet order -last-)
      double           coeffs_outer[POL_DEGREE+1] = {4.46492043e+05, -4.04156579e-02,  1.58636639e-09, -3.09364148e-17,  2.93955326e-25, -1.06998148e-33};
      double           coeffs_inner[POL_DEGREE+1] = {8.32889327e+05, -8.93829326e-02,  3.90153820e-09, -8.23245656e-17,  8.35181243e-25, -3.25202160e-33};
      // lba_mode can be: lba_inner or lba_outer. Read from the parset file
      string           lba_mode;
      // system efficiency: for the moment set to 1.0
      double           eta = 1.0;
      // exposure: read from the measurement set: for the moment set to 1.0
      float            exposure = 1.0;
      // frequencies: 
       
      // channel width: read from the masurement set: for the moment set to 1.0



    };

  } // end namespace
}

#endif
