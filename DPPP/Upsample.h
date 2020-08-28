// Upsample.h: DPPP step class to upsample visibilities
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
/// @brief DPPP step class to Upsample visibilities
/// @author Tammo Jan Dijkema

#ifndef DPPP_Upsample_H
#define DPPP_Upsample_H

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    /// @brief DPPP step class to Upsample visibilities

    /// This class is an empty DPStep subclass to use as implementation template

    class Upsample: public DPStep
    {
    public:
      /// Construct the object.
      /// Parameters are obtained from the parset using the given prefix.
      Upsample (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~Upsample();

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

    private:
      string                itsName;
      double                itsOldTimeInterval;
      unsigned int                  itsTimeStep;

      std::vector<DPBuffer> itsPrevBuffers;
      std::vector<DPBuffer> itsBuffers;
      unsigned int                  itsFirstToFlush;

      NSTimer               itsTimer;
    };

  } // end namespace
}

#endif
