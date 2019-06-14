//# Split.h: DPPP step class to Split visibilities from a source model
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

#ifndef DPPP_Split_H
#define DPPP_Split_H

// @file
// @brief DPPP step class to Split visibilities from a source model

#include "DPInput.h"
#include "DPBuffer.h"

#include <utility>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is an empty DPStep subclass to use as implementation template

    class Split: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Split (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~Split();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      virtual void addToMS (const string&);

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      //# Data members.
      string           itsName;

      std::vector<std::string>   itsReplaceParms; // The names of the parameters that differ along the substeps
      std::vector<DPStep::ShPtr> itsSubsteps;
      bool             itsAddedToMS; // Used in addToMS to prevent recursion
    };

  } //# end namespace
}

#endif
