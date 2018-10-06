//# H5ParmPredict.h: DPPP step class to H5ParmPredict visibilities from a source model
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

#ifndef DPPP_H5ParmPredict_H
#define DPPP_H5ParmPredict_H

// @file
// @brief DPPP step class to H5ParmPredict visibilities from a source model

#include "DPInput.h"
#include "DPBuffer.h"
#include "Predict.h"
#include "H5Parm.h"

#include "../Common/ThreadPool.h"

#include <utility>

namespace DP3 {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to H5ParmPredict visibilities with optionally beam

    typedef std::pair<size_t, size_t> Baseline;
    typedef std::pair<ModelComponent::ConstPtr, Patch::ConstPtr> Source;

    class H5ParmPredict: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      H5ParmPredict (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~H5ParmPredict();

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
      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      DPBuffer         itsBuffer;

      std::vector<Predict::ShPtr> itsPredictSteps;
      ResultStep*      itsResultStep;

      std::string      itsH5ParmName;
      std::vector<std::string> itsDirections;

      NSTimer          itsTimer;
      ThreadPool itsThreadPool;
      std::mutex itsMeasuresMutex;
    };

  } //# end namespace
}

#endif
