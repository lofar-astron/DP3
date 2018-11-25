//# OneApplyCal.h: DPPP step class to apply a calibration correction to the data
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
//# $Id: OneApplyCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_ONEAPPLYCAL_H
#define DPPP_ONEAPPLYCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include "DPInput.h"
#include "DPBuffer.h"
#include "H5Parm.h"
#include "FlagCounter.h"

#include "../ParmDB/ParmFacade.h"
#include "../ParmDB/ParmSet.h"
#include "../ParmDB/Parm.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <mutex>

namespace DP3 {
  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class applying calibration parameters to the data.
    // It only applies one correction.

    class OneApplyCal: public DPStep
    {
    public:
      // Define the shared pointer for this type.
      typedef std::shared_ptr<OneApplyCal> ShPtr;

      enum class InterpolationType {NEAREST, LINEAR};

      enum CorrectType {GAIN, FULLJONES, TEC, CLOCK, ROTATIONANGLE, SCALARPHASE, PHASE,
                        ROTATIONMEASURE, SCALARAMPLITUDE, AMPLITUDE};

      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      OneApplyCal (DPInput*, const ParameterSet&, const std::string& prefix,
                const std::string& defaultPrefix,
                bool substep=false, std::string predictDirection="");

      virtual ~OneApplyCal();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer& buffer)
      {
        return process(buffer, nullptr);
      }

      bool process (const DPBuffer&, std::mutex* hdf5Mutex);
      
      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      bool invert() {
        return itsInvert;
      }

    private:
      // Read parameters from the associated parmdb and store them in itsParms
      void updateParms (const double bufStartTime, std::mutex* hdf5Mutex);

      // If needed, show the flag counts.
      virtual void showCounts (std::ostream&) const;

      void initDataArrays();

      // Check the number of polarizations in the parmdb or h5parm
      uint nPol(const std::string& parmName);

      // Replace values by NaN on places where weight is zero
      static void applyFlags(std::vector<double>& values,
                             const std::vector<double>& weights);

      static std::string correctTypeToString(CorrectType);
      static CorrectType stringToCorrectType(const string&);

      //# Data members.
      DPInput*         itsInput;
      DPBuffer         itsBuffer;
      string           itsName;
      string           itsParmDBName;
      bool             itsUseH5Parm;
      string           itsSolSetName;
      std::shared_ptr<BBS::ParmFacade> itsParmDB;
      H5Parm           itsH5Parm;
      string           itsSolTabName;
      H5Parm::SolTab   itsSolTab;
      H5Parm::SolTab   itsSolTab2; // in the case of full Jones, amp and phase table need to be open
      CorrectType      itsCorrectType;
      bool             itsInvert;
      InterpolationType itsInterpolationType;
      uint             itsTimeSlotsPerParmUpdate;
      double           itsSigmaMMSE;
      bool             itsUpdateWeights;

      uint             itsCount; // number of steps

      // Expressions to search for in itsParmDB
      std::vector<casacore::String>   itsParmExprs;

      // parameters, numparms, antennas, time x frequency
      casacore::Cube<casacore::DComplex> itsParms;
      uint            itsTimeStep; // time step within current chunk
      uint            itsNCorr;
      double          itsTimeInterval;
      double          itsLastTime; // last time of current chunk
      FlagCounter     itsFlagCounter;
      bool            itsUseAP;      //# use ampl/phase or real/imag
      hsize_t         itsDirection;
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
