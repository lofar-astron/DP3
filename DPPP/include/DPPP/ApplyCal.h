//# ApplyCal.h: DPPP step class to apply a calibration correction to the data
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
//# $Id: ApplyCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_APPLYCAL_H
#define DPPP_APPLYCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <ParmDB/ParmFacade.h>
#include <ParmDB/ParmSet.h>
#include <ParmDB/Parm.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <DPPP/FlagCounter.h>

namespace LOFAR {
  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class applying calibration parameters to the data.

    class ApplyCal: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      ApplyCal (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~ApplyCal();

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
      // Apply a diagonal Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      void applyDiag (casa::Complex* vis, float* weight, int antA, int antB,
          int chan, int time);

      // Apply a full Jones matrix to the 2x2 visibilities matrix: A.V.B^H
      void applyFull (casa::Complex* vis, float* weight, int antA, int antB,
          int chan, int time);

      // Read parameters from the associated parmdb and store them in itsParms
      void updateParms (const double bufStartTime);

      // Invert a 2x2 matrix in place
      void invert (casa::DComplex* v, double sigmaMMSE=0) const;

      void initDataArrays();

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      string           itsParmDBName;
      boost::shared_ptr<BBS::ParmFacade> itsParmDB;
      string           itsCorrectType;
      uint             itsTimeSlotsPerParmUpdate;
      double           itsSigmaMMSE;
      bool             itsUpdateWeights;

      // Expressions to search for in itsParmDB
      vector<casa::String>   itsParmExprs;

      // parameters, numparms, antennas, time x frequency
      vector<vector<vector<casa::DComplex> > > itsParms;
      uint            itsTimeStep;
      uint            itsNCorr;
      double          itsTimeInterval;
      double          itsLastTime;
      FlagCounter     itsFlagCounter;
      bool            itsUseAP;      //# use ampl/phase or real/imag
      NSTimer         itsTimer;
    };

  } //# end namespace
}

#endif
