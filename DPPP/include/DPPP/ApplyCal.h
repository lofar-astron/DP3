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
//# @author Ger van Diepen

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

namespace LOFAR {

  class ParameterSet;

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
      void applyPhase (casa::Complex* vis, int antA, int antB, int chan,
          int time);
      void applyClock (casa::Complex* vis, int antA, int antB, int chan,
          int time);
      void applyGain (casa::Complex* vis, float* weight, int antA, int antB,
          int chan, int time);

      void updateParms (const double bufStartTime);

      void initDataArrays();

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      string           itsParmDBName;
      boost::shared_ptr<BBS::ParmFacade> itsParmDB;
      string           itsCorrectType;
      uint             itsTimeSlotsPerParmUpdate;
      // Expressions to search for in itsParmDB

      vector<casa::String>   itsParmExprs;

      vector<vector<casa::DComplex> > itsParms0;
      vector<vector<casa::DComplex> > itsParms1;
      uint             itsTimeStep;
      uint             itsNCorr;
      double          itsTimeInterval;
      double          itsLastTime;
      bool            itsUseAP;      //# use ampl/phase or real/imag
      NSTimer          itsTimer;
    };

  } //# end namespace
}

#endif
