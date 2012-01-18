//# Demixer.h: DPPP step class to subtract A-team sources
//# Copyright (C) 2011
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

#ifndef DPPP_DEMIXER_H
#define DPPP_DEMIXER_H

// @file
// @brief DPPP step class to average in time and/or freq

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/PhaseShift.h>

#include <casa/Arrays/Cube.h>
#include <measures/Measures/MDirection.h>

namespace LOFAR {

  namespace DPPP {
    class ParSet;

    // @ingroup NDPPP

    // This class is a DPStep class to subtract the strong A-team sources.
    // It is based on the demixing.py script made by Bas vd Tol and operates
    // per time chunk as follows:
    // <ul>
    //  <li> The data are phase-shifted and averaged for each source.
    //  <li> Demixing is done using the combined results.
    //  <li> For each source a BBS solve, smooth, and predict is done.
    //  <li> The predicted results are subtracted from the averaged data.
    // </ul>

    class Demixer: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Demixer (DPInput*, const ParSet&, const string& prefix);

      virtual ~Demixer();

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
      // Demix and return the result.
      DPBuffer demix() const;

      // Subtract sources and return the result.
      DPBuffer subtract() const;

      // Fill the direction rotation matrices.
      void fillMatrices();

      void addFactors (const DPBuffer& newBuf);
      void averageFactors();

      void calcDirs();

      void demix();

      // Convert dir index to a linear index.
      int toIndex (int ndir, int i, int j) const
        { return (ndir + ndir + 1 - i)/2 + j - i; }

      //# Data members.
      DPInput*                 itsInput;
      string                   itsName;
      vector<PhaseShift*>      itsPhaseShifts;
      vector<DPStep::ShPtr>    itsFirstSteps;   //# phaseshift/average chain
      vector<ResultStep*>      itsDemixResults;
      vector<DPStep::ShPtr>    itsSecondSteps;
      vector<ResultStep*>      itsSubtractInputs;
      vector<string>           itsSources;
      ///      vector<string>        itsFieldSources;
      ///      vector<string>        itsMovingSources;
      vector<DPBuffer>         itsBuf;
      bool                     itsJointSolve;
      uint                     itsNrDir;
      uint                     itsNrBl;
      uint                     itsNrCorr;
      uint                     itsNrChanIn;
      uint                     itsNrChanOut;
      uint                     itsNChanAvg;
      uint                     itsNTimeAvg;
      uint                     itsResChanAvg;
      uint                     itsResTimeAvg;
      uint                     itsNTimeChunk;
      uint                     itsNTimeIn;
      uint                     itsNTimeOut;
      casa::Array<casa::DComplex> itsFactorBuf; //# ncorr,nchan,nbl,ndir*ndir
      vector<casa::Array<casa::DComplex> > itsFactors; //# demix factors/time
      //# each Array is basically cube(ncorr,nchan,nbl) of matrix(ndir,ndir)
      NSTimer                  itsTimer;
      NSTimer                  itsTimerPhaseShift;
      NSTimer                  itsTimerDemix;
      NSTimer                  itsTimerBBS;
      NSTimer                  itsTimerSubtract;
    };

  } //# end namespace
}

#endif
