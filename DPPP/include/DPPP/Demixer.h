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
#include <DPPP/BBSExpr.h>

#include <casa/Arrays/Cube.h>
#include <measures/Measures/MDirection.h>

#include <BBSKernel/MeasurementExprLOFAR.h>
#include <BBSKernel/BaselineMask.h>
#include <BBSKernel/CorrelationMask.h>
#include <BBSKernel/ParmManager.h>
#include <BBSKernel/VisBuffer.h>
#include <ParmDB/SourceDB.h>
#include <DPPP/EstimateNDPPP.h>

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
      casa::MDirection handleCenter(const vector<string> &center) const;

      // Solve gains and subtract sources.
      void demix();

//      // Do the subtraction.
//      void subtract();

      // Add the decorrelation factor contribution for each time slot.
      void addFactors (const DPBuffer& newBuf);

      // Calculate the decorrelation factors by averaging them.
      void averageFactors();

      //# Data members.
      DPInput*                 itsInput;
      string                   itsName;
      vector<PhaseShift*>      itsPhaseShifts;
      vector<DPStep::ShPtr>    itsFirstSteps;   //# phaseshift/average steps
      vector<MultiResultStep*> itsAvgResults;
      vector<BBSExpr::ShPtr>   itsBBSExpr;
      vector<BBS::MeasurementExpr::Ptr> itsModels;
      string                   itsTarget;
      vector<string>           itsSources;
      vector<string>           itsExtraSources;
      vector<string>           itsAllSources;
//      vector<DPBuffer>         itsBuf;
      double                   itsTimeStart;
      double                   itsTimeInterval;
      vector<double>           itsTimeCenters;
      vector<double>           itsTimeWidths;
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
      double                   itsTimeIntervalAvg;
      double                   itsTimeIntervalRes;
      casa::Array<casa::DComplex> itsFactorBuf; //# ncorr,nchan,nbl,ndir*ndir
      vector<casa::Array<casa::DComplex> > itsFactors; //# demix factors/time
      //# each Array is basically cube(ncorr,nchan,nbl) of matrix(ndir,ndir)
      NSTimer                  itsTimer;
      NSTimer                  itsTimerPhaseShift;
      NSTimer                  itsTimerDemix;
      NSTimer                  itsTimerBBS;
      NSTimer                  itsTimerSubtract;

      boost::shared_ptr<BBS::SourceDB> itsSourceDB;
      vector<BBS::ParmGroup>           itsModelParms;
      BBS::ParmGroup                   itsParms;
      BBS::BaselineSeq                 itsBaselines;
      BBS::BaselineMask                itsBaselineMask;
      BBS::CorrelationSeq              itsCorrelations;
      BBS::CorrelationMask             itsCorrelationMask;
      BBS::Axis::ShPtr                 itsFreqAxisAvg;
      BBS::EstimateOptions             itsOptions;
    };

  } //# end namespace
}

#endif
