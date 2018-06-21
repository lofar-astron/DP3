//# DDE.h: DPPP step class to calibrate direction dependent gains
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
//# $Id: DDECal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_DDECAL_H
#define DPPP_DDECAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include <DPPP/DPInput.h>
#include <DPPP/GainCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/H5Parm.h>
#include <DPPP/BaselineSelection.h>
#include <DPPP/Patch.h>
#include <DPPP/UVWFlagger.h>
#include <DPPP/Predict.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/ApplyBeam.h>
#include <DPPP_DDECal/MultiDirSolver.h>
#include <DPPP_DDECal/Constraint.h>
#include <StationResponse/Station.h>
#include <StationResponse/Types.h>
#include <ParmDB/Parm.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Arrays/ArrayMath.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to calibrate (direction independent) gains.

    typedef vector<Patch::ConstPtr> PatchList;
    typedef std::pair<size_t, size_t> Baseline;

    class DDECal: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      DDECal (DPInput*, const ParameterSet&, const std::string& prefix);

      virtual ~DDECal();

      // Create an DDECal object using the given parset.
      static DPStep::ShPtr makeStep (DPInput*, const ParameterSet&,
                                     const std::string&);

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Call the actual solver (called once per solution interval)
      void doSolve();

      // Initialize H5parm-file
      void initH5parm();

      // Write out the solutions
      void writeSolutions();

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;


    private:
      // Initialize solutions
      void initializeScalarSolutions();
      
      void initializeFullMatrixSolutions();

      // Convert itsDirections to a vector of strings like "[Patch1, Patch2]"
      // Used for setting source names.
      std::vector<std::string> getDirectionNames();

      //# Data members.
      DPInput*         itsInput;
      std::string      itsName;
      vector<DPBuffer> itsBufs;

      // The time of the current buffer (in case of solint, average time)
      double           itsAvgTime;
      std::vector<casacore::Complex*> itsDataPtrs;

      // For each timeslot, a vector of nDir buffers
      std::vector<std::vector<casacore::Complex*> > itsModelDataPtrs;

      // For each time, for each channel block, a vector of size nAntennas * nDirections
      std::vector<std::vector<std::vector<casacore::DComplex> > > itsSols;
      std::vector<uint>
        itsNIter, // Number of iterations taken
        itsNApproxIter;

      // For each time, for each constraint, a vector of results (e.g. tec and phase)
      std::vector<std::vector<std::vector<Constraint::Result> > > itsConstraintSols;

      casacore::Cube<casacore::Complex> itsModelData;
      std::string      itsH5ParmName;
      H5Parm           itsH5Parm;
      std::string      itsParsetString; // Parset, for logging in H5Parm

      GainCal::CalType itsMode;
      bool             itsPropagateSolutions;
      uint             itsTimeStep;
      uint             itsSolInt;
      uint             itsStepInSolInt;
      uint             itsNChan;
      vector<size_t>   itsChanBlockStart;    // For each channel block, the index in the channels at which this channel block starts
      vector<double>   itsChanBlockFreqs;
      vector<vector<string> > itsDirections; // For each direction, a vector of patches
      vector<casacore::CountedPtr<Constraint> > itsConstraints;

      vector<double>   itsWeights;

      UVWFlagger       itsUVWFlagStep;
      ResultStep::ShPtr itsDataResultStep; // Result step for data after UV-flagging
      vector<Predict>     itsPredictSteps;
      vector<MultiResultStep::ShPtr> itsResultSteps; // For each directions, a multiresultstep with all times

      NSTimer          itsTimer;
      NSTimer          itsTimerPredict;
      NSTimer          itsTimerSolve;
      NSTimer          itsTimerWrite;
      double           itsCoreConstraint;
      double           itsScreenCoreConstraint;

      MultiDirSolver   itsMultiDirSolver;
      bool itsFullMatrixMinimalization;
      bool itsApproximateTEC;
      std::string itsStatFilename;
      std::unique_ptr<std::ofstream> itsStatStream;
    };

  } //# end namespace
}

#endif
