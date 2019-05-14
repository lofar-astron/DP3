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

#include "../DPPP/DPInput.h"
#include "../DPPP/GainCal.h"
#include "../DPPP/DPBuffer.h"
#include "../DPPP/H5Parm.h"
#include "../DPPP/BaselineSelection.h"
#include "../DPPP/Patch.h"
#include "../DPPP/UVWFlagger.h"
#include "../DPPP/Predict.h"
#include "../DPPP/SourceDBUtil.h"
#include "../DPPP/ApplyBeam.h"

#include "MultiDirSolver.h"
#include "Constraint.h"

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/Station.h>
#include <StationResponse/Types.h>
#endif

#include "../ParmDB/Parm.h"

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <string>
#include <vector>

class FacetPredict;

namespace DP3 {

  class ParameterSet;
	class ThreadPool;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to calibrate (direction independent) gains.

    typedef std::vector<Patch::ConstPtr> PatchList;
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

      void checkMinimumVisibilities();

      void flagChannelBlock(size_t cbIndex);

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
      void initializeConstraints(const ParameterSet& parset, const string& prefix);
      void initializeIDG(const ParameterSet& parset, const string& prefix);
      void initializePredictSteps(const ParameterSet& parset, const string& prefix);

      // Initialize solutions
      void initializeScalarSolutions();

      void initializeFullMatrixSolutions();

      // Convert itsDirections to a vector of strings like "[Patch1, Patch2]"
      // Used for setting source names.
      std::vector<std::string> getDirectionNames();

      void subtractCorrectedModel(bool fullJones);
      
      void idgCallback(size_t row, size_t direction, size_t dataDescId, const std::complex<float>* values);

      //# Data members.
      DPInput*         itsInput;
      std::string      itsName;
      std::vector<DPBuffer> itsBufs;
      std::vector<casacore::Cube<bool>> itsOriginalFlags;
      std::vector<casacore::Cube<float>> itsOriginalWeights;

      bool itsUseModelColumn;
      std::vector<casacore::Cube<casacore::Complex>> itsModelData;

      // The time of the current buffer (in case of solint, average time)
      double           itsAvgTime;
      std::vector<casacore::Complex*> itsDataPtrs;
      std::vector<float*> itsWeightPtrs;

      // For each timeslot, a vector of nDir buffers, each of size nbl x nch x npol
      std::vector<std::vector<casacore::Complex*> > itsModelDataPtrs;
      
      std::vector<std::vector<std::vector<casacore::Complex>>> itsIDGBuffers;

      // For each time, for each channel block, a vector of size nAntennas * nDirections
      std::vector<std::vector<std::vector<casacore::DComplex> > > itsSols;
      std::vector<size_t>
        itsNIter, // Number of iterations taken
        itsNApproxIter;

      // For each time, for each constraint, a vector of results (e.g. tec and phase)
      std::vector<std::vector<std::vector<Constraint::Result> > > itsConstraintSols;

      std::string      itsH5ParmName;
      H5Parm           itsH5Parm;
      std::string      itsParsetString; // Parset, for logging in H5Parm

      GainCal::CalType itsMode;
      bool itsPropagateSolutions;
      bool itsPropagateConvergedOnly;
      bool itsFlagUnconverged;
      bool itsFlagDivergedOnly;
      bool itsUseIDG;
      bool itsOnlyPredict;
      size_t itsTimeStep;
      size_t itsSolInt;
      double itsMinVisRatio;
      size_t itsStepInSolInt;
      size_t itsNChan;
      // For each channel block, the nr of unflagged vis and the total nr of vis.
      std::vector<std::pair<size_t, size_t>> itsVisInInterval;
      std::vector<size_t> itsChanBlockStart;    // For each channel block, the index in the channels at which this channel block starts
      std::vector<double> itsChanBlockFreqs;
      std::vector<std::vector<string> > itsDirections; // For each direction, a vector of patches
      std::vector<std::unique_ptr<Constraint> > itsConstraints;

      std::vector<double>   itsWeightsPerAntenna;

      UVWFlagger       itsUVWFlagStep;
      ResultStep::ShPtr itsDataResultStep; // Result step for data after UV-flagging
      std::vector<Predict>     itsPredictSteps;
      std::vector<MultiResultStep::ShPtr> itsResultSteps; // For each directions, a multiresultstep with all times

      NSTimer          itsTimer;
      NSTimer          itsTimerPredict;
      NSTimer          itsTimerSolve;
      NSTimer          itsTimerWrite;
      double           itsCoreConstraint;
      std::vector<std::set<std::string>> itsAntennaConstraint;
      double           itsSmoothnessConstraint;
      double           itsScreenCoreConstraint;
      MultiDirSolver   itsMultiDirSolver;
      bool itsFullMatrixMinimalization;
      bool itsApproximateTEC;
			bool itsSubtract;
      bool itsSaveFacets;
      std::string itsStatFilename;
			std::unique_ptr<ThreadPool> itsThreadPool;
      std::unique_ptr<FacetPredict> itsFacetPredictor;
      std::unique_ptr<std::ofstream> itsStatStream;
    };

  } //# end namespace
}

#endif
