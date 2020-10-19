// DDE.h: DPPP step class to calibrate direction dependent gains
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
/// @brief DPPP step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DPPP_DDECAL_H
#define DPPP_DDECAL_H

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
#include "../DPPP/SolutionInterval.h"

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

namespace aocommon {
class ThreadPool;
}  // namespace aocommon

namespace DP3 {

class ParameterSet;

namespace DPPP {

class IDGPredict;

typedef std::vector<Patch::ConstPtr> PatchList;
typedef std::pair<size_t, size_t> Baseline;

/// @brief This class is a DPStep class to calibrate (direction independent)
/// gains.
class DDECal : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DDECal(DPInput*, const ParameterSet&, const std::string& prefix);

  virtual ~DDECal();

  /// Create an DDECal object using the given parset.
  static DPStep::ShPtr makeStep(DPInput*, const ParameterSet&,
                                const std::string&);

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer&);

  void checkMinimumVisibilities(size_t bufferIndex);

  void flagChannelBlock(size_t cbIndex, size_t bufferIndex);

  /// Call the actual solver (called once per solution interval)
  void doSolve();

  /// Initialize H5parm-file
  void initH5parm();

  /// Write out the solutions
  void writeSolutions();

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  void initializeConstraints(const ParameterSet& parset, const string& prefix);
  void initializeIDG(const ParameterSet& parset, const string& prefix);
  void initializePredictSteps(const ParameterSet& parset, const string& prefix);

  void doPrepare(const DPBuffer& bufin, size_t sol_int, size_t step);

  /// Initialize solutions
  void initializeScalarSolutions(size_t);

  void initializeFullMatrixSolutions(size_t);

  /// Convert itsDirections to a vector of strings like "[Patch1, Patch2]"
  /// Used for setting source names.
  std::vector<std::string> getDirectionNames();

  void subtractCorrectedModel(bool fullJones, size_t bufferIndex);

  DPInput* itsInput;
  std::string itsName;
  /// The solution intervals that are buffered, limited by solintcount
  std::vector<SolutionInterval> sol_ints_;

  bool itsUseModelColumn;

  /// The time of the current buffer (in case of solint, average time)
  double itsAvgTime;

  /// For each time, for each channel block, a vector of size nAntennas *
  /// nDirections
  std::vector<std::vector<std::vector<casacore::DComplex>>> itsSols;
  std::vector<size_t> itsNIter,  // Number of iterations taken
      itsNApproxIter;

  /// For each time, for each constraint, a vector of results (e.g. tec and
  /// phase)
  std::vector<std::vector<std::vector<Constraint::Result>>> itsConstraintSols;

  std::string itsH5ParmName;
  H5Parm itsH5Parm;
  std::string itsParsetString;  ///< Parset, for logging in H5Parm

  GainCal::CalType itsMode;
  bool itsPropagateSolutions;
  bool itsPropagateConvergedOnly;
  bool itsFlagUnconverged;
  bool itsFlagDivergedOnly;
  bool itsUseIDG;
  bool itsOnlyPredict;
  size_t itsTimeStep;
  size_t itsSolInt;  ///< Number of timeslots to store per solution interval
  size_t itsSolIntCount;  ///< Number of solution intervals to buffer
  size_t itsNSolInts;     ///< Total number of created solution intervals
  double itsMinVisRatio;
  /// The current amount of timeslots on the solution interval
  size_t itsStepInSolInt;
  /// The current amount of solution intervals in sol_ints_
  size_t itsStepInSolInts;
  size_t itsNChan;
  /// For each channel block, the nr of unflagged vis and the total nr of vis.
  std::vector<std::pair<size_t, size_t>> itsVisInInterval;
  /// For each channel block, the index in the channels at which this channel
  /// block starts.
  std::vector<size_t> itsChanBlockStart;
  std::vector<double> itsChanBlockFreqs;
  /// For each direction, a vector of patches.
  std::vector<std::vector<string>> itsDirections;
  std::vector<std::unique_ptr<Constraint>> itsConstraints;

  std::vector<double> itsWeightsPerAntenna;

  UVWFlagger itsUVWFlagStep;
  /// Result step for data after UV-flagging
  ResultStep::ShPtr itsDataResultStep;
  std::vector<DPStep*> itsSteps;
  /// For each directions, a multiresultstep with all times.
  std::vector<MultiResultStep::ShPtr> itsResultSteps;

  NSTimer itsTimer;
  NSTimer itsTimerPredict;
  NSTimer itsTimerSolve;
  NSTimer itsTimerWrite;
  double itsCoreConstraint;
  std::vector<std::set<std::string>> itsAntennaConstraint;
  double itsSmoothnessConstraint;
  double itsScreenCoreConstraint;
  MultiDirSolver itsMultiDirSolver;
  bool itsFullMatrixMinimalization;
  bool itsApproximateTEC;
  bool itsSubtract;
  std::string itsStatFilename;
  std::unique_ptr<aocommon::ThreadPool> itsThreadPool;
  std::unique_ptr<std::ofstream> itsStatStream;
};

}  // namespace DPPP
}  // namespace DP3

#endif
