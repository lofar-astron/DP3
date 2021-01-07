// DDE.h: DPPP step class to calibrate direction dependent gains
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DPPP_DDECAL_H
#define DPPP_DDECAL_H

#include "../DPPP/DPInput.h"
#include "../DPPP/GainCal.h"
#include "../DPPP/DPBuffer.h"
#include "../DPPP/BaselineSelection.h"
#include "../DPPP/Patch.h"
#include "../DPPP/UVWFlagger.h"
#include "../DPPP/Predict.h"
#include "../DPPP/SourceDBUtil.h"
#include "../DPPP/ApplyBeam.h"
#include "../DPPP/SolutionInterval.h"

#include "SolverBase.h"
#include "Constraint.h"

#include <schaapcommon/h5parm/h5parm.h>

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
  DDECal(DPInput*, const ParameterSet& parameterSet, const std::string& prefix);

  virtual ~DDECal();

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

  virtual void updateInfo(const DPInfo&);

  virtual void show(std::ostream&) const;

  virtual void showTimings(std::ostream&, double duration) const;

  virtual bool modifiesData() const override {
    return itsSubtract || itsOnlyPredict;
  }

 private:
  void initializeSolver(const ParameterSet& parset, const string& prefix);
  void initializeColumnReaders(const ParameterSet&, const string& prefix);
  void initializeIDG(const ParameterSet& parset, const string& prefix);
  void initializePredictSteps(const ParameterSet& parset, const string& prefix);

  void setModelNextSteps(std::shared_ptr<DPStep>, const std::string direction,
                         const ParameterSet& parset, const string prefix);

  void doPrepare(const DPBuffer& bufin, size_t sol_int, size_t step);

  /// Initialize solutions
  void initializeScalarSolutions(size_t);

  void initializeFullMatrixSolutions(size_t);

  /// Convert itsDirections to a vector of strings like "[Patch1, Patch2]"
  /// Used for setting source names.
  std::vector<std::string> getDirectionNames();

  void storeModelData(
      const std::vector<std::vector<std::complex<float>*>>& input_model_data);
  void subtractCorrectedModel(bool fullJones, size_t bufferIndex);

  DPInput* itsInput;
  std::string itsName;
  /// The solution intervals that are buffered, limited by solintcount
  std::vector<SolutionInterval> sol_ints_;

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
  schaapcommon::h5parm::H5Parm itsH5Parm;
  std::string itsParsetString;  ///< Parset, for logging in H5Parm

  GainCal::CalType itsMode;
  bool itsPropagateSolutions;
  bool itsPropagateConvergedOnly;
  bool itsFlagUnconverged;
  bool itsFlagDivergedOnly;
  bool itsOnlyPredict;
  size_t itsTimeStep;
  size_t itsSolInt;  ///< Number of timeslots to store per solution interval
  size_t itsSolIntCount;  ///< Number of solution intervals to buffer
  size_t itsNSolInts;     ///< Total number of created solution intervals
  double itsMinVisRatio;
  /// The current amount of timeslots on the solution interval
  size_t itsStepInSolInt;
  /// The current amount of solution intervals in sol_ints_
  size_t itsBufferedSolInts;
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
  /// Normally, the solver takes the model data and modifies it, thereby
  /// destroying the original model data. This model data is used to store
  /// the result of the predictions when they are still required after
  /// solving (e.g. for subtraction).
  std::vector<std::vector<std::vector<std::complex<float>>>> itsModelData;

  std::vector<double> itsWeightsPerAntenna;

  UVWFlagger itsUVWFlagStep;
  /// Result step for data after UV-flagging
  ResultStep::ShPtr itsDataResultStep;
  std::vector<std::shared_ptr<DPStep>> itsSteps;
  /// For each directions, a multiresultstep with all times.
  std::vector<MultiResultStep::ShPtr> itsResultSteps;

  NSTimer itsTimer;
  NSTimer itsTimerPredict;
  NSTimer itsTimerSolve;
  NSTimer itsTimerWrite;
  std::mutex itsMeasuresMutex;
  double itsCoreConstraint;
  std::vector<std::set<std::string>> itsAntennaConstraint;
  double itsSmoothnessConstraint;
  double itsScreenCoreConstraint;
  std::unique_ptr<SolverBase> itsSolver;
  size_t itsPolsInSolutions;
  bool itsApproximateTEC;
  bool itsSubtract;
  std::string itsStatFilename;
  std::unique_ptr<aocommon::ThreadPool> itsThreadPool;
  std::unique_ptr<std::ofstream> itsStatStream;
};

}  // namespace DPPP
}  // namespace DP3

#endif
