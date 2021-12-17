// DDE.h: DPPP step class to calibrate direction dependent gains
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DP3_DDECAL_H
#define DP3_DDECAL_H

#include "ApplyBeam.h"
#include "GainCal.h"
#include "InputStep.h"
#include "Predict.h"
#include "UVWFlagger.h"

#include "../base/DPBuffer.h"
#include "../base/BaselineSelection.h"
#include "../base/Patch.h"
#include "../base/SourceDBUtil.h"
#include "../base/SolutionInterval.h"

#include "../ddecal/Settings.h"
#include "../ddecal/SolutionWriter.h"
#include "../ddecal/constraints/Constraint.h"
#include "../ddecal/gain_solvers/SolverBase.h"

#include "../parmdb/Parm.h"

#include <schaapcommon/h5parm/h5parm.h>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <string>
#include <vector>

namespace aocommon {
class ThreadPool;
}  // namespace aocommon

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {

class IDGPredict;

typedef std::vector<base::Patch::ConstPtr> PatchList;
typedef std::pair<size_t, size_t> Baseline;

/// @brief This class is a Step class to calibrate (direction independent)
/// gains.
class DDECal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DDECal(InputStep*, const common::ParameterSet& parameterSet,
         const std::string& prefix);

  virtual ~DDECal();

  virtual bool process(const base::DPBuffer&);

  void checkMinimumVisibilities(size_t bufferIndex);

  void flagChannelBlock(size_t cbIndex, size_t bufferIndex);

  /// Call the actual solver (called once per solution interval)
  void doSolve();

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  virtual void updateInfo(const base::DPInfo&);

  virtual void show(std::ostream&) const;

  virtual void showTimings(std::ostream&, double duration) const;

  bool modifiesData() const override {
    return itsSettings.subtract || itsSettings.only_predict;
  }

 private:
  void initializeColumnReaders(const common::ParameterSet&,
                               const string& prefix);
  void initializeIDG(const common::ParameterSet& parset, const string& prefix);
  void initializePredictSteps(const common::ParameterSet& parset,
                              const string& prefix);

  void setModelNextSteps(Step&, const std::string& direction,
                         const common::ParameterSet& parset,
                         const string& prefix) const;

  void doPrepare(const base::DPBuffer& bufin, size_t sol_int, size_t step);

  /// Initialize solutions
  void InitializeScalarOrDiagonalSolutions(size_t);

  void initializeFullMatrixSolutions(size_t);

  /// Write all solutions to an H5Parm file using itsSolutionWriter.
  void WriteSolutions();

  void storeModelData(
      const std::vector<std::vector<base::DPBuffer>>& input_model_buffers);
  void subtractCorrectedModel(size_t bufferIndex);

  InputStep& itsInput;
  const ddecal::Settings itsSettings;

  /// The solution intervals that are buffered, limited by solintcount
  std::vector<base::SolutionInterval> itsSolIntBuffers;

  /// The time of the current buffer (in case of solint, average time)
  double itsAvgTime;

  /// For each time, for each channel block, a vector of size nAntennas *
  /// SolverBase::NSolutions() * nPolarizations, with nPolarizations changing
  /// fastest.
  std::vector<std::vector<std::vector<casacore::DComplex>>> itsSols;
  std::vector<size_t> itsNIter,  // Number of iterations taken
      itsNApproxIter;

  /// For each time, for each constraint, a vector of results (e.g. tec and
  /// phase)
  std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>
      itsConstraintSols;

  ddecal::SolutionWriter itsSolutionWriter;

  size_t itsTimeStep;
  /// Number of timeslots to store per solution interval as requested
  /// by the user in the parset.
  size_t itsRequestedSolInt;
  /// For each direction, a number of solutions per solution interval
  std::vector<size_t> itsSolutionsPerDirection;
  size_t itsSolIntCount;  ///< Number of solution intervals to buffer
  size_t itsNSolInts;     ///< Total number of created solution intervals
  /// The current amount of solution intervals in itsSolInts
  size_t itsBufferedSolInts;
  size_t itsNChan;
  /// For each channel block, the nr of unflagged vis and the total nr of vis.
  std::vector<std::pair<size_t, size_t>> itsVisInInterval;
  /// For each channel block, the index in the channels at which this channel
  /// block starts.
  std::vector<size_t> itsChanBlockStart;
  std::vector<double> itsChanBlockFreqs;
  /// For each direction, a vector of patches.
  std::vector<std::vector<std::string>> itsDirections;
  /// Maps direction indices to the cluster central direction.
  std::vector<base::Direction> itsSourceDirections;
  /// Normally, the solver takes the model data and modifies it, thereby
  /// destroying the original model data. This model data is used to store
  /// the result of the predictions when they are still required after
  /// solving (e.g. for subtraction).
  std::vector<std::vector<std::vector<std::complex<float>>>> itsModelData;

  /// First antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas1;
  /// Second antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas2;
  std::vector<double> itsWeightsPerAntenna;

  UVWFlagger itsUVWFlagStep;
  /// Result step for data after UV-flagging
  std::shared_ptr<ResultStep> itsDataResultStep;
  std::vector<std::shared_ptr<ModelDataStep>> itsSteps;
  /// For each directions, a multiresultstep with all times.
  std::vector<std::shared_ptr<MultiResultStep>> itsResultSteps;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPredict;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerWrite;
  std::mutex itsMeasuresMutex;
  std::unique_ptr<ddecal::SolverBase> itsSolver;
  std::unique_ptr<aocommon::ThreadPool> itsThreadPool;
  std::unique_ptr<std::ofstream> itsStatStream;
};

}  // namespace steps
}  // namespace dp3

#endif
