// DDECal.h: DP3 step class to calibrate direction dependent gains
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_DDECAL_H_
#define DP3_STEPS_DDECAL_H_

#include <fstream>
#include <string>
#include <vector>

#include <aocommon/recursivefor.h>

#include <schaapcommon/h5parm/h5cache.h>
#include <schaapcommon/h5parm/jonesparameters.h>

#include "../common/ParameterSet.h"

#include "../ddecal/Settings.h"
#include "../ddecal/SolutionWriter.h"
#include "../ddecal/constraints/Constraint.h"
#include "../ddecal/gain_solvers/SolverBase.h"

#include "MultiResultStep.h"
#include "ResultStep.h"
#include "UVWFlagger.h"

namespace dp3 {
namespace steps {

/// @brief This class is a Step class to calibrate (direction dependent) gains.
class DDECal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DDECal(const common::ParameterSet& parameterSet, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  common::Fields getProvidedFields() const override {
    return (itsSettings.subtract ||
            (itsSettings.only_predict && !itsSettings.keep_model_data))
               ? kDataField
               : common::Fields();
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void checkMinimumVisibilities(size_t bufferIndex);

  void flagChannelBlock(size_t cbIndex, size_t bufferIndex);

  /// Call the actual solver (called once per solution interval)
  void doSolve();

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

 private:
  void initializeColumnReaders(const common::ParameterSet&,
                               const std::string& prefix);
  void initializeModelReuse();
  void initializeInitialSolutionsH5Parm(const common::ParameterSet& parset,
                                        const std::string& prefix);
  void initializeIDG(const common::ParameterSet& parset,
                     const std::string& prefix);
  void initializePredictSteps(const common::ParameterSet& parset,
                              const std::string& prefix);

  void setModelNextSteps(Step&, const std::string& direction,
                         const common::ParameterSet& parset,
                         const std::string& prefix) const;

  void doPrepare();

  /// Initializes solutions for a new solution interval.
  /// Based on progation settings, either copies the previous solution or
  /// writes default values to the new solution.
  /// @param buffer_index Index within the current solution interval set.
  void InitializeSolutions(size_t buffer_index);

  /// Write all solutions to an H5Parm file using itsSolutionWriter.
  void WriteSolutions();

  /// Sums all model data buffers into the main data buffer.
  /// Removes the model data buffers if 'keepmodel' is false.
  /// This function implements the behavior for the "onlypredict" setting.
  void SumModels(size_t buffer_index);

  /// Applies a single solution to all directions.
  /// (Helper function for CorrectAndSubtractModels.)
  void ApplySolution(
      base::DPBuffer& buffer, size_t baseline, size_t channel,
      const std::vector<std::complex<double>>& channel_block_solutions) const;

  /// Applies the solutions to the model data for all directions.
  /// If "keepmodel" is true, overwrites the model data with the corrected model
  /// data. If "keepmodel" is false, removes the model data buffers.
  /// If "subtract" is true, subtracts all corrected model data from the main
  /// input data buffer.
  void CorrectAndSubtractModels(size_t buffer_index);

  /// Read the Jones matrix for a single time step and a single direction from
  /// one or two solution tables.
  xt::xtensor<std::complex<float>, 3> ReadJonesMatrixFromH5Parm(
      const base::Direction& direction, double timestamp,
      schaapcommon::h5parm::GainType gain_type,
      schaapcommon::h5parm::SolTab* first_soltab,
      schaapcommon::h5parm::SolTab* second_soltab);

  ddecal::Settings itsSettings;

  /// The input data buffers for the current set of solution intervals.
  /// Maximum dimensions: itsSolIntCount x itsRequestedSolInt
  std::vector<std::vector<std::unique_ptr<base::DPBuffer>>> itsInputBuffers;
  /// Original flags of the input buffers for the current solution interval.
  /// This member is only used if itsUVWFlagger is active.
  /// Dimensions: ( solution_interval x step_within_interval x baseline x
  /// channel x correlation )
  xt::xtensor<bool, 5> itsOriginalFlags;

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

  std::unique_ptr<ddecal::SolutionWriter> itsSolutionWriter;

  /// Number of timeslots to store per solution interval as requested
  /// by the user in the parset.
  size_t itsRequestedSolInt;
  size_t itsSolIntCount;  ///< Number of solution intervals to buffer
  /// Index of the first solution in the current solution interval set.
  size_t itsFirstSolutionIndex;
  size_t itsNChan;
  /// For each channel block, the nr of unflagged vis and the total nr of vis.
  std::vector<std::pair<size_t, size_t>> itsVisInInterval;
  /// For each channel block, the index in the channels at which this channel
  /// block starts.
  std::vector<size_t> itsChanBlockStart;
  std::vector<double> itsChanBlockFreqs;
  /// For each direction, a vector of patches.
  std::vector<std::vector<std::string>> itsDirections;
  /// For each direction, the name for the model data in DPBuffer.
  std::vector<std::string> itsDirectionNames;
  /// Expanded version of reusemodel patterns.
  std::vector<std::string> itsReusedDirectionNames;
  /// Maps direction indices to the cluster central direction.
  std::vector<base::Direction> itsSourceDirections;

  /// First antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas1;
  /// Second antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas2;
  std::vector<double> itsWeightsPerAntenna;

  UVWFlagger itsUVWFlagStep;
  /// Result step for data after UV-flagging
  std::shared_ptr<ResultStep> itsDataResultStep;
  /// For each direction, the first step in the chain that computes the model.
  /// When reusing model data, the step for that direction is empty/null.
  std::vector<std::shared_ptr<ModelDataStep>> itsSteps;
  /// For each direction, the required fields of the step chain.
  std::vector<common::Fields> itsRequiredFields;
  /// For each directions, a multiresultstep with all times.
  /// When reusing model data, the result step for that direction is empty/null.
  std::vector<std::shared_ptr<MultiResultStep>> itsResultSteps;

  /// Store the solution for later steps of processing in DPBuffer. Note: only
  /// works for 1 direction.
  bool itsStoreSolutionInBuffer;

  /// Stores the H5Parm file and loads all solutions into memory when the user
  /// requests the solver to use initial solutions.
  std::unique_ptr<schaapcommon::h5parm::H5Parm> itsInitialSolutions;
  std::string itsInitialSolutionsH5ParmName;
  std::vector<std::string> itsInitialSolutionsSolTab;
  std::vector<schaapcommon::h5parm::SolTab> itsSolutionTables;
  bool itsInitialSolutionsIsFullJones;
  /// Specifies the InterpolationType, MissingAntennaBehavior, and GainType for
  /// extracting the Jones parameters from itsInitialSolutions.
  /// @{
  schaapcommon::h5parm::JonesParameters::InterpolationType itsInterpolationType;
  schaapcommon::h5parm::JonesParameters::MissingAntennaBehavior
      itsMissingAntennaBehavior;
  std::vector<schaapcommon::h5parm::GainType> itsGainTypes;
  /// @}

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPredict;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerWrite;
  std::mutex itsMeasuresMutex;
  std::unique_ptr<ddecal::SolverBase> itsSolver;
  std::unique_ptr<std::ofstream> itsStatStream;
};

}  // namespace steps
}  // namespace dp3

#endif
