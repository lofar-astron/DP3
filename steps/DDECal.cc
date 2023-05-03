// DDECal.cc: DP3 step class to do a direction dependent gain calibration
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema & André Offringa

#include "DDECal.h"

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <utility>

#include <casacore/casa/Quanta/Quantum.h>

#include <xtensor/xview.hpp>

#include <aocommon/matrix2x2.h>

#include <schaapcommon/facets/facet.h>

#include <Version.h>

#include <dp3/base/DP3.h>

#include "../base/SourceDBUtil.h"

#include "../common/StreamUtil.h"

#include "../ddecal/SolverFactory.h"
#ifdef ENABLE_SCREENFITTER
#include "../ddecal/constraints/ScreenConstraint.h"
#endif
#include "../ddecal/SolutionResampler.h"
#include "../ddecal/constraints/SmoothnessConstraint.h"
#include "../ddecal/gain_solvers/SolveData.h"
#include "../ddecal/gain_solvers/SolverTools.h"
#include "../ddecal/linear_solvers/LLSSolver.h"

#include "IDGPredict.h"
#include "MsColumnReader.h"
#include "Predict.h"

using aocommon::FitsReader;

using schaapcommon::facets::Facet;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::ddecal::LLSSolver;
using dp3::ddecal::LLSSolverType;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

DDECal::DDECal(const common::ParameterSet& parset, const std::string& prefix)
    : itsSettings(parset, prefix),
      itsAvgTime(0),
      itsSols(),
      itsSolutionWriter(itsSettings.h5parm_name),
      itsRequestedSolInt(itsSettings.solution_interval),
      itsSolutionsPerDirection(itsSettings.n_solutions_per_direction),
      itsSolIntCount(1),
      itsFirstSolutionIndex(0),
      itsBufferedSolInts(0),
      itsNChan(itsSettings.n_channels),
      itsUVWFlagStep(parset, prefix, Step::MsType::kRegular),
      itsStoreSolutionInBuffer(parset.getBool(prefix + "storebuffer", false)),
      itsSolver(ddecal::CreateSolver(itsSettings, parset, prefix)),
      itsStatStream() {
  if (!itsSettings.stat_filename.empty()) {
    itsStatStream = std::make_unique<std::ofstream>(itsSettings.stat_filename);
  }

  // Initialize steps
  initializeColumnReaders(parset, prefix);
  initializeIDG(parset, prefix);
  initializePredictSteps(parset, prefix);

  if (itsDirections.size() == 0) {
    throw std::runtime_error(
        "DDECal initialized with 0 directions: something is wrong with your "
        "parset or your sourcedb");
  }

  // To be compatible with IDGPredict and a predict with a
  // provided source_db, fill itsSolutionPerDirection if it's
  // empty at this stage.
  if (itsSolutionsPerDirection.empty()) {
    itsSolutionsPerDirection.assign(itsDirections.size(), 1);
  }

  const size_t max_n_solutions_per_direction = *std::max_element(
      itsSolutionsPerDirection.begin(), itsSolutionsPerDirection.end());

  for (const size_t val : itsSolutionsPerDirection) {
    if (itsRequestedSolInt % val) {
      throw std::runtime_error(
          "Values in ddecal.solutions_per_direction should be integer divisors "
          "of solint");
    }
  }

  // Since info().ntime() might not be set at this stage, also throw an error
  // in case itsRequestedSolInt equals 0, and any of the solves per direction >
  // 1
  const size_t actual_solution_interval =
      (max_n_solutions_per_direction > 1)
          ? itsRequestedSolInt / max_n_solutions_per_direction
          : max_n_solutions_per_direction;
  if (actual_solution_interval == 0) {
    throw std::runtime_error(
        "Maximum value in ddecal.solutions_per_direction is larger than "
        "dde.solint value OR provided dde.solint equals 0 in combination with "
        "dde.solutions_per_direction entries > 1.");
  }
}

void DDECal::initializeColumnReaders(const common::ParameterSet& parset,
                                     const string& prefix) {
  for (const std::string& col : itsSettings.model_data_columns) {
    if (itsSettings.model_data_columns.size() == 1) {
      itsDirections.emplace_back(1, "pointing");
    } else {
      itsDirections.emplace_back(1, col);
    }
    itsSteps.push_back(std::make_shared<MsColumnReader>(parset, prefix, col));
    setModelNextSteps(*itsSteps.back(), col, parset, prefix);
  }
}

void DDECal::initializeIDG(const common::ParameterSet& parset,
                           const string& prefix) {
  // TODO it would be nicer to get a new method in the DS9 reader to first get
  // names of directions and pass that to idgpredict. It will then read it
  // itself instead of DDECal having to do everything. It is better to do it all
  // in IDGPredict, so we can also make it
  if (itsSettings.idg_region_filename.empty() &&
      itsSettings.idg_image_filenames.empty()) {
    return;
  }

  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<float>>>
      readers = IDGPredict::GetReaders(itsSettings.idg_image_filenames);
  std::vector<Facet> facets = IDGPredict::GetFacets(
      itsSettings.idg_region_filename, readers.first.front());

  for (size_t i = 0; i < facets.size(); ++i) {
    std::string dir_name = "dir" + std::to_string(i);
    if (!facets[i].DirectionLabel().empty()) {
      dir_name = facets[i].DirectionLabel();
    }
    itsDirections.emplace_back(1, dir_name);

    itsSteps.push_back(std::make_shared<IDGPredict>(
        parset, prefix, readers, std::vector<Facet>{facets[i]}));
    setModelNextSteps(*itsSteps.back(), facets[i].DirectionLabel(), parset,
                      prefix);
  }
}

void DDECal::initializePredictSteps(const common::ParameterSet& parset,
                                    const string& prefix) {
  std::vector<std::vector<std::string>> directions =
      base::MakeDirectionList(itsSettings.directions, itsSettings.source_db);

  for (std::vector<std::string>& direction : directions) {
    itsSteps.push_back(std::make_shared<Predict>(parset, prefix, direction));
    setModelNextSteps(*itsSteps.back(), direction.front(), parset, prefix);
    itsDirections.push_back(std::move(direction));
  }
}

void DDECal::setModelNextSteps(Step& step, const std::string& direction,
                               const common::ParameterSet& parset,
                               const string& prefix) const {
  std::string step_names_key = prefix + "modelnextsteps." + direction;
  if (!parset.isDefined(step_names_key)) {
    step_names_key = prefix + "modelnextsteps";  // Fall back setting.
  }

  if (parset.isDefined(step_names_key)) {
    Step::ShPtr first_step = base::MakeStepsFromParset(
        parset, "", step_names_key, "", false, steps::Step::MsType::kRegular);

    if (first_step) {
      step.setNextStep(first_step);
    }
  }
}

void DDECal::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  itsUVWFlagStep.updateInfo(infoIn);
  itsThreadPool = std::make_unique<aocommon::ThreadPool>(getInfo().nThreads());

  if (itsRequestedSolInt == 0) {
    itsRequestedSolInt = info().ntime();
  }

  // Update info for substeps and set other required parameters
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsSteps[dir]->setInfo(infoIn);

    if (auto s = std::dynamic_pointer_cast<Predict>(itsSteps[dir])) {
      s->SetThreadData(*itsThreadPool, &itsMeasuresMutex);
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(itsSteps[dir])) {
      itsSolIntCount =
          std::max(itsSolIntCount,
                   s->GetBufferSize() / itsSteps.size() / itsRequestedSolInt);
      // We increment by one so the IDGPredict will not flush in its process
      s->SetBufferSize(itsRequestedSolInt * itsSolIntCount + 1);
    } else if (!std::dynamic_pointer_cast<MsColumnReader>(itsSteps[dir])) {
      throw std::runtime_error("DDECal received an invalid first model step");
    }
  }

  itsSolver->SetNThreads(getInfo().nThreads());

  if (!itsUVWFlagStep.isDegenerate()) {
    itsDataResultStep = std::make_shared<ResultStep>();
    itsUVWFlagStep.setNextStep(itsDataResultStep);
  }
  if (!itsUVWFlagStep.isDegenerate() || itsSettings.min_vis_ratio > 0.0) {
    // The UVWFlagger and/or flagChannelBlock() may modify the flags.
    // Create a buffer for storing the original flags.
    itsOriginalFlags.resize({itsSolIntCount, itsRequestedSolInt,
                             infoIn.nbaselines(), infoIn.nchan(),
                             infoIn.ncorr()});
  }

  // For each sub step chain, add a resultstep and get required fieilds.
  itsRequiredFields.reserve(itsSteps.size());
  itsResultSteps.resize(itsSteps.size());
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsRequiredFields.emplace_back(base::GetChainRequiredFields(itsSteps[dir]));

    itsResultSteps[dir] =
        std::make_shared<MultiResultStep>(itsRequestedSolInt * itsSolIntCount);

    // Add the resultstep to the end of the model next steps
    std::shared_ptr<Step> step = itsSteps[dir];
    while (step->getNextStep()) {
      step = step->getNextStep();
    }
    step->setNextStep(itsResultSteps[dir]);
  }

  if (itsNChan == 0 || itsNChan > info().nchan()) {
    itsNChan = info().nchan();
  }

  // Create lists with used antenna indices, similarly to
  // DPInfo::removeUnusedAnt.
  itsAntennas1.resize(info().getAnt1().size());
  itsAntennas2.resize(info().getAnt2().size());
  for (size_t i = 0; i < itsAntennas1.size(); ++i) {
    itsAntennas1[i] = info().antennaMap()[info().getAnt1()[i]];
    itsAntennas2[i] = info().antennaMap()[info().getAnt2()[i]];
  }

  // Fill antenna info in H5Parm, need to convert from casa types to std types
  // Fill in metadata for all antennas, also those that may be filtered out.
  std::vector<std::string> antennaNames(info().antennaNames().size());
  std::vector<std::array<double, 3>> antennaPos(info().antennaPos().size());
  for (unsigned int i = 0; i < info().antennaNames().size(); ++i) {
    antennaNames[i] = info().antennaNames()[i];
    casacore::Quantum<casacore::Vector<double>> pos =
        info().antennaPos()[i].get("m");
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  itsSolutionWriter.AddAntennas(antennaNames, antennaPos);

  size_t nSolTimes =
      (info().ntime() + itsRequestedSolInt - 1) / itsRequestedSolInt;
  size_t nChannelBlocks = info().nchan() / itsNChan;
  itsSols.resize(nSolTimes);
  itsNIter.resize(nSolTimes);
  itsNApproxIter.resize(nSolTimes);
  itsConstraintSols.resize(nSolTimes);
  itsVisInInterval.assign(nChannelBlocks, std::pair<size_t, size_t>(0, 0));

  itsChanBlockStart.resize(nChannelBlocks + 1);
  itsChanBlockFreqs.resize(nChannelBlocks);
  itsChanBlockStart.front() = 0;
  for (size_t chBlock = 0; chBlock != nChannelBlocks; ++chBlock) {
    itsChanBlockStart[chBlock + 1] =
        (chBlock + 1) * info().nchan() / nChannelBlocks;
    const size_t blockSize =
        itsChanBlockStart[chBlock + 1] - itsChanBlockStart[chBlock];
    const double* freqStart =
        info().chanFreqs().data() + itsChanBlockStart[chBlock];
    const double meanFreq =
        std::accumulate(freqStart, freqStart + blockSize, 0.0) / blockSize;
    itsChanBlockFreqs[chBlock] = meanFreq;
  }

  itsWeightsPerAntenna.assign(
      itsChanBlockFreqs.size() * info().antennaUsed().size(), 0.0);

  itsSourceDirections.reserve(itsSteps.size());
  for (const std::shared_ptr<ModelDataStep>& s : itsSteps) {
    itsSourceDirections.push_back(s->GetFirstDirection());
  }

  // Prepare positions and names for the used antennas only.
  std::vector<std::string> used_antenna_names;
  std::vector<std::array<double, 3>> used_antenna_positions;
  used_antenna_names.reserve(info().antennaUsed().size());
  used_antenna_positions.reserve(info().antennaUsed().size());
  for (const int& ant : info().antennaUsed()) {
    used_antenna_names.push_back(info().antennaNames()[ant]);
    used_antenna_positions.push_back(antennaPos[ant]);
  }

  for (ddecal::SolverBase* solver : itsSolver->ConstraintSolvers()) {
    InitializeSolverConstraints(*solver, itsSettings, used_antenna_positions,
                                used_antenna_names, itsSolutionsPerDirection,
                                itsSourceDirections, itsChanBlockFreqs);
  }

  size_t nSt = info().antennaUsed().size();
  // Give renumbered antennas to solver
  itsSolver->Initialize(nSt, itsSolutionsPerDirection, nChannelBlocks);

  for (size_t i = 0; i < nSolTimes; ++i) {
    itsSols[i].resize(nChannelBlocks);
  }
}

void DDECal::show(std::ostream& os) const {
  os << "DDECal " << itsSettings.name << '\n'
     << "  mode (constraints):  " << ToString(itsSettings.mode) << '\n'
     << "  algorithm:           "
     << ddecal::ToString(itsSettings.solver_algorithm) << '\n'
     << "  H5Parm:              " << itsSettings.h5parm_name << '\n'
     << "  write sol to buffer: " << std::boolalpha << itsStoreSolutionInBuffer
     << '\n'
     << "  solint:              " << itsRequestedSolInt << '\n'
     << "  nchan:               " << itsNChan << '\n'
     << "  directions:          " << itsDirections << '\n'
     << "  sols per direction:  " << itsSolutionsPerDirection << '\n';
  if (itsSettings.min_vis_ratio != 0.0) {
    os << "  min visib. ratio:    " << itsSettings.min_vis_ratio << '\n';
  }
  os << "  tolerance:           " << itsSolver->GetAccuracy() << '\n'
     << "  max iter:            " << itsSolver->GetMaxIterations() << '\n'
     << "  flag unconverged:    " << std::boolalpha
     << itsSettings.flag_unconverged << '\n'
     << "     diverged only:    " << std::boolalpha
     << itsSettings.flag_diverged_only << '\n'
     << "  propagate solutions: " << std::boolalpha
     << itsSettings.propagate_solutions << '\n'
     << "       converged only: " << std::boolalpha
     << itsSettings.propagate_converged_only << '\n'
     << "  detect stalling:     " << std::boolalpha
     << itsSolver->GetDetectStalling() << '\n'
     << "  step size:           " << itsSolver->GetStepSize() << '\n';
  ShowConstraintSettings(os, itsSettings);
  os << "  approximate fitter:  " << itsSettings.approximate_tec << '\n'
     << "  only predict:        " << itsSettings.only_predict << '\n'
     << "  subtract model:      " << itsSettings.subtract << '\n';
  for (unsigned int i = 0; i < itsSteps.size(); ++i) {
    std::shared_ptr<Step> step = itsSteps[i];
    os << "Model steps for direction " << itsDirections[i][0] << '\n';
    do {
      step->show(os);
    } while (step = step->getNextStep());
    os << '\n';
  }
  itsUVWFlagStep.show(os);
}

void DDECal::showTimings(std::ostream& os, double duration) const {
  double totaltime = itsTimer.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " DDECal " << itsSettings.name << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerPredict.getElapsed(), totaltime);
  os << " of it spent in predict" << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerSolve.getElapsed(), totaltime);
  os << " of it spent in estimating gains and computing residuals" << '\n';

  itsSolver->GetTimings(os, itsTimerSolve.getElapsed());

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerWrite.getElapsed(), totaltime);
  os << " of it spent in writing gain solutions to disk" << '\n';

  os << "          ";
  os << "Substeps taken:" << '\n';
  for (auto& step : itsSteps) {
    os << "          ";
    step->showTimings(os, duration);
  }

  os << "Iterations taken: [";
  for (size_t i = 0; i < itsNIter.size() - 1; ++i) {
    os << itsNIter[i];
    if (itsNApproxIter[i] != 0) os << '|' << itsNApproxIter[i];
    os << ",";
  }
  os << itsNIter[itsNIter.size() - 1];
  if (itsNApproxIter[itsNIter.size() - 1] != 0)
    os << '|' << itsNApproxIter[itsNIter.size() - 1];
  os << "]" << '\n';
}

void DDECal::InitializeSolutions(size_t buffer_index) {
  const size_t solution_index = itsFirstSolutionIndex + buffer_index;
  assert(solution_index < itsSols.size());

  bool propagate_solutions =
      solution_index > 0 && itsSettings.propagate_solutions;
  if (propagate_solutions &&
      itsNIter[solution_index - 1] > itsSolver->GetMaxIterations() &&
      itsSettings.propagate_converged_only) {
    propagate_solutions = false;
  }

  if (propagate_solutions) {
    // Initialize solutions with those of the previous step.
    itsSols[solution_index] = itsSols[solution_index - 1];
  } else {
    const size_t n_solutions = std::accumulate(
        itsSolutionsPerDirection.begin(), itsSolutionsPerDirection.end(), 0u);
    const size_t n_solution_values = n_solutions * info().antennaUsed().size() *
                                     itsSolver->NSolutionPolarizations();

    if (itsSolver->NSolutionPolarizations() == 4) {
      // Initialize solutions with unity matrix [1 0 ; 0 1].
      for (std::vector<casacore::DComplex>& solution_vector :
           itsSols[solution_index]) {
        solution_vector.resize(n_solution_values);
        for (size_t i = 0; i < n_solution_values; i += 4) {
          solution_vector[i + 0] = 1.0;
          solution_vector[i + 1] = 0.0;
          solution_vector[i + 2] = 0.0;
          solution_vector[i + 3] = 1.0;
        }
      }
    } else {
      // Initialize solutions with 1.
      for (std::vector<casacore::DComplex>& solution_vector :
           itsSols[solution_index]) {
        solution_vector.assign(n_solution_values, 1.0);
      }
    }
  }
}

void DDECal::flagChannelBlock(size_t cbIndex, size_t bufferIndex) {
  const size_t nBl = info().nbaselines();
  const size_t nChanBlocks = itsChanBlockFreqs.size();
  const size_t channel_begin = itsChanBlockStart[cbIndex];
  const size_t channel_end = itsChanBlockStart[cbIndex + 1];
  // Set the antenna-based weights to zero
  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t ant1 = info().antennaMap()[info().getAnt1()[bl]];
    size_t ant2 = info().antennaMap()[info().getAnt2()[bl]];
    itsWeightsPerAntenna[ant1 * nChanBlocks + cbIndex] = 0.0;
    itsWeightsPerAntenna[ant2 * nChanBlocks + cbIndex] = 0.0;
  }
  // Set the visibility flags to true. SolverTools::AssignAndWeight will write
  // zeroes to the weighted data if it is flagged.
  for (std::unique_ptr<DPBuffer>& buffer :
       itsSolIntBuffers[bufferIndex].DataBuffers()) {
    xt::view(buffer->GetFlags(), xt::all(),
             xt::range(channel_begin, channel_end), xt::all()) = true;
  }
}

void DDECal::checkMinimumVisibilities(size_t bufferIndex) {
  for (size_t cb = 0; cb != itsChanBlockFreqs.size(); ++cb) {
    double fraction =
        double(itsVisInInterval[cb].first) / itsVisInInterval[cb].second;
    if (fraction < itsSettings.min_vis_ratio) flagChannelBlock(cb, bufferIndex);
  }
}

void DDECal::doSolve() {
  std::vector<std::string> direction_names;

  for (size_t dir = 0; dir < itsDirections.size(); ++dir) {
    if (auto s = dynamic_cast<IDGPredict*>(itsSteps[dir].get())) {
      itsTimerPredict.start();
      s->flush();
      itsTimerPredict.stop();
    }
    // TODO(AST-1240): Use real names instead of increasing numbers.
    direction_names.push_back(std::to_string(dir));
    for (size_t i = 0; i < itsResultSteps[dir]->size(); ++i) {
      const size_t sol_int = i / itsRequestedSolInt;
      const size_t timestep = i % itsRequestedSolInt;
      itsSolIntBuffers[sol_int].DataBuffers()[timestep]->MoveData(
          itsResultSteps[dir]->get()[i], "", direction_names.back());
    }
  }

  std::vector<ddecal::SolverBase*> solvers = itsSolver->ConstraintSolvers();
  const size_t n_channel_blocks = itsChanBlockFreqs.size();
  const size_t n_antennas = info().antennaUsed().size();

  // Declare weighted_buffers outside the loop, so we can reuse its memory.
  std::vector<base::DPBuffer> weighted_buffers(itsRequestedSolInt);

  // DDECal requires the unweighted model model when the model is subtracted
  // after calibration. Since the model data can be large, memory allocation
  // for this optional feature is done conditionally.
  const bool subtract_corrected_model =
      itsSettings.subtract || itsSettings.only_predict;

  for (size_t i = 0; i < itsSolIntBuffers.size(); ++i) {
    const size_t solution_index = itsFirstSolutionIndex + i;
    assert(solution_index < itsSols.size());

    ddecal::SolverBase::SolveResult solveResult;
    if (!itsSettings.only_predict) {
      if (itsSettings.min_vis_ratio > 0.0) {
        checkMinimumVisibilities(i);
      }

      for (ddecal::SolverBase* solver : solvers) {
        for (const std::unique_ptr<ddecal::Constraint>& constraint :
             solver->GetConstraints()) {
          constraint->SetWeights(itsWeightsPerAntenna);
        }
      }

      std::vector<std::unique_ptr<DPBuffer>>& unweighted_buffers =
          itsSolIntBuffers[i].DataBuffers();

      // The last solution interval can be smaller.
      weighted_buffers.resize(unweighted_buffers.size());

      ddecal::AssignAndWeight(unweighted_buffers, direction_names,
                              weighted_buffers, subtract_corrected_model);

      InitializeSolutions(i);

      itsTimerSolve.start();

      const ddecal::SolveData solve_data(
          weighted_buffers, direction_names, n_channel_blocks, n_antennas,
          itsSolutionsPerDirection, itsAntennas1, itsAntennas2);

      solveResult = itsSolver->Solve(solve_data, itsSols[solution_index],
                                     itsAvgTime / itsRequestedSolInt,
                                     itsStatStream.get());

      itsTimerSolve.stop();

      itsNIter[solution_index] = solveResult.iterations;
      itsNApproxIter[solution_index] = solveResult.constraint_iterations;
    }

    if (subtract_corrected_model) {
      subtractCorrectedModel(i);
    }

    // Check for nonconvergence and flag if desired. Unconverged solutions are
    // identified by the number of iterations being one more than the max
    // allowed number
    if (solveResult.iterations > itsSolver->GetMaxIterations() &&
        itsSettings.flag_unconverged) {
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          if (itsSettings.flag_diverged_only) {
            // Set weights with negative values (indicating unconverged
            // solutions that diverged) to zero (all other unconverged
            // solutions remain unflagged)
            for (double& weight : result.weights) {
              if (weight < 0.) weight = 0.;
            }
          } else {
            // Set all weights to zero
            result.weights.assign(result.weights.size(), 0.);
          }
        }
      }
    } else {
      // Set any negative weights (indicating unconverged solutions that
      // diverged) to one (all other unconverged solutions are unflagged
      // already)
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          for (double& weight : result.weights) {
            if (weight < 0.0) weight = 1.0;
          }
        }
      }
    }

    // Store constraint solutions if any constaint has a non-empty result
    bool someConstraintHasResult = false;
    for (const auto& constraint_results : solveResult.results) {
      if (!constraint_results.empty()) {
        someConstraintHasResult = true;
        break;
      }
    }
    if (someConstraintHasResult) {
      itsConstraintSols[solution_index] = solveResult.results;
    }

    // Store calibration solution for later calibration application steps.
    if (itsStoreSolutionInBuffer) {
      itsSolIntBuffers[i].DataBuffers().front()->SetSolution(
          itsSols[solution_index]);
    }
  }

  itsTimer.stop();

  for (size_t i = 0; i < itsSolIntBuffers.size(); ++i) {
    for (size_t step = 0; step < itsSolIntBuffers[i].Size(); ++step) {
      if (itsOriginalFlags.size() > 0) {
        // Restore original flags, if the UVWFlagger or flagChannelBlock ran.
        itsSolIntBuffers[i].DataBuffers()[step]->GetFlags() = xt::view(
            itsOriginalFlags, i, step, xt::all(), xt::all(), xt::all());
      }
      // Push data (possibly changed) to next step
      getNextStep()->process(itsSolIntBuffers[i][step]);
    }
  }

  itsTimer.start();
}

bool DDECal::process(std::unique_ptr<DPBuffer> bufin) {
  itsTimer.start();

  // Create a new solution interval if needed
  if (itsSolIntBuffers.empty() ||
      itsSolIntBuffers.back().Size() == itsRequestedSolInt) {
    itsSolIntBuffers.emplace_back(itsRequestedSolInt);
  }

  const size_t currentIntervalIndex = itsSolIntBuffers.back().Size();
  itsSolIntBuffers[itsBufferedSolInts].PushBack(std::move(bufin));
  doPrepare(itsSolIntBuffers.back().DataBuffers()[currentIntervalIndex],
            itsBufferedSolInts, currentIntervalIndex);

  if (currentIntervalIndex + 1 == itsRequestedSolInt) {
    ++itsBufferedSolInts;
  }

  if (itsBufferedSolInts == itsSolIntCount) {
    doSolve();

    // Clean up, prepare for next iteration
    itsFirstSolutionIndex += itsSolIntBuffers.size();
    itsAvgTime = 0;
    itsBufferedSolInts = 0;
    itsVisInInterval.assign(itsVisInInterval.size(),
                            std::pair<size_t, size_t>(0, 0));
    itsWeightsPerAntenna.assign(itsWeightsPerAntenna.size(), 0.0);

    for (size_t dir = 0; dir < itsResultSteps.size(); ++dir) {
      itsResultSteps[dir]->clear();
    }
    itsSolIntBuffers.clear();
  }

  itsTimer.stop();

  return false;
}

void DDECal::doPrepare(std::unique_ptr<DPBuffer>& bufin, size_t sol_int,
                       size_t step) {
  if (itsOriginalFlags.size() > 0) {
    // Save the original flags, so DDECal can restore any flags changed by the
    // UVWFlagger/flagChannelBlock before passing the buffer to the next step.
    xt::view(itsOriginalFlags, sol_int, step, xt::all(), xt::all(), xt::all()) =
        bufin->GetFlags();
  }
  if (!itsUVWFlagStep.isDegenerate()) {
    itsUVWFlagStep.process(std::move(bufin));
    bufin = itsDataResultStep->extract();
  }

  itsTimerPredict.start();

  itsThreadPool->For(0, itsSteps.size(), [&](size_t dir, size_t) {
    itsSteps[dir]->process(
        std::make_unique<DPBuffer>(*bufin, itsRequiredFields[dir]));
  });

  // Handle weights and flags
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = 4;

  size_t nchanblocks = itsChanBlockFreqs.size();

  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t chanblock = 0, ant1 = info().antennaMap()[info().getAnt1()[bl]],
           ant2 = info().antennaMap()[info().getAnt2()[bl]];
    for (size_t ch = 0; ch < nCh; ++ch) {
      if (ch == itsChanBlockStart[chanblock + 1]) {
        chanblock++;
      }
      for (size_t cr = 0; cr < nCr; ++cr) {
        itsVisInInterval[chanblock].second++;  // total nr of vis
        if (!bufin->GetFlags()(bl, ch, cr)) {
          // Add this weight to both involved antennas
          const double weight = bufin->GetWeights()(bl, ch, cr);
          if (weight != 0.0) {
            itsWeightsPerAntenna[ant1 * nchanblocks + chanblock] += weight;
            itsWeightsPerAntenna[ant2 * nchanblocks + chanblock] += weight;
            itsVisInInterval[chanblock].first++;  // unflagged nr of vis
          }
        }
      }
    }
  }

  const double weightFactor =
      1. / (nCh * (info().antennaUsed().size() - 1) * nCr * itsRequestedSolInt);
  for (double& weight : itsWeightsPerAntenna) {
    weight *= weightFactor;
  }

  itsTimerPredict.stop();

  itsAvgTime += itsAvgTime + bufin->getTime();
}

void DDECal::WriteSolutions() {
  itsTimerWrite.start();

  // Create antenna info for H5Parm, used antennas only.
  std::vector<std::string> used_antenna_names;
  used_antenna_names.reserve(info().antennaUsed().size());
  for (size_t used_antenna : info().antennaUsed()) {
    used_antenna_names.emplace_back(info().antennaNames()[used_antenna]);
  }

  const std::string history = "CREATE by " + DP3Version::AsString() + "\n" +
                              "step " + itsSettings.name + " in parset: \n" +
                              itsSettings.parset_string;

  const size_t n_directions = itsSolutionsPerDirection.size();
  const size_t n_solutions = std::accumulate(
      itsSolutionsPerDirection.begin(), itsSolutionsPerDirection.end(), 0u);
  if (n_solutions == n_directions) {
    itsSolutionWriter.Write(itsSols, itsConstraintSols, info().startTime(),
                            info().timeInterval() * itsRequestedSolInt,
                            itsSettings.mode, used_antenna_names,
                            itsSourceDirections, itsDirections,
                            info().chanFreqs(), itsChanBlockFreqs, history);
  } else {
    const size_t n_antennas = used_antenna_names.size();
    ddecal::SolutionResampler resampler(itsSolutionsPerDirection, n_antennas,
                                        itsSolver->NSolutionPolarizations(),
                                        itsRequestedSolInt);
    const size_t solution_interval =
        itsRequestedSolInt / resampler.GetNrSubSteps();

    std::vector<std::vector<std::vector<casacore::DComplex>>>
        upsampled_solutions = resampler.Upsample(itsSols);
    itsSolutionWriter.Write(
        upsampled_solutions, itsConstraintSols, info().startTime(),
        info().timeInterval() * solution_interval, itsSettings.mode,
        used_antenna_names, itsSourceDirections, itsDirections,
        info().chanFreqs(), itsChanBlockFreqs, history);
  }

  itsTimerWrite.stop();
}

void DDECal::finish() {
  itsTimer.start();

  if (itsSolIntBuffers.size() > 0) {
    doSolve();
  }

  if (!itsSettings.only_predict) WriteSolutions();

  itsSolIntBuffers.clear();
  itsTimer.stop();

  // Let the next steps finish.
  getNextStep()->finish();
}

void DDECal::subtractCorrectedModel(size_t bufferIndex) {
  // The original unweighted data and model data are still in the solution
  // interval. Here we apply the solutions to all the model data directions and
  // subtract them from the data.

  const size_t n_directions = itsSolutionsPerDirection.size();
  const size_t n_solutions = std::accumulate(
      itsSolutionsPerDirection.begin(), itsSolutionsPerDirection.end(), 0u);
  if (n_solutions != n_directions) {
    throw std::runtime_error(
        "Subtraction is not implemented for DDECal with direction dependent "
        "solution intervals");
  }

  const size_t solution_index = itsFirstSolutionIndex + bufferIndex;
  assert(solution_index < itsSols.size());
  std::vector<std::vector<std::complex<double>>>& solutions =
      itsSols[solution_index];

  assert(bufferIndex < itsSolIntBuffers.size());
  base::SolutionInterval& solution_interval = itsSolIntBuffers[bufferIndex];

  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nDir = itsDirections.size();
  for (size_t time = 0; time != solution_interval.Size(); ++time) {
    DPBuffer& data_buffer = *solution_interval.DataBuffers()[time];
    for (size_t bl = 0; bl < nBl; ++bl) {
      const size_t ant1 = info().getAnt1()[bl];
      const size_t ant2 = info().getAnt2()[bl];
      size_t chanblock = 0;

      for (size_t ch = 0; ch < nCh; ++ch) {
        if (ch == itsChanBlockStart[chanblock + 1]) {
          chanblock++;
        }

        std::complex<float>* data = &data_buffer.GetData()(bl, ch, 0);
        if (itsSettings.only_predict) {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());

          for (size_t dir = 0; dir != nDir; ++dir) {
            // TODO(AST-1240): Use real direction name.
            const std::string name = std::to_string(dir);
            const std::complex<float>* model_data =
                &data_buffer.GetData(name)(bl, ch, 0);
            value += aocommon::MC2x2(model_data);
          }
          for (size_t cr = 0; cr < 4; ++cr) data[cr] = value[cr];
        } else {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());
          for (size_t dir = 0; dir != nDir; ++dir) {
            // TODO(AST-1240): Use real direction name.
            const std::string name = std::to_string(dir);
            const std::complex<float>* model_data =
                &data_buffer.GetData(name)(bl, ch, 0);
            if (itsSolver->NSolutionPolarizations() == 4) {
              const aocommon::MC2x2 sol1(
                  &solutions[chanblock][(ant1 * nDir + dir) * 4]);
              const aocommon::MC2x2 sol2(
                  &solutions[chanblock][(ant2 * nDir + dir) * 4]);
              value +=
                  sol1.Multiply(aocommon::MC2x2(model_data)).MultiplyHerm(sol2);
            } else if (itsSolver->NSolutionPolarizations() == 2) {
              const aocommon::MC2x2 sol1(
                  solutions[chanblock][(ant1 * nDir + dir) * 2 + 0], 0.0, 0.0,
                  solutions[chanblock][(ant1 * nDir + dir) * 2 + 1]);
              const aocommon::MC2x2 sol2(
                  solutions[chanblock][(ant2 * nDir + dir) * 2 + 0], 0.0, 0.0,
                  solutions[chanblock][(ant2 * nDir + dir) * 2 + 1]);
              value +=
                  sol1.Multiply(aocommon::MC2x2(model_data)).MultiplyHerm(sol2);
            } else {
              std::complex<double> solfactor(
                  solutions[chanblock][ant1 * nDir + dir] *
                  std::conj(solutions[chanblock][ant2 * nDir + dir]));
              for (size_t cr = 0; cr < 4; ++cr)
                value[cr] += solfactor * std::complex<double>(model_data[cr]);
            }
          }
          for (size_t cr = 0; cr < 4; ++cr) data[cr] -= value[cr];
        }
      }  // channel loop
    }    // bl loop

    // Remove unweighted model data, since DDECal no longer needs it.
    // TODO(AST-1243): Making discarding model data optional.
    for (size_t dir = 0; dir != nDir; ++dir) {
      // TODO(AST-1240): Use real direction name.
      const std::string name = std::to_string(dir);
      data_buffer.RemoveData(name);
    }
  }  // time loop
}

}  // namespace steps
}  // namespace dp3
