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

#include <aocommon/logger.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

#include <schaapcommon/facets/facet.h>

#include <Version.h>

#include <dp3/base/DP3.h>

#include "../model/SourceDBUtil.h"

#include "../common/StreamUtil.h"

#include "../ddecal/SolverFactory.h"
#ifdef ENABLE_SCREENFITTER
#include "../ddecal/constraints/ScreenConstraint.h"
#endif
#include "../ddecal/constraints/SmoothnessConstraint.h"
#include "../ddecal/gain_solvers/SolveData.h"
#include "../ddecal/gain_solvers/SolverTools.h"
#include "../ddecal/linear_solvers/LLSSolver.h"

#include "IDGPredict.h"
#include "MsColumnReader.h"
#include "Predict.h"
#include "SagecalPredict.h"

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
      itsSolIntCount(1),
      itsFirstSolutionIndex(0),
      itsNChan(itsSettings.n_channels),
      itsUVWFlagStep(parset, prefix, Step::MsType::kRegular),
      itsStoreSolutionInBuffer(parset.getBool(prefix + "storebuffer", false)),
      itsStatStream() {
  if (!itsSettings.stat_filename.empty()) {
    itsStatStream = std::make_unique<std::ofstream>(itsSettings.stat_filename);
  }

  // Initialize steps
  initializeColumnReaders(parset, prefix);
  initializeIDG(parset, prefix);
  initializePredictSteps(parset, prefix);
  initializeInitialSolutionsH5Parm(parset, prefix);
}

void DDECal::initializeInitialSolutionsH5Parm(
    const common::ParameterSet& parset, const std::string& prefix) {
  itsInitialSolutionsH5ParmName =
      parset.getString(prefix + "initialsolutions.h5parm", "");
  if (itsInitialSolutionsH5ParmName.empty()) {
    return;
  }

  const std::vector<std::string> default_solution_tables = {"amplitude000",
                                                            "phase000"};
  itsInitialSolutionsSolTab = parset.getStringVector(
      prefix + "initialsolutions.soltab", default_solution_tables);
  itsInitialSolutions = std::make_unique<schaapcommon::h5parm::H5Parm>(
      itsInitialSolutionsH5ParmName, itsInitialSolutionsSolTab);

  const std::string interpolation_type =
      parset.getString(prefix + "initialsolutions.interpolation", "nearest");
  if (interpolation_type == "nearest") {
    itsInterpolationType = JonesParameters::InterpolationType::NEAREST;
  } else if (interpolation_type == "linear") {
    itsInterpolationType = JonesParameters::InterpolationType::LINEAR;
  } else {
    throw std::runtime_error("Unsupported interpolation method: " +
                             interpolation_type);
  }

  std::string gain_type =
      parset.getString(prefix + "initialsolutions.gaintype", "");
  if (gain_type.empty()) {
    gain_type =
        itsInitialSolutions->GetSolTab(itsInitialSolutionsSolTab[0]).GetType();
  }
  itsGainType = JonesParameters::H5ParmTypeStringToGainType(gain_type);

  itsMissingAntennaBehavior =
      JonesParameters::StringToMissingAntennaBehavior(parset.getString(
          prefix + "initialsolutions.missingantennabehavior", "error"));
}

void DDECal::initializeColumnReaders(const common::ParameterSet& parset,
                                     const std::string& prefix) {
  for (const std::string& col : itsSettings.model_data_columns) {
    itsDirections.emplace_back(1, col);
    itsDirectionNames.emplace_back(prefix + col);
    itsSteps.push_back(std::make_shared<MsColumnReader>(parset, prefix, col));
    setModelNextSteps(*itsSteps.back(), col, parset, prefix);
  }
}

void DDECal::initializeModelReuse() {
  const std::vector<std::pair<std::string, std::string>> reused_directions =
      itsSettings.GetReusedDirections(getInfoIn().GetDirections());

  for (const auto& [name, name_without_prefix] : reused_directions) {
    // Keep using the original name with prefix in DPBuffers.
    itsDirectionNames.emplace_back(name);
    itsReusedDirectionNames.emplace_back(name);

    itsDirections.emplace_back(1, name_without_prefix);

    // Add a nullptr step for this direction, so there is still an entry in
    // itsSteps for each direction.
    itsSteps.emplace_back();
  }
}

void DDECal::initializeIDG(const common::ParameterSet& parset,
                           const std::string& prefix) {
  // TODO it would be nicer to get a new method in the DS9 reader to first get
  // names of directions and pass that to idgpredict. It will then read it
  // itself instead of DDECal having to do everything. It is better to do it all
  // in IDGPredict, so we can also make it
  if (itsSettings.idg_region_filename.empty() &&
      itsSettings.idg_image_filenames.empty()) {
    return;
  }

  const std::vector<FitsReader> readers =
      IDGPredict::GetReaders(itsSettings.idg_image_filenames);
  std::vector<Facet> facets =
      IDGPredict::GetFacets(itsSettings.idg_region_filename, readers.front());

  for (size_t i = 0; i < facets.size(); ++i) {
    std::string facet_dir_label = facets[i].DirectionLabel();
    std::string dir_name =
        facet_dir_label.empty() ? ("dir" + std::to_string(i)) : facet_dir_label;

    itsDirections.emplace_back(1, dir_name);
    itsDirectionNames.emplace_back(prefix + dir_name);

    itsSteps.push_back(std::make_shared<IDGPredict>(
        parset, prefix, readers, std::vector<Facet>{facets[i]}));
    setModelNextSteps(*itsSteps.back(), facet_dir_label, parset, prefix);
  }
}

void DDECal::initializePredictSteps(const common::ParameterSet& parset,
                                    const std::string& prefix) {
  std::vector<std::vector<std::string>> directions =
      model::MakeDirectionList(itsSettings.directions, itsSettings.source_db);

  for (std::vector<std::string>& direction : directions) {
    if (itsSettings.use_sagecal_predict) {
      itsSteps.push_back(
          std::make_shared<SagecalPredict>(parset, prefix, direction));
    } else {
      itsSteps.push_back(std::make_shared<Predict>(parset, prefix, direction));
    }
    setModelNextSteps(*itsSteps.back(), direction.front(), parset, prefix);
    itsDirectionNames.push_back(prefix + direction.front());
    itsDirections.push_back(std::move(direction));
  }
}

void DDECal::setModelNextSteps(Step& step, const std::string& direction,
                               const common::ParameterSet& parset,
                               const std::string& prefix) const {
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

  itsSolver = ddecal::CreateSolver(itsSettings, infoIn.antennaNames());

  initializeModelReuse();

  if (itsDirections.size() == 0) {
    throw std::runtime_error(
        "DDECal initialized with 0 directions: something is wrong with your "
        "parset or your sourcedb");
  }

  itsSettings.PrepareSolutionsPerDirection(itsDirections.size());

  itsUVWFlagStep.updateInfo(infoIn);

  if (itsRequestedSolInt == 0) {
    itsRequestedSolInt = getInfoOut().ntime();
  }

  // Update info for substeps and set other required parameters
  for (std::shared_ptr<ModelDataStep>& step : itsSteps) {
    if (!step) continue;

    step->setInfo(infoIn);

    if (auto s = std::dynamic_pointer_cast<Predict>(step)) {
      s->SetThreadData(&itsMeasuresMutex);
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(step)) {
      itsSolIntCount =
          std::max(itsSolIntCount,
                   s->GetBufferSize() / itsSteps.size() / itsRequestedSolInt);
      // We increment by one so the IDGPredict will not flush in its process
      s->SetBufferSize(itsRequestedSolInt * itsSolIntCount + 1);
    } else if (!std::dynamic_pointer_cast<MsColumnReader>(step) &&
               !std::dynamic_pointer_cast<SagecalPredict>(step)) {
      throw std::runtime_error("DDECal received an invalid first model step");
    }
  }

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

  // For each sub step chain, add a resultstep and get required fields.
  itsRequiredFields.resize(itsSteps.size());
  itsResultSteps.resize(itsSteps.size());
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    if (!itsSteps[dir]) continue;

    itsRequiredFields[dir] = base::GetChainRequiredFields(itsSteps[dir]);

    itsResultSteps[dir] =
        std::make_shared<MultiResultStep>(itsRequestedSolInt * itsSolIntCount);

    // Add the resultstep to the end of the model next steps
    std::shared_ptr<Step> step = itsSteps[dir];
    while (step->getNextStep()) {
      step = step->getNextStep();
    }
    step->setNextStep(itsResultSteps[dir]);
  }

  if (itsNChan == 0 || itsNChan > getInfoOut().nchan()) {
    itsNChan = getInfoOut().nchan();
  }

  // Create lists with used antenna indices, similarly to
  // DPInfo::RemoveUnusedAntennas.
  itsAntennas1.resize(getInfoOut().getAnt1().size());
  itsAntennas2.resize(getInfoOut().getAnt2().size());
  for (size_t i = 0; i < itsAntennas1.size(); ++i) {
    itsAntennas1[i] = getInfoOut().antennaMap()[getInfoOut().getAnt1()[i]];
    itsAntennas2[i] = getInfoOut().antennaMap()[getInfoOut().getAnt2()[i]];
  }

  // Fill antenna info in H5Parm, need to convert from casa types to std types
  // Fill in metadata for all antennas, also those that may be filtered out.
  std::vector<std::string> antennaNames(getInfoOut().antennaNames().size());
  std::vector<std::array<double, 3>> antennaPos(
      getInfoOut().antennaPos().size());
  for (unsigned int i = 0; i < getInfoOut().antennaNames().size(); ++i) {
    antennaNames[i] = getInfoOut().antennaNames()[i];
    casacore::Quantum<casacore::Vector<double>> pos =
        getInfoOut().antennaPos()[i].get("m");
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  itsSolutionWriter.AddAntennas(antennaNames, antennaPos);

  size_t nSolTimes =
      (getInfoOut().ntime() + itsRequestedSolInt - 1) / itsRequestedSolInt;
  size_t nChannelBlocks = getInfoOut().nchan() / itsNChan;
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
        (chBlock + 1) * getInfoOut().nchan() / nChannelBlocks;
    const size_t blockSize =
        itsChanBlockStart[chBlock + 1] - itsChanBlockStart[chBlock];
    const double* freqStart =
        getInfoOut().chanFreqs().data() + itsChanBlockStart[chBlock];
    const double meanFreq =
        std::accumulate(freqStart, freqStart + blockSize, 0.0) / blockSize;
    itsChanBlockFreqs[chBlock] = meanFreq;
  }

  itsWeightsPerAntenna.assign(
      itsChanBlockFreqs.size() * getInfoOut().antennaUsed().size(), 0.0);

  itsSourceDirections.reserve(itsSteps.size());
  const std::map<std::string, dp3::base::Direction>& directions =
      getInfoOut().GetDirections();
  for (unsigned int i = 0; i < itsSteps.size(); ++i) {
    const std::shared_ptr<ModelDataStep>& step = itsSteps[i];
    base::Direction direction = getInfoOut().phaseCenterDirection();

    if (step) {
      direction = step->GetFirstDirection();
    } else {
      auto search_result = directions.find(itsDirectionNames[i]);
      if (search_result != directions.end()) {
        direction = search_result->second;
      }
    }

    itsSourceDirections.push_back(direction);
  }

  // Save directions of model data for next steps.
  if (itsSettings.keep_model_data) {
    assert(itsSourceDirections.size() == itsDirectionNames.size());

    for (size_t i = 0; i < itsDirectionNames.size(); ++i) {
      GetWritableInfoOut().GetDirections()[itsDirectionNames[i]] =
          itsSourceDirections[i];
    }
  }

  // Prepare positions and names for the used antennas only.
  const std::vector<std::string> used_antenna_names =
      getInfoOut().GetUsedAntennaNames();
  std::vector<std::array<double, 3>> used_antenna_positions;
  used_antenna_positions.reserve(getInfoOut().antennaUsed().size());
  for (const int& ant : getInfoOut().antennaUsed()) {
    used_antenna_positions.push_back(antennaPos[ant]);
  }

  for (ddecal::SolverBase* solver : itsSolver->ConstraintSolvers()) {
    InitializeSolverConstraints(*solver, itsSettings, used_antenna_positions,
                                used_antenna_names, itsSourceDirections,
                                itsChanBlockFreqs);
  }

  size_t nSt = getInfoOut().antennaUsed().size();
  // Give renumbered antennas to solver
  itsSolver->Initialize(nSt, itsSettings.solutions_per_direction,
                        nChannelBlocks);

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
     << "  solution interval:   " << itsRequestedSolInt << '\n'
     << "  nchan:               " << itsNChan << '\n'
     << "  direction count:     " << itsDirections.size() << '\n'
     << "  directions:          " << itsDirections << '\n'
     << "  sols per direction:  " << itsSettings.solutions_per_direction
     << '\n';
  if (!itsInitialSolutionsH5ParmName.empty()) {
    os << "  initial sols H5Parm: " << itsInitialSolutionsH5ParmName << '\n'
       << "               soltab: ";
    for (std::string table : itsInitialSolutionsSolTab) {
      os << table << " ";
    }
    os << '\n'
       << "                 type: "
       << JonesParameters::GainTypeToHumanReadableString(itsGainType) << '\n'
       << "        interp method: "
       << (itsInterpolationType == JonesParameters::InterpolationType::NEAREST
               ? "nearest"
               : "linear")
       << '\n'
       << "          missing ant: "
       << JonesParameters::MissingAntennaBehaviorToString(
              itsMissingAntennaBehavior)
       << '\n';
  }
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
     << "  subtract model:      " << itsSettings.subtract << '\n'
     << "  keep model:          " << itsSettings.keep_model_data << '\n';
  for (size_t i = 0; i < itsSteps.size(); ++i) {
    std::shared_ptr<Step> step = itsSteps[i];
    if (step) {
      os << "Model steps for direction " << itsDirections[i][0] << '\n';
      do {
        step->show(os);
        step = step->getNextStep();
      } while (step);
    } else {
      os << "Direction " << itsDirections[i][0] << " reuses data from "
         << itsDirectionNames[i];
    }
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
    if (!step) continue;
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

  bool use_initial_solutions = !itsInitialSolutionsH5ParmName.empty();

  if (propagate_solutions) {
    // Initialize solutions with those of the previous step.
    itsSols[solution_index] = itsSols[solution_index - 1];
  } else if (use_initial_solutions) {
    // Initialize solutions with those from an existing H5Parm, stored in memory
    // with H5Cache
    const std::vector<double> solution_timestamp = {itsAvgTime};
    const std::vector<std::string> used_antenna_names =
        getInfoOut().GetUsedAntennaNames();
    schaapcommon::h5parm::SolTab first_soltab =
        itsInitialSolutions->GetSolTab(itsInitialSolutionsSolTab[0]);
    schaapcommon::h5parm::SolTab second_soltab;
    if (itsInitialSolutionsSolTab.size() == 2) {
      // Only attempt to read a 2nd soltab when it's requested
      second_soltab =
          itsInitialSolutions->GetSolTab(itsInitialSolutionsSolTab[1]);
    }
    std::vector<std::unique_ptr<schaapcommon::h5parm::JonesParameters>>
        jones_parameters_per_direction(itsDirections.size());
    for (size_t dir = 0; dir < itsDirections.size(); ++dir) {
      // Retrieve initial solutions for the patch that's closest to the
      // directions in the sourcedb, since the patch names and the number of
      // patches in the provided skymodel might not be the same as the
      // directions in the H5Parm.
      const std::string closest_patch = itsInitialSolutions->GetNearestSource(
          itsSourceDirections[dir].ra, itsSourceDirections[dir].dec);
      const hsize_t direction_index = first_soltab.GetDirIndex(closest_patch);
      jones_parameters_per_direction[dir] =
          std::make_unique<schaapcommon::h5parm::JonesParameters>(
              itsChanBlockFreqs, solution_timestamp, used_antenna_names,
              itsGainType, itsInterpolationType, direction_index, &first_soltab,
              &second_soltab, false, 0, itsMissingAntennaBehavior);
    }

    const size_t n_directions = itsDirections.size();
    const size_t n_subsolutions = itsSettings.GetNSolutions();
    const size_t n_antennas_used = getInfoOut().antennaUsed().size();
    const size_t n_polarization_parameters_per_solution =
        itsSolver->NSolutionPolarizations();
    const size_t n_values_per_channel_block =
        n_subsolutions * n_antennas_used *
        n_polarization_parameters_per_solution;

    for (size_t channel_block = 0; channel_block < itsChanBlockFreqs.size();
         ++channel_block) {
      itsSols[solution_index][channel_block].resize(n_values_per_channel_block);
      for (size_t antenna_index = 0; antenna_index < n_antennas_used;
           ++antenna_index) {
        size_t n_assigned_subsolutions = 0;
        for (size_t direction_index = 0; direction_index < n_directions;
             ++direction_index) {
          const casacore::Cube<std::complex<float>>& jones_parameters =
              jones_parameters_per_direction[direction_index]->GetParms();
          size_t n_subsolutions_per_direction =
              itsSettings.solutions_per_direction[direction_index];
          for (size_t direction_solution_index = 0;
               direction_solution_index < n_subsolutions_per_direction;
               ++direction_solution_index) {
            size_t direction_dependent_solution_index =
                n_assigned_subsolutions + direction_solution_index;
            for (size_t polarization_index = 0;
                 polarization_index < n_polarization_parameters_per_solution;
                 ++polarization_index) {
              const size_t flattened_index =
                  antenna_index * n_subsolutions *
                      n_polarization_parameters_per_solution +
                  direction_dependent_solution_index *
                      n_polarization_parameters_per_solution +
                  polarization_index;
              itsSols[solution_index][channel_block][flattened_index] =
                  jones_parameters(polarization_index, antenna_index,
                                   channel_block);
            }
          }
          n_assigned_subsolutions += n_subsolutions_per_direction;
        }
      }
    }
  } else {
    const size_t n_solutions = itsSettings.GetNSolutions();
    const size_t n_solution_values = n_solutions *
                                     getInfoOut().antennaUsed().size() *
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
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nChanBlocks = itsChanBlockFreqs.size();
  const size_t channel_begin = itsChanBlockStart[cbIndex];
  const size_t channel_end = itsChanBlockStart[cbIndex + 1];
  // Set the antenna-based weights to zero
  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t ant1 = getInfoOut().antennaMap()[getInfoOut().getAnt1()[bl]];
    size_t ant2 = getInfoOut().antennaMap()[getInfoOut().getAnt2()[bl]];
    itsWeightsPerAntenna[ant1 * nChanBlocks + cbIndex] = 0.0;
    itsWeightsPerAntenna[ant2 * nChanBlocks + cbIndex] = 0.0;
  }
  // Set the visibility flags to true. SolverTools::AssignAndWeight will write
  // zeroes to the weighted data if it is flagged.
  for (std::unique_ptr<DPBuffer>& buffer : itsInputBuffers[bufferIndex]) {
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
  for (size_t dir = 0; dir < itsDirections.size(); ++dir) {
    // For directions that reuse model data, the model data of the various
    // time steps should already be in the input buffers.
    if (!itsSteps[dir]) continue;

    // For directions that do not, the model data has been computed in the
    // substeps now and is waiting in the tailing MultiResultStep of each
    // substep chain. Move the model data from there to the input buffers.
    if (auto s = dynamic_cast<IDGPredict*>(itsSteps[dir].get())) {
      itsTimerPredict.start();
      s->flush();
      itsTimerPredict.stop();
    }
    for (size_t i = 0; i < itsResultSteps[dir]->size(); ++i) {
      // In the MultiResultStep, DPBuffers are stored flattened (1D vector) as
      // compared to the target shape of the input buffers (vector of vectors):
      const size_t sol_int = i / itsRequestedSolInt;
      const size_t timestep = i % itsRequestedSolInt;
      itsInputBuffers[sol_int][timestep]->MoveData(
          *itsResultSteps[dir]->get()[i], "", itsDirectionNames[dir]);
    }
  }

  std::vector<ddecal::SolverBase*> solvers = itsSolver->ConstraintSolvers();
  const size_t n_channel_blocks = itsChanBlockFreqs.size();
  const size_t n_antennas = getInfoOut().antennaUsed().size();

  // DDECal requires the unweighted model model when the model is subtracted
  // after calibration. Since the model data can be large, memory allocation
  // for this optional feature is done conditionally.
  const bool keep_model_data = itsSettings.only_predict ||
                               itsSettings.subtract ||
                               itsSettings.keep_model_data;

  for (size_t i = 0; i < itsInputBuffers.size(); ++i) {
    const size_t solution_index = itsFirstSolutionIndex + i;

    // Keep the original size of itsInputBuffers[i] in case it is temporaryly
    // padded with empty buffers.
    const size_t original_size = itsInputBuffers[i].size();

    // To learn how this condition could be met see AST-1589.
    if (solution_index >= itsSols.size()) {
      throw std::runtime_error(
          "The number of computed time slots is smaller than "
          "the actual number of time intervals. Most likely, "
          "this is a consequence of irregular time intervals.");
    }

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

      aocommon::Logger::Debug
          << "Initializing DDECal solver for current calibration interval.\n";

      // If the input buffers are not filled up to the requested solution
      // interval, fill them with padding buffers. This padding is used by the
      // solver algorithms, and they are removed once the solver execution in
      // completed, so they are not passed to the next steps.
      // The padding buffers are filled with the data of the first buffer, but
      // they have all the data invalidated by setting all the flags to true.
      while (itsInputBuffers[i].size() < itsRequestedSolInt) {
        std::unique_ptr<DPBuffer> padding_buffer = std::make_unique<DPBuffer>();
        const common::Fields fields =
            kDataField | kFlagsField | kWeightsField | kUvwField;
        padding_buffer->Copy(*itsInputBuffers[i][0], fields, true);
        padding_buffer->GetFlags().fill(true);
        itsInputBuffers[i].push_back(std::move(padding_buffer));
      }

      // The last solution interval can be smaller.
      std::vector<base::DPBuffer> weighted_buffers(itsInputBuffers[i].size());

      const bool linear_mode =
          itsSettings.solver_algorithm == ddecal::SolverAlgorithm::kLowRank;
      ddecal::AssignAndWeight(itsInputBuffers[i], itsDirectionNames,
                              weighted_buffers, keep_model_data, linear_mode);

      InitializeSolutions(i);

      itsTimerSolve.start();

      switch (itsSettings.solver_data_use) {
        case ddecal::SolverDataUse::kSingle: {
          const ddecal::UniSolveData solve_data(
              weighted_buffers, itsDirectionNames, n_channel_blocks, n_antennas,
              itsSettings.solutions_per_direction, itsAntennas1, itsAntennas2);
          weighted_buffers.clear();
          if (itsSettings.model_weighted_constraints) {
            itsSolver->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }
          aocommon::Logger::Debug << "Running DDECal single-visibility solver "
                                     "for current calibration interval.\n";

          solveResult = itsSolver->Solve(solve_data, itsSols[solution_index],
                                         itsAvgTime / itsRequestedSolInt,
                                         itsStatStream.get());
        } break;
        case ddecal::SolverDataUse::kDual: {
          const ddecal::DuoSolveData solve_data(
              weighted_buffers, itsDirectionNames, n_channel_blocks, n_antennas,
              itsSettings.solutions_per_direction, itsAntennas1, itsAntennas2);
          weighted_buffers.clear();
          if (itsSettings.model_weighted_constraints) {
            itsSolver->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }

          aocommon::Logger::Debug << "Running DDECal dual-visibility solver "
                                     "for current calibration interval.\n";

          solveResult = itsSolver->Solve(solve_data, itsSols[solution_index],
                                         itsAvgTime / itsRequestedSolInt,
                                         itsStatStream.get());
        } break;
        case ddecal::SolverDataUse::kFull: {
          const ddecal::FullSolveData solve_data(
              weighted_buffers, itsDirectionNames, n_channel_blocks, n_antennas,
              itsSettings.solutions_per_direction, itsAntennas1, itsAntennas2);
          weighted_buffers.clear();
          if (itsSettings.model_weighted_constraints) {
            itsSolver->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }

          aocommon::Logger::Debug
              << "Running DDECal solver for current calibration interval.\n";

          solveResult = itsSolver->Solve(solve_data, itsSols[solution_index],
                                         itsAvgTime / itsRequestedSolInt,
                                         itsStatStream.get());
        } break;
      }

      itsTimerSolve.stop();

      itsNIter[solution_index] = solveResult.iterations;
      itsNApproxIter[solution_index] = solveResult.constraint_iterations;
    }

    // Restoring itsInputBuffers[i] to its original size to remove any padding
    // that may have been added.
    itsInputBuffers[i].resize(original_size);

    if (itsSettings.only_predict) {
      SumModels(i);
    } else if (itsSettings.subtract || itsSettings.keep_model_data) {
      CorrectAndSubtractModels(i);
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

    // Store constraint solutions if any constraint has a non-empty result
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
      itsInputBuffers[i].front()->SetSolution(itsSols[solution_index]);
    }
  }

  itsTimer.stop();

  for (size_t sol_int = 0; sol_int < itsInputBuffers.size(); ++sol_int) {
    for (size_t timestep = 0; timestep < itsInputBuffers[sol_int].size();
         ++timestep) {
      if (itsOriginalFlags.size() > 0) {
        // Restore original flags, if the UVWFlagger or flagChannelBlock ran.
        itsInputBuffers[sol_int][timestep]->GetFlags() =
            xt::view(itsOriginalFlags, sol_int, timestep, xt::all(), xt::all(),
                     xt::all());
      }
      // Push data (possibly changed) to next step
      getNextStep()->process(std::move(itsInputBuffers[sol_int][timestep]));
    }
  }

  itsTimer.start();
}

bool DDECal::process(std::unique_ptr<DPBuffer> bufin) {
  itsTimer.start();

  // Check that all extra input data is there.
  // TODO(AST-1241): Handle these dependencies using Fields.
  for (const std::string& name : itsReusedDirectionNames) {
    if (!bufin->HasData(name)) {
      throw std::runtime_error("DDECal '" + itsSettings.name +
                               "' did not receive model data named '" + name +
                               "'.");
    }
    assert(bufin->GetData(name).shape() == bufin->GetData().shape());
  }

  // Create a new solution interval if needed
  if (itsInputBuffers.empty() ||
      itsInputBuffers.back().size() == itsRequestedSolInt) {
    itsInputBuffers.emplace_back();
  }

  itsInputBuffers.back().push_back(std::move(bufin));
  doPrepare();

  if (itsInputBuffers.size() == itsSolIntCount &&
      itsInputBuffers.back().size() == itsRequestedSolInt) {
    doSolve();

    // Clean up, prepare for next iteration
    itsFirstSolutionIndex += itsInputBuffers.size();
    itsAvgTime = 0;
    itsVisInInterval.assign(itsVisInInterval.size(),
                            std::pair<size_t, size_t>(0, 0));
    itsWeightsPerAntenna.assign(itsWeightsPerAntenna.size(), 0.0);

    for (std::shared_ptr<MultiResultStep>& result_step : itsResultSteps) {
      if (result_step) result_step->clear();
    }
    itsInputBuffers.clear();
  }

  itsTimer.stop();

  return false;
}

void DDECal::doPrepare() {
  // When UVW flagging is enabled, this input buffer is passed through the
  // UVWFlagger by moving it to itsUVWFlagStep and then extracting it again
  // from itsDataResultStep. itsInputBuffers then holds the updated buffer.
  std::unique_ptr<DPBuffer>& input_buffer = itsInputBuffers.back().back();

  if (itsOriginalFlags.size() > 0) {
    // Save the original flags, so DDECal can restore any flags changed by the
    // UVWFlagger/flagChannelBlock before passing the buffer to the next step.
    const size_t solution_interval_index = itsInputBuffers.size() - 1;
    const size_t step_in_solution_interval = itsInputBuffers.back().size() - 1;
    xt::view(itsOriginalFlags, solution_interval_index,
             step_in_solution_interval, xt::all(), xt::all(), xt::all()) =
        input_buffer->GetFlags();
  }

  if (!itsUVWFlagStep.isDegenerate()) {
    itsUVWFlagStep.process(std::move(input_buffer));
    input_buffer = itsDataResultStep->take();
  }

  itsTimerPredict.start();
  aocommon::Logger::Debug
      << "Acquiring one timestep of model data for DDECal.\n";
  // Enclose the recursive_for
  {
    aocommon::RecursiveFor recursive_for;
    recursive_for.Run(0, itsSteps.size(), [&](size_t direction) {
      if (itsSteps[direction]) {  // When reusing model data, there is no step.
        // Don't process column readers yet; they need to be run serially (see
        // further below)
        const bool is_column_reader =
            dynamic_cast<MsColumnReader*>(itsSteps[direction].get());
        if (!is_column_reader)
          itsSteps[direction]->process(std::make_unique<DPBuffer>(
              *input_buffer, itsRequiredFields[direction]));
      }
    });
  }
  // Call column readers serially, since CasaCore does not support reading
  // multiple columns in parallel.
  for (size_t direction = 0; direction != itsSteps.size(); ++direction) {
    if (itsSteps[direction] &&
        dynamic_cast<MsColumnReader*>(itsSteps[direction].get())) {
      itsSteps[direction]->process(std::make_unique<DPBuffer>(
          *input_buffer, itsRequiredFields[direction]));
    }
  }

  // Handle weights and flags
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nCh = getInfoOut().nchan();
  const size_t nCr = 4;

  size_t nchanblocks = itsChanBlockFreqs.size();

  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t chanblock = 0;
    size_t ant1 = getInfoOut().antennaMap()[getInfoOut().getAnt1()[bl]];
    size_t ant2 = getInfoOut().antennaMap()[getInfoOut().getAnt2()[bl]];
    for (size_t ch = 0; ch < nCh; ++ch) {
      if (ch == itsChanBlockStart[chanblock + 1]) {
        chanblock++;
      }
      for (size_t cr = 0; cr < nCr; ++cr) {
        // Add this weight to both involved antennas
        const double weight = input_buffer->GetWeights()(bl, ch, cr) *
                              !input_buffer->GetFlags()(bl, ch, cr);
        itsWeightsPerAntenna[ant1 * nchanblocks + chanblock] += weight;
        itsWeightsPerAntenna[ant2 * nchanblocks + chanblock] += weight;

        itsVisInInterval[chanblock].first +=
            (weight > 0);                      // unflagged nr of vis
        itsVisInInterval[chanblock].second++;  // total nr of vis
      }
    }
  }

  const double weightFactor =
      1. / (nCh * (getInfoOut().antennaUsed().size() - 1) * nCr *
            itsRequestedSolInt);
  for (double& weight : itsWeightsPerAntenna) {
    weight *= weightFactor;
  }

  itsTimerPredict.stop();

  itsAvgTime += itsAvgTime + input_buffer->GetTime();
}

void DDECal::WriteSolutions() {
  itsTimerWrite.start();

  // Create antenna info for H5Parm, used antennas only.
  const std::vector<std::string> used_antenna_names =
      getInfoOut().GetUsedAntennaNames();

  const std::string history = "CREATE by " + DP3Version::AsString() + "\n" +
                              "step " + itsSettings.name + " in parset: \n" +
                              itsSettings.parset_string;

  itsSolutionWriter.Write(
      itsSols, itsConstraintSols, getInfoOut().startTime(),
      getInfoOut().lastTime(), getInfoOut().timeInterval(), itsRequestedSolInt,
      itsSettings.solutions_per_direction, itsSettings.mode, used_antenna_names,
      itsSourceDirections, itsDirections, getInfoOut().chanFreqs(),
      itsChanBlockFreqs, history);

  itsTimerWrite.stop();
}

void DDECal::finish() {
  itsTimer.start();

  if (!itsInputBuffers.empty()) {
    doSolve();
  }

  if (!itsSettings.only_predict) WriteSolutions();

  itsInputBuffers.clear();
  itsTimer.stop();

  // Let the next steps finish.
  getNextStep()->finish();
}

namespace {
aocommon::MC2x2 ApplyFullJonesSolution(
    const std::vector<std::complex<double>>& solutions, size_t solution_index1,
    size_t solution_index2, const std::complex<float>* model_data) {
  const aocommon::MC2x2 solution1(&solutions[solution_index1 * 4]);
  const aocommon::MC2x2 solution2(&solutions[solution_index2 * 4]);
  return solution1.Multiply(aocommon::MC2x2(model_data))
      .MultiplyHerm(solution2);
}

aocommon::MC2x2 ApplyDiagonalSolution(
    const std::vector<std::complex<double>>& solutions, size_t solution_index1,
    size_t solution_index2, const std::complex<float>* model_data) {
  const aocommon::MC2x2Diag solution1(&solutions[solution_index1 * 2]);
  const aocommon::MC2x2Diag solution2(&solutions[solution_index2 * 2]);
  return solution1 * aocommon::MC2x2(model_data) * solution2.HermTranspose();
}

aocommon::MC2x2 ApplyScalarSolution(
    const std::vector<std::complex<double>>& solutions, size_t solution_index1,
    size_t solution_index2, const std::complex<float>* model_data) {
  const std::complex<double> solution_factor =
      solutions[solution_index1] * std::conj(solutions[solution_index2]);
  return aocommon::MC2x2(model_data) * solution_factor;
}
}  // namespace

void DDECal::ApplySolution(
    DPBuffer& buffer, size_t baseline, size_t channel,
    const std::vector<std::complex<double>>& solutions) const {
  const size_t antenna1 = getInfoOut().getAnt1()[baseline];
  const size_t antenna2 = getInfoOut().getAnt2()[baseline];
  const size_t n_directions = itsDirectionNames.size();

  aocommon::MC2x2 sum_over_directions = aocommon::MC2x2::Zero();

  for (size_t direction = 0; direction != n_directions; ++direction) {
    const size_t solution_index1 = antenna1 * n_directions + direction;
    const size_t solution_index2 = antenna2 * n_directions + direction;

    std::complex<float>* model_data =
        &buffer.GetData(itsDirectionNames[direction])(baseline, channel, 0);

    aocommon::MC2x2 corrected_model_data;
    if (itsSolver->NSolutionPolarizations() == 4) {
      corrected_model_data = ApplyFullJonesSolution(
          solutions, solution_index1, solution_index2, model_data);
    } else if (itsSolver->NSolutionPolarizations() == 2) {
      corrected_model_data = ApplyDiagonalSolution(solutions, solution_index1,
                                                   solution_index2, model_data);
    } else {
      assert(itsSolver->NSolutionPolarizations() == 1);
      corrected_model_data = ApplyScalarSolution(solutions, solution_index1,
                                                 solution_index2, model_data);
    }

    // Always update sum_over_directions since 'if (itsSettings.subtract)'
    // is probably more expensive than adding 8 aligned doubles.
    sum_over_directions += corrected_model_data;

    if (itsSettings.keep_model_data) {
      corrected_model_data.AssignTo(model_data);
    }
  }

  if (itsSettings.subtract) {
    for (size_t correlation = 0; correlation < 4; ++correlation) {
      buffer.GetData()(baseline, channel, correlation) -=
          sum_over_directions.Get(correlation);
    }
  }
}

void DDECal::CorrectAndSubtractModels(size_t buffer_index) {
  // The original unweighted data and model data are still in the solution
  // interval. Here we apply the solutions to all the model data directions and
  // subtract them from the data if the "subtract" setting is true.

  const size_t n_solutions = itsSettings.GetNSolutions();
  if (n_solutions != itsSettings.solutions_per_direction.size()) {
    throw std::runtime_error(
        "Model correction is not implemented for DDECal with direction "
        "dependent solution intervals");
  }

  const size_t solution_index = itsFirstSolutionIndex + buffer_index;
  assert(solution_index < itsSols.size());
  std::vector<std::vector<std::complex<double>>>& solutions =
      itsSols[solution_index];

  assert(buffer_index < itsInputBuffers.size());
  std::vector<std::unique_ptr<DPBuffer>>& solution_interval =
      itsInputBuffers[buffer_index];

  for (std::unique_ptr<DPBuffer>& data_buffer : solution_interval) {
    for (size_t bl = 0; bl < getInfoOut().nbaselines(); ++bl) {
      size_t chanblock = 0;

      for (size_t ch = 0; ch < getInfoOut().nchan(); ++ch) {
        if (ch == itsChanBlockStart[chanblock + 1]) {
          chanblock++;
        }

        ApplySolution(*data_buffer, bl, ch, solutions[chanblock]);
      }
    }
    if (!itsSettings.keep_model_data) {
      for (const std::string& name : itsDirectionNames) {
        data_buffer->RemoveData(name);
      }
    }
  }
}

void DDECal::SumModels(size_t buffer_index) {
  assert(buffer_index < itsInputBuffers.size());
  std::vector<std::unique_ptr<DPBuffer>>& solution_interval =
      itsInputBuffers[buffer_index];

  for (std::unique_ptr<DPBuffer>& data_buffer : solution_interval) {
    for (auto name_iterator = itsDirectionNames.begin();
         name_iterator != itsDirectionNames.end(); ++name_iterator) {
      if (itsDirectionNames.begin() == name_iterator) {
        data_buffer->GetData().assign(data_buffer->GetData(*name_iterator));
      } else {
        data_buffer->GetData() += data_buffer->GetData(*name_iterator);
      }
      if (!itsSettings.keep_model_data) {
        data_buffer->RemoveData(*name_iterator);
      }
    }
  }
}

}  // namespace steps
}  // namespace dp3
