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

#include <xtensor/views/xview.hpp>

#include <aocommon/logger.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

#include <schaapcommon/facets/facet.h>

#include "base/DP3.h"
#include "base/Version.h"

#include "model/SourceDBUtil.h"

#include "common/StreamUtil.h"

#include "ddecal/SolverFactory.h"
#include "ddecal/constraints/SmoothnessConstraint.h"
#include "ddecal/gain_solvers/SolveData.h"
#include "ddecal/gain_solvers/SolverTools.h"
#include "ddecal/linear_solvers/LLSSolver.h"

#include "IDGPredict.h"
#include "MsColumnReader.h"
#include "Predict.h"
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
#include "SagecalPredict.h"
#endif

using aocommon::FitsReader;

using schaapcommon::facets::Facet;
using schaapcommon::h5parm::GainType;
using schaapcommon::h5parm::JonesParameters;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::ddecal::LLSSolver;
using dp3::ddecal::LLSSolverType;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

DDECal::DDECal(const common::ParameterSet& parset, const std::string& prefix)
    : settings_(parset, prefix),
      average_time_(0),
      solutions_(),
      requested_solution_interval_(settings_.solution_interval),
      n_solution_intervals_(1),
      first_solution_index_(0),
      n_channels_(settings_.n_channels),
      uvw_flag_step_(parset, prefix, Step::MsType::kRegular),
      store_solution_in_buffer_(parset.getBool(prefix + "storebuffer", false)),
      statistics_stream_() {
  if (!settings_.stat_filename.empty()) {
    statistics_stream_ =
        std::make_unique<std::ofstream>(settings_.stat_filename);
  }

  if (!settings_.h5parm_name.empty()) {
    solution_writer_ =
        std::make_unique<ddecal::SolutionWriter>(settings_.h5parm_name);
  }

  // Initialize steps
  initializeColumnReaders(parset, prefix);
  initializeIDG(parset, prefix);
  initializePredictSteps(parset, prefix);
  initializeInitialSolutionsH5Parm(parset, prefix);
}

void DDECal::initializeInitialSolutionsH5Parm(
    const common::ParameterSet& parset, const std::string& prefix) {
  initial_solutions_h5_parm_name =
      parset.getString(prefix + "initialsolutions.h5parm", "");
  if (initial_solutions_h5_parm_name.empty()) {
    return;
  }

  const std::vector<std::string> default_solution_tables = {"amplitude000",
                                                            "phase000"};
  initial_solutions_table_ = parset.getStringVector(
      prefix + "initialsolutions.soltab", default_solution_tables);
  initial_solutions_ = std::make_unique<schaapcommon::h5parm::H5Parm>(
      initial_solutions_h5_parm_name, initial_solutions_table_);

  for (const std::string& soltab_name : initial_solutions_table_) {
    solution_tables_.push_back(initial_solutions_->GetSolTab(soltab_name));
  }

  // Check if H5Parm stores full-Jones solutions.
  initial_solutions_are_full_jones_ = false;
  if (solution_tables_[0].HasAxis("pol") &&
      initial_solutions_table_.size() == 2 &&
      initial_solutions_table_[0].find("amplitude") != std::string::npos &&
      initial_solutions_table_[1].find("phase") != std::string::npos) {
    initial_solutions_are_full_jones_ =
        solution_tables_[0].GetAxis("pol").size == 4;
  }

  if (initial_solutions_are_full_jones_) {
    gain_types_.resize(1);
    gain_types_[0] = GainType::kFullJones;
  } else {
    gain_types_.reserve(initial_solutions_table_.size());
    for (const std::string& soltab_name : initial_solutions_table_) {
      gain_types_.push_back(JonesParameters::H5ParmTypeStringToGainType(
          initial_solutions_->GetSolTab(soltab_name).GetType()));
    }
  }

  const std::string interpolation_type =
      parset.getString(prefix + "initialsolutions.interpolation", "nearest");
  if (interpolation_type == "nearest") {
    interpolation_type_ = JonesParameters::InterpolationType::NEAREST;
  } else if (interpolation_type == "linear") {
    interpolation_type_ = JonesParameters::InterpolationType::LINEAR;
  } else {
    throw std::runtime_error("Unsupported interpolation method: " +
                             interpolation_type);
  }

  missing_antenna_behavior_ =
      JonesParameters::StringToMissingAntennaBehavior(parset.getString(
          prefix + "initialsolutions.missingantennabehavior", "error"));
}

void DDECal::initializeColumnReaders(const common::ParameterSet& parset,
                                     const std::string& prefix) {
  for (const std::string& col : settings_.model_data_columns) {
    patches_per_direction_.emplace_back(1, col);
    direction_names_.emplace_back(prefix + col);
    steps_.push_back(std::make_shared<MsColumnReader>(parset, prefix,
                                                      MsType::kRegular, col));
    setModelNextSteps(*steps_.back(), col, parset, prefix);
  }
}

void DDECal::initializeModelReuse() {
  const std::vector<std::pair<std::string, std::string>> reused_directions =
      settings_.GetReusedDirections(getInfoIn().GetDirections());

  for (const auto& [name, name_without_prefix] : reused_directions) {
    // Keep using the original name with prefix in DPBuffers.
    direction_names_.emplace_back(name);
    reused_direction_names_.emplace_back(name);

    patches_per_direction_.emplace_back(1, name_without_prefix);

    // Add a nullptr step for this direction, so there is still an entry in
    // itsSteps for each direction.
    steps_.emplace_back();
  }
}

void DDECal::initializeIDG(const common::ParameterSet& parset,
                           const std::string& prefix) {
  // TODO it would be nicer to get a new method in the DS9 reader to first get
  // names of directions and pass that to idgpredict. It will then read it
  // itself instead of DDECal having to do everything. It is better to do it all
  // in IDGPredict, so we can also make it
  if (settings_.idg_region_filename.empty() &&
      settings_.idg_image_filenames.empty()) {
    return;
  }

  const std::vector<FitsReader> readers =
      IDGPredict::GetReaders(settings_.idg_image_filenames);
  std::vector<Facet> facets =
      IDGPredict::GetFacets(settings_.idg_region_filename, readers.front());

  for (size_t i = 0; i < facets.size(); ++i) {
    std::string facet_dir_label = facets[i].DirectionLabel();
    std::string dir_name =
        facet_dir_label.empty() ? ("dir" + std::to_string(i)) : facet_dir_label;

    patches_per_direction_.emplace_back(1, dir_name);
    direction_names_.emplace_back(prefix + dir_name);

    steps_.push_back(std::make_shared<IDGPredict>(
        parset, prefix, readers, std::vector<Facet>{facets[i]}));
    setModelNextSteps(*steps_.back(), facet_dir_label, parset, prefix);
  }
}

void DDECal::initializePredictSteps(const common::ParameterSet& parset,
                                    const std::string& prefix) {
  std::vector<std::vector<std::string>> directions =
      model::MakeDirectionList(settings_.directions, settings_.source_db);

  for (std::vector<std::string>& direction : directions) {
    if (settings_.use_sagecal_predict) {
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
      steps_.push_back(
          std::make_shared<SagecalPredict>(parset, prefix, direction));
#else
      throw std::runtime_error(
          "use_sagecal_predict setting is not supported since DP3 was built "
          "without "
          "SAGECal support.");
#endif
    } else {
      steps_.push_back(std::make_shared<Predict>(parset, prefix, direction));
    }
    setModelNextSteps(*steps_.back(), direction.front(), parset, prefix);
    direction_names_.push_back(prefix + direction.front());
    patches_per_direction_.push_back(std::move(direction));
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

  solver_ = ddecal::CreateSolver(settings_, infoIn.antennaNames());

  initializeModelReuse();

  if (patches_per_direction_.size() == 0) {
    throw std::runtime_error(
        "DDECal initialized with 0 directions: something is wrong with your "
        "parset or your sourcedb");
  }

  settings_.PrepareSubSolutionsPerDirection(patches_per_direction_.size());

  uvw_flag_step_.updateInfo(infoIn);

  if (requested_solution_interval_ == 0) {
    requested_solution_interval_ = getInfoOut().ntime();
  }

  // Update info for substeps and set other required parameters
  for (std::shared_ptr<ModelDataStep>& step : steps_) {
    if (!step) continue;

    step->setInfo(infoIn);

    if (auto s = std::dynamic_pointer_cast<Predict>(step)) {
      s->SetThreadData(&measures_mutex_);
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(step)) {
      n_solution_intervals_ =
          std::max(n_solution_intervals_, s->GetBufferSize() / steps_.size() /
                                              requested_solution_interval_);
      // We increment by one so the IDGPredict will not flush in its process
      s->SetBufferSize(requested_solution_interval_ * n_solution_intervals_ +
                       1);
    } else if (std::dynamic_pointer_cast<MsColumnReader>(step)) {
      // Step is valid. There's nothing to set.
#if defined(HAVE_LIBDIRAC) || defined(HAVE_LIBDIRAC_CUDA)
    } else if (std::dynamic_pointer_cast<SagecalPredict>(step)) {
      // Step is valid. There's nothing to set.
#endif
    } else {
      throw std::runtime_error("DDECal received an invalid first model step");
    }
  }

  if (!uvw_flag_step_.isDegenerate()) {
    data_result_step_ = std::make_shared<ResultStep>();
    uvw_flag_step_.setNextStep(data_result_step_);
  }
  if (!uvw_flag_step_.isDegenerate() || settings_.min_vis_ratio > 0.0) {
    // The UVWFlagger and/or flagChannelBlock() may modify the flags.
    // Create a buffer for storing the original flags.
    original_flags_.resize({n_solution_intervals_, requested_solution_interval_,
                            infoIn.nbaselines(), infoIn.nchan(),
                            infoIn.ncorr()});
  }

  // For each sub step chain, add a resultstep and get required fields.
  required_fields_.resize(steps_.size());
  result_steps_.resize(steps_.size());
  for (size_t dir = 0; dir < steps_.size(); ++dir) {
    if (!steps_[dir]) continue;

    required_fields_[dir] = base::GetChainRequiredFields(steps_[dir]);

    result_steps_[dir] = std::make_shared<MultiResultStep>(
        requested_solution_interval_ * n_solution_intervals_);

    // Add the resultstep to the end of the model next steps
    std::shared_ptr<Step> step = steps_[dir];
    while (step->getNextStep()) {
      step = step->getNextStep();
    }
    step->setNextStep(result_steps_[dir]);
  }

  if (n_channels_ == 0 || n_channels_ > getInfoOut().nchan()) {
    n_channels_ = getInfoOut().nchan();
  }

  // Create lists with used antenna indices, similarly to
  // DPInfo::RemoveUnusedAntennas.
  antennas1_.resize(getInfoOut().getAnt1().size());
  antennas2_.resize(getInfoOut().getAnt2().size());
  for (size_t i = 0; i < antennas1_.size(); ++i) {
    antennas1_[i] = getInfoOut().antennaMap()[getInfoOut().getAnt1()[i]];
    antennas2_[i] = getInfoOut().antennaMap()[getInfoOut().getAnt2()[i]];
  }

  // Convert antenna positions from casacore types to std types.
  std::vector<std::array<double, 3>> antennaPos(
      getInfoOut().antennaPos().size());
  for (unsigned int i = 0; i < getInfoOut().antennaPos().size(); ++i) {
    const casacore::Quantum<casacore::Vector<double>> pos =
        getInfoOut().antennaPos()[i].get("m");
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  if (solution_writer_) {
    // Fill antenna info in H5Parm.
    // Fill in metadata for all antennas, also those that may be filtered out.
    solution_writer_->AddAntennas(getInfoOut().antennaNames(), antennaPos);
  }

  size_t nSolTimes = (getInfoOut().ntime() + requested_solution_interval_ - 1) /
                     requested_solution_interval_;
  size_t nChannelBlocks = getInfoOut().nchan() / n_channels_;
  solutions_.resize(nSolTimes);
  n_iterations_.resize(nSolTimes);
  n_approximating_iterations_.resize(nSolTimes);
  constraint_solutions_.resize(nSolTimes);
  visibilities_in_interval_.assign(nChannelBlocks,
                                   std::pair<size_t, size_t>(0, 0));

  channel_block_start_.resize(nChannelBlocks + 1);
  channel_block_frequencies_.resize(nChannelBlocks);
  channel_block_start_.front() = 0;
  for (size_t chBlock = 0; chBlock != nChannelBlocks; ++chBlock) {
    channel_block_start_[chBlock + 1] =
        (chBlock + 1) * getInfoOut().nchan() / nChannelBlocks;
    const size_t blockSize =
        channel_block_start_[chBlock + 1] - channel_block_start_[chBlock];
    const double* freqStart =
        getInfoOut().chanFreqs().data() + channel_block_start_[chBlock];
    const double meanFreq =
        std::accumulate(freqStart, freqStart + blockSize, 0.0) / blockSize;
    channel_block_frequencies_[chBlock] = meanFreq;
  }

  weights_per_antenna_.assign(
      channel_block_frequencies_.size() * getInfoOut().antennaUsed().size(),
      0.0);

  source_directions_.reserve(steps_.size());
  std::map<std::string, dp3::base::Direction>& directions =
      GetWritableInfoOut().GetDirections();
  for (unsigned int i = 0; i < steps_.size(); ++i) {
    const std::shared_ptr<ModelDataStep>& step = steps_[i];
    base::Direction direction = getInfoOut().phaseCenterDirection();

    if (step) {
      direction = step->GetFirstDirection();
    } else {
      auto search_result = directions.find(direction_names_[i]);
      if (search_result != directions.end()) {
        direction = search_result->second;
        if (!settings_.keep_model_data) {
          // Remove the consumed directions from the output directions list.
          directions.erase(search_result);
        }
      }
    }

    source_directions_.push_back(direction);
  }

  // Save directions of model data for next steps.
  if (settings_.keep_model_data) {
    assert(source_directions_.size() == direction_names_.size());

    for (size_t i = 0; i < direction_names_.size(); ++i) {
      GetWritableInfoOut().GetDirections()[direction_names_[i]] =
          source_directions_[i];
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

  for (ddecal::SolverBase* solver : solver_->ConstraintSolvers()) {
    InitializeSolverConstraints(*solver, settings_, used_antenna_positions,
                                used_antenna_names, source_directions_,
                                channel_block_frequencies_);
  }

  size_t nSt = getInfoOut().antennaUsed().size();
  // Give renumbered antennas to solver
  solver_->Initialize(nSt, settings_.sub_solutions_per_direction,
                      nChannelBlocks);

  for (size_t i = 0; i < nSolTimes; ++i) {
    solutions_[i].resize(nChannelBlocks);
  }
}

void DDECal::show(std::ostream& os) const {
  os << "DDECal " << settings_.name << '\n'
     << "  mode (constraints):  " << ToString(settings_.mode) << '\n'
     << "  algorithm:           "
     << ddecal::ToString(settings_.solver_algorithm) << '\n'
     << "  H5Parm:              " << settings_.h5parm_name << '\n'
     << "  write sol to buffer: " << std::boolalpha << store_solution_in_buffer_
     << '\n'
     << "  solution interval:   " << requested_solution_interval_ << '\n'
     << "  nchan:               " << n_channels_ << '\n'
     << "  direction count:     " << patches_per_direction_.size() << '\n'
     << "  directions:          " << patches_per_direction_ << '\n'
     << "  sols per direction:  " << settings_.sub_solutions_per_direction
     << '\n';
  if (!initial_solutions_h5_parm_name.empty()) {
    os << "  initial sols H5Parm: " << initial_solutions_h5_parm_name << '\n'
       << "               soltab: ";
    for (size_t i = 0; i < initial_solutions_table_.size(); ++i) {
      os << initial_solutions_table_[i]
         << (i == initial_solutions_table_.size() - 1 ? " " : ", ");
    }
    os << '\n' << "                 type: ";
    for (size_t i = 0; i < gain_types_.size(); ++i) {
      os << JonesParameters::GainTypeToHumanReadableString(gain_types_[i])
         << (i == gain_types_.size() - 1 ? " " : ", ");
    }
    os << '\n'
       << "        interp method: "
       << (interpolation_type_ == JonesParameters::InterpolationType::NEAREST
               ? "nearest"
               : "linear")
       << '\n'
       << "          missing ant: "
       << JonesParameters::MissingAntennaBehaviorToString(
              missing_antenna_behavior_)
       << '\n';
  }
  if (settings_.min_vis_ratio != 0.0) {
    os << "  min visib. ratio:    " << settings_.min_vis_ratio << '\n';
  }
  os << "  tolerance:           " << solver_->GetAccuracy() << '\n'
     << "  max iter:            " << solver_->GetMaxIterations() << '\n'
     << "  flag unconverged:    " << std::boolalpha
     << settings_.flag_unconverged << '\n'
     << "     diverged only:    " << std::boolalpha
     << settings_.flag_diverged_only << '\n'
     << "  propagate solutions: " << std::boolalpha
     << settings_.propagate_solutions << '\n'
     << "       converged only: " << std::boolalpha
     << settings_.propagate_converged_only << '\n'
     << "  detect stalling:     " << std::boolalpha
     << solver_->GetDetectStalling() << '\n'
     << "  step size:           " << solver_->GetStepSize() << '\n';
  ShowConstraintSettings(os, settings_);
  os << "  approximate fitter:  " << settings_.approximate_tec << '\n'
     << "  only predict:        " << settings_.only_predict << '\n'
     << "  subtract model:      " << settings_.subtract << '\n'
     << "  keep model:          " << settings_.keep_model_data << '\n';
  for (size_t i = 0; i < steps_.size(); ++i) {
    std::shared_ptr<Step> step = steps_[i];
    if (step) {
      os << "Model steps for direction " << patches_per_direction_[i][0]
         << '\n';
      do {
        step->show(os);
        step = step->getNextStep();
      } while (step);
    } else {
      os << "Direction " << patches_per_direction_[i][0] << " reuses data from "
         << direction_names_[i];
    }
    os << '\n';
  }
  uvw_flag_step_.show(os);
}

void DDECal::showTimings(std::ostream& os, double duration) const {
  double totaltime = timer_.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " DDECal " << settings_.name << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, predict_timer_.getElapsed(), totaltime);
  os << " of it spent in predict" << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, solve_timer_.getElapsed(), totaltime);
  os << " of it spent in estimating gains and computing residuals" << '\n';

  solver_->GetTimings(os, solve_timer_.getElapsed());

  os << "          ";
  FlagCounter::showPerc1(os, write_timer_.getElapsed(), totaltime);
  os << " of it spent in writing gain solutions to disk" << '\n';

  os << "          ";
  os << "Substeps taken:" << '\n';
  for (const std::shared_ptr<ModelDataStep>& step : steps_) {
    if (!step) continue;

    std::ostringstream step_stream;
    step->showTimings(step_stream, duration);
    const std::string step_string = std::move(step_stream).str();
    if (!step_string.empty()) os << "          " << step_string;
  }

  os << "Iterations taken: [";
  for (size_t i = 0; i < n_iterations_.size(); ++i) {
    if (i > 0) os << ",";
    os << n_iterations_[i];
    if (n_approximating_iterations_[i] != 0)
      os << '|' << n_approximating_iterations_[i];
  }
  os << "]" << '\n';
}

xt::xtensor<std::complex<float>, 3> DDECal::ReadJonesMatrixFromH5Parm(
    const base::Direction& direction, double timestamp,
    schaapcommon::h5parm::GainType gain_type,
    schaapcommon::h5parm::SolTab* first_soltab,
    schaapcommon::h5parm::SolTab* second_soltab) {
  // Retrieve initial solutions for the patch that's closest to the
  // directions in the sourcedb, since the patch names and the number of
  // patches in the provided skymodel might not be the same as the
  // directions in the H5Parm.
  std::vector<double> timestamps = {timestamp};
  const std::string closest_patch =
      initial_solutions_->GetNearestSource(direction.ra, direction.dec);
  const hsize_t soltab_direction_index =
      first_soltab->GetDirIndex(closest_patch);
  size_t n_parm_values = 0;
  if (gain_type == GainType::kTec || gain_type == GainType::kClock) {
    n_parm_values =
        first_soltab->HasAxis("pol") ? first_soltab->GetAxis("pol").size : 1;
  }
  auto jones_parameters =
      std::make_unique<schaapcommon::h5parm::JonesParameters>(
          channel_block_frequencies_, timestamps,
          getInfoOut().GetUsedAntennaNames(), gain_type, interpolation_type_,
          soltab_direction_index, first_soltab, second_soltab, false,
          n_parm_values, missing_antenna_behavior_);

  // Write casacore Cube with Jones parameters to xtensor.
  const casacore::Cube<std::complex<float>>& jones_matrix =
      jones_parameters->GetParms();
  // Ensure that the xtensor has shape: frequency x antenna x polarization.
  const std::array<size_t, 3> xt_shape = {
      static_cast<size_t>(jones_matrix.shape()[2]),
      static_cast<size_t>(jones_matrix.shape()[1]),
      static_cast<size_t>(jones_matrix.shape()[0])};
  xt::xtensor<std::complex<float>, 3> jones_tensor =
      xt::adapt(jones_matrix.tovector(), xt_shape);

  return jones_tensor;
}

void DDECal::InitializeSolutions(size_t buffer_index) {
  const size_t solution_index = first_solution_index_ + buffer_index;
  assert(solution_index < solutions_.size());

  bool propagate_solutions =
      solution_index > 0 && settings_.propagate_solutions;
  if (propagate_solutions &&
      n_iterations_[solution_index - 1] > solver_->GetMaxIterations() &&
      settings_.propagate_converged_only) {
    propagate_solutions = false;
  }

  bool use_initial_solutions = !initial_solutions_h5_parm_name.empty();

  if (propagate_solutions) {
    // Initialize solutions with those of the previous step.
    solutions_[solution_index] = solutions_[solution_index - 1];
  } else if (use_initial_solutions) {
    // Initialize solutions with those from an existing H5Parm, stored in memory
    // with H5Cache
    const std::vector<std::string> used_antenna_names =
        getInfoOut().GetUsedAntennaNames();

    const size_t n_directions = patches_per_direction_.size();
    const size_t n_subsolutions = settings_.GetNSolutions();
    const size_t n_antennas_used = getInfoOut().antennaUsed().size();
    const size_t n_polarization_parameters_per_solution =
        solver_->NSolutionPolarizations();
    const size_t n_frequencies = channel_block_frequencies_.size();
    const size_t n_values_per_channel_block =
        n_subsolutions * n_antennas_used *
        n_polarization_parameters_per_solution;

    const std::array<size_t, 4> xt_shape = {
        n_directions, n_frequencies, n_antennas_used,
        n_polarization_parameters_per_solution};
    xt::xtensor<std::complex<float>, 4> jones_parameters_per_direction(
        xt_shape);
    if (initial_solutions_are_full_jones_) {
      for (size_t direction_index = 0; direction_index < n_directions;
           ++direction_index) {
        xt::view(jones_parameters_per_direction, direction_index, xt::all(),
                 xt::all(), xt::all()) =
            ReadJonesMatrixFromH5Parm(
                source_directions_[direction_index], average_time_,
                gain_types_[0], &solution_tables_[0], &solution_tables_[1]);
      }
    } else {
      const size_t n_soltabs = initial_solutions_table_.size();
      // Load solutions per soltab and multiply to get the full Jones matrix for
      // each direction.
      const std::array<size_t, 5> xt_shape = {
          n_soltabs, n_directions, n_frequencies, n_antennas_used,
          n_polarization_parameters_per_solution};
      xt::xtensor<std::complex<float>, 5> jones_parameters_per_soltab(xt_shape);
      for (size_t soltab_index = 0; soltab_index < n_soltabs; ++soltab_index) {
        for (size_t direction_index = 0; direction_index < n_directions;
             ++direction_index) {
          xt::view(jones_parameters_per_soltab, soltab_index, direction_index,
                   xt::all(), xt::all(), xt::all()) =
              ReadJonesMatrixFromH5Parm(
                  source_directions_[direction_index], average_time_,
                  gain_types_[soltab_index], &solution_tables_[soltab_index],
                  nullptr);
        }
      }

      jones_parameters_per_direction =
          xt::prod(jones_parameters_per_soltab, {0});
    }

    for (size_t channel_block = 0; channel_block < n_frequencies;
         ++channel_block) {
      solutions_[solution_index][channel_block].resize(
          n_values_per_channel_block);
      for (size_t antenna_index = 0; antenna_index < n_antennas_used;
           ++antenna_index) {
        size_t n_assigned_subsolutions = 0;
        for (size_t direction_index = 0; direction_index < n_directions;
             ++direction_index) {
          const xt::xtensor<std::complex<float>, 3>& jones_parameters =
              xt::view(jones_parameters_per_direction, direction_index,
                       xt::all(), xt::all(), xt::all());
          size_t n_subsolutions_per_direction =
              settings_.sub_solutions_per_direction[direction_index];
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
              solutions_[solution_index][channel_block][flattened_index] =
                  jones_parameters(channel_block, antenna_index,
                                   polarization_index);
            }
          }
          n_assigned_subsolutions += n_subsolutions_per_direction;
        }
      }
    }
  } else {
    const size_t n_solutions = settings_.GetNSolutions();
    const size_t n_solution_values = n_solutions *
                                     getInfoOut().antennaUsed().size() *
                                     solver_->NSolutionPolarizations();

    if (solver_->NSolutionPolarizations() == 4) {
      // Initialize solutions with unity matrix [1 0 ; 0 1].
      for (std::vector<casacore::DComplex>& solution_vector :
           solutions_[solution_index]) {
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
           solutions_[solution_index]) {
        solution_vector.assign(n_solution_values, 1.0);
      }
    }
  }
}

void DDECal::flagChannelBlock(size_t cbIndex, size_t bufferIndex) {
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nChanBlocks = channel_block_frequencies_.size();
  const size_t channel_begin = channel_block_start_[cbIndex];
  const size_t channel_end = channel_block_start_[cbIndex + 1];
  // Set the antenna-based weights to zero
  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t ant1 = getInfoOut().antennaMap()[getInfoOut().getAnt1()[bl]];
    size_t ant2 = getInfoOut().antennaMap()[getInfoOut().getAnt2()[bl]];
    weights_per_antenna_[ant1 * nChanBlocks + cbIndex] = 0.0;
    weights_per_antenna_[ant2 * nChanBlocks + cbIndex] = 0.0;
  }
  // Set the visibility flags to true. SolverTools::AssignAndWeight will write
  // zeroes to the weighted data if it is flagged.
  for (std::unique_ptr<DPBuffer>& buffer : input_buffers_[bufferIndex]) {
    xt::view(buffer->GetFlags(), xt::all(),
             xt::range(channel_begin, channel_end), xt::all()) = true;
  }
}

void DDECal::checkMinimumVisibilities(size_t bufferIndex) {
  for (size_t cb = 0; cb != channel_block_frequencies_.size(); ++cb) {
    double fraction = double(visibilities_in_interval_[cb].first) /
                      visibilities_in_interval_[cb].second;
    if (fraction < settings_.min_vis_ratio) flagChannelBlock(cb, bufferIndex);
  }
}

void DDECal::doSolve() {
  for (size_t dir = 0; dir < patches_per_direction_.size(); ++dir) {
    // For directions that reuse model data, the model data of the various
    // time steps should already be in the input buffers.
    if (!steps_[dir]) continue;

    // For directions that do not, the model data has been computed in the
    // substeps now and is waiting in the tailing MultiResultStep of each
    // substep chain. Move the model data from there to the input buffers.
    if (auto s = dynamic_cast<IDGPredict*>(steps_[dir].get())) {
      predict_timer_.start();
      s->flush();
      predict_timer_.stop();
    }
    for (size_t i = 0; i < result_steps_[dir]->size(); ++i) {
      // In the MultiResultStep, DPBuffers are stored flattened (1D vector) as
      // compared to the target shape of the input buffers (vector of vectors):
      const size_t sol_int = i / requested_solution_interval_;
      const size_t timestep = i % requested_solution_interval_;
      input_buffers_[sol_int][timestep]->MoveData(*result_steps_[dir]->get()[i],
                                                  "", direction_names_[dir]);
    }
  }

  std::vector<ddecal::SolverBase*> solvers = solver_->ConstraintSolvers();
  const size_t n_channel_blocks = channel_block_frequencies_.size();
  const size_t n_antennas = getInfoOut().antennaUsed().size();

  // DDECal requires the unweighted model model when the model is subtracted
  // after calibration. Since the model data can be large, memory allocation
  // for this optional feature is done conditionally.
  const bool keep_model_data =
      settings_.only_predict || settings_.subtract || settings_.keep_model_data;

  for (size_t i = 0; i < input_buffers_.size(); ++i) {
    const size_t solution_index = first_solution_index_ + i;

    // Keep the original size of itsInputBuffers[i] in case it is temporaryly
    // padded with empty buffers.
    const size_t original_size = input_buffers_[i].size();

    // To learn how this condition could be met see AST-1589.
    if (solution_index >= solutions_.size()) {
      throw std::runtime_error(
          "The number of computed time slots is smaller than "
          "the actual number of time intervals. Most likely, "
          "this is a consequence of irregular time intervals.");
    }

    ddecal::SolverBase::SolveResult solveResult;
    if (!settings_.only_predict) {
      if (settings_.min_vis_ratio > 0.0) {
        checkMinimumVisibilities(i);
      }

      for (ddecal::SolverBase* solver : solvers) {
        for (const std::unique_ptr<ddecal::Constraint>& constraint :
             solver->GetConstraints()) {
          constraint->SetWeights(weights_per_antenna_);
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
      while (input_buffers_[i].size() < requested_solution_interval_) {
        std::unique_ptr<DPBuffer> padding_buffer = std::make_unique<DPBuffer>();
        const common::Fields fields =
            kDataField | kFlagsField | kWeightsField | kUvwField;
        padding_buffer->Copy(*input_buffers_[i][0], fields, true);
        padding_buffer->GetFlags().fill(true);
        input_buffers_[i].push_back(std::move(padding_buffer));
      }

      // The last solution interval can be smaller.
      std::vector<base::DPBuffer> weighted_buffers(input_buffers_[i].size());

      const bool linear_mode =
          settings_.solver_algorithm == ddecal::SolverAlgorithm::kLowRank;
      ddecal::AssignAndWeight(input_buffers_[i], direction_names_,
                              weighted_buffers, keep_model_data, linear_mode);

      InitializeSolutions(i);

      solve_timer_.start();

      switch (settings_.solver_data_use) {
        case ddecal::SolverDataUse::kSingle: {
          const ddecal::UniSolveData solve_data(
              weighted_buffers, direction_names_, n_channel_blocks, n_antennas,
              settings_.sub_solutions_per_direction, antennas1_, antennas2_);
          weighted_buffers.clear();
          if (settings_.model_weighted_constraints) {
            solver_->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }
          aocommon::Logger::Debug << "Running DDECal single-visibility solver "
                                     "for current calibration interval.\n";

          solveResult =
              solver_->Solve(solve_data, solutions_[solution_index],
                             average_time_ / requested_solution_interval_,
                             statistics_stream_.get());
        } break;
        case ddecal::SolverDataUse::kDual: {
          const ddecal::DuoSolveData solve_data(
              weighted_buffers, direction_names_, n_channel_blocks, n_antennas,
              settings_.sub_solutions_per_direction, antennas1_, antennas2_);
          weighted_buffers.clear();
          if (settings_.model_weighted_constraints) {
            solver_->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }

          aocommon::Logger::Debug << "Running DDECal dual-visibility solver "
                                     "for current calibration interval.\n";

          solveResult =
              solver_->Solve(solve_data, solutions_[solution_index],
                             average_time_ / requested_solution_interval_,
                             statistics_stream_.get());
        } break;
        case ddecal::SolverDataUse::kFull: {
          const ddecal::FullSolveData solve_data(
              weighted_buffers, direction_names_, n_channel_blocks, n_antennas,
              settings_.sub_solutions_per_direction, antennas1_, antennas2_);
          weighted_buffers.clear();
          if (settings_.model_weighted_constraints) {
            solver_->SetDdConstraintWeights(solve_data.GetSolutionWeights());
          }

          aocommon::Logger::Debug
              << "Running DDECal solver for current calibration interval.\n";

          solveResult =
              solver_->Solve(solve_data, solutions_[solution_index],
                             average_time_ / requested_solution_interval_,
                             statistics_stream_.get());
        } break;
      }

      solve_timer_.stop();

      n_iterations_[solution_index] = solveResult.iterations;
      n_approximating_iterations_[solution_index] =
          solveResult.constraint_iterations;
    }

    // Restoring itsInputBuffers[i] to its original size to remove any padding
    // that may have been added.
    input_buffers_[i].resize(original_size);

    if (settings_.only_predict) {
      if (!settings_.keep_model_data) SumModels(i);
    } else if (settings_.subtract || settings_.keep_model_data) {
      CorrectAndSubtractModels(i);
    }

    // Check for nonconvergence and flag if desired. Unconverged solutions are
    // identified by the number of iterations being one more than the max
    // allowed number
    if (solveResult.iterations > solver_->GetMaxIterations() &&
        settings_.flag_unconverged) {
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          if (settings_.flag_diverged_only) {
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
      constraint_solutions_[solution_index] = solveResult.results;
    }

    // Store calibration solution for later calibration application steps.
    if (store_solution_in_buffer_) {
      input_buffers_[i].front()->SetSolution(solutions_[solution_index]);
    }
  }

  timer_.stop();

  for (size_t sol_int = 0; sol_int < input_buffers_.size(); ++sol_int) {
    for (size_t timestep = 0; timestep < input_buffers_[sol_int].size();
         ++timestep) {
      if (original_flags_.size() > 0) {
        // Restore original flags, if the UVWFlagger or flagChannelBlock ran.
        input_buffers_[sol_int][timestep]->GetFlags() =
            xt::view(original_flags_, sol_int, timestep, xt::all(), xt::all(),
                     xt::all());
      }
      // Push data (possibly changed) to next step
      getNextStep()->process(std::move(input_buffers_[sol_int][timestep]));
    }
  }

  timer_.start();
}

bool DDECal::process(std::unique_ptr<DPBuffer> bufin) {
  timer_.start();

  // Check that all extra input data is there.
  // TODO(AST-1241): Handle these dependencies using Fields.
  for (const std::string& name : reused_direction_names_) {
    if (!bufin->HasData(name)) {
      throw std::runtime_error("DDECal '" + settings_.name +
                               "' did not receive model data named '" + name +
                               "'.");
    }
    assert(bufin->GetData(name).shape() == bufin->GetData().shape());
  }

  // Create a new solution interval if needed
  if (input_buffers_.empty() ||
      input_buffers_.back().size() == requested_solution_interval_) {
    input_buffers_.emplace_back();
  }

  input_buffers_.back().push_back(std::move(bufin));
  doPrepare();

  if (input_buffers_.size() == n_solution_intervals_ &&
      input_buffers_.back().size() == requested_solution_interval_) {
    doSolve();

    // Clean up, prepare for next iteration
    first_solution_index_ += input_buffers_.size();
    average_time_ = 0;
    visibilities_in_interval_.assign(visibilities_in_interval_.size(),
                                     std::pair<size_t, size_t>(0, 0));
    weights_per_antenna_.assign(weights_per_antenna_.size(), 0.0);

    for (std::shared_ptr<MultiResultStep>& result_step : result_steps_) {
      if (result_step) result_step->clear();
    }
    input_buffers_.clear();
  }

  timer_.stop();

  return false;
}

void DDECal::doPrepare() {
  // When UVW flagging is enabled, this input buffer is passed through the
  // UVWFlagger by moving it to itsUVWFlagStep and then extracting it again
  // from itsDataResultStep. itsInputBuffers then holds the updated buffer.
  std::unique_ptr<DPBuffer>& input_buffer = input_buffers_.back().back();

  if (original_flags_.size() > 0) {
    // Save the original flags, so DDECal can restore any flags changed by the
    // UVWFlagger/flagChannelBlock before passing the buffer to the next step.
    const size_t solution_interval_index = input_buffers_.size() - 1;
    const size_t step_in_solution_interval = input_buffers_.back().size() - 1;
    xt::view(original_flags_, solution_interval_index,
             step_in_solution_interval, xt::all(), xt::all(), xt::all()) =
        input_buffer->GetFlags();
  }

  if (!uvw_flag_step_.isDegenerate()) {
    uvw_flag_step_.process(std::move(input_buffer));
    input_buffer = data_result_step_->take();
  }

  predict_timer_.start();
  aocommon::Logger::Debug
      << "Acquiring one timestep of model data for DDECal.\n";
  // Enclose the recursive_for
  {
    aocommon::RecursiveFor recursive_for;
    recursive_for.Run(0, steps_.size(), [&](size_t direction) {
      if (steps_[direction]) {  // When reusing model data, there is no step.
        // Don't process column readers yet; they need to be run serially (see
        // further below)
        const bool is_column_reader =
            dynamic_cast<MsColumnReader*>(steps_[direction].get());
        if (!is_column_reader)
          steps_[direction]->process(std::make_unique<DPBuffer>(
              *input_buffer, required_fields_[direction]));
      }
    });
  }
  // Call column readers serially, since CasaCore does not support reading
  // multiple columns in parallel.
  for (size_t direction = 0; direction != steps_.size(); ++direction) {
    if (steps_[direction] &&
        dynamic_cast<MsColumnReader*>(steps_[direction].get())) {
      steps_[direction]->process(std::make_unique<DPBuffer>(
          *input_buffer, required_fields_[direction]));
    }
  }

  // Handle weights and flags
  const size_t nBl = getInfoOut().nbaselines();
  const size_t nCh = getInfoOut().nchan();
  const size_t nCr = 4;

  size_t nchanblocks = channel_block_frequencies_.size();

  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t chanblock = 0;
    size_t ant1 = getInfoOut().antennaMap()[getInfoOut().getAnt1()[bl]];
    size_t ant2 = getInfoOut().antennaMap()[getInfoOut().getAnt2()[bl]];
    for (size_t ch = 0; ch < nCh; ++ch) {
      if (ch == channel_block_start_[chanblock + 1]) {
        chanblock++;
      }
      for (size_t cr = 0; cr < nCr; ++cr) {
        // Add this weight to both involved antennas
        const double weight = input_buffer->GetWeights()(bl, ch, cr) *
                              !input_buffer->GetFlags()(bl, ch, cr);
        weights_per_antenna_[ant1 * nchanblocks + chanblock] += weight;
        weights_per_antenna_[ant2 * nchanblocks + chanblock] += weight;

        visibilities_in_interval_[chanblock].first +=
            (weight > 0);                               // unflagged nr of vis
        visibilities_in_interval_[chanblock].second++;  // total nr of vis
      }
    }
  }

  const double weightFactor =
      1. / (nCh * (getInfoOut().antennaUsed().size() - 1) * nCr *
            requested_solution_interval_);
  for (double& weight : weights_per_antenna_) {
    weight *= weightFactor;
  }

  predict_timer_.stop();

  average_time_ += average_time_ + input_buffer->GetTime();
}

void DDECal::WriteSolutions() {
  write_timer_.start();

  const std::string history = "CREATE by " + DP3Version::AsString() + "\n" +
                              "step " + settings_.name + " in parset: \n" +
                              settings_.parset_string;

  solution_writer_->Write(
      solutions_, constraint_solutions_, getInfoOut().startTime(),
      getInfoOut().lastTime(), getInfoOut().timeInterval(),
      requested_solution_interval_, settings_.sub_solutions_per_direction,
      settings_.mode, getInfoOut().GetUsedAntennaNames(), source_directions_,
      patches_per_direction_, getInfoOut().chanFreqs(),
      channel_block_frequencies_, history);

  write_timer_.stop();
}

void DDECal::finish() {
  timer_.start();

  if (!input_buffers_.empty()) {
    doSolve();
  }

  if (!settings_.only_predict) WriteSolutions();

  input_buffers_.clear();
  timer_.stop();

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
  const size_t n_directions = direction_names_.size();

  aocommon::MC2x2 sum_over_directions = aocommon::MC2x2::Zero();

  for (size_t direction = 0; direction != n_directions; ++direction) {
    const size_t solution_index1 = antenna1 * n_directions + direction;
    const size_t solution_index2 = antenna2 * n_directions + direction;

    std::complex<float>* model_data =
        &buffer.GetData(direction_names_[direction])(baseline, channel, 0);

    aocommon::MC2x2 corrected_model_data;
    if (solver_->NSolutionPolarizations() == 4) {
      corrected_model_data = ApplyFullJonesSolution(
          solutions, solution_index1, solution_index2, model_data);
    } else if (solver_->NSolutionPolarizations() == 2) {
      corrected_model_data = ApplyDiagonalSolution(solutions, solution_index1,
                                                   solution_index2, model_data);
    } else {
      assert(solver_->NSolutionPolarizations() == 1);
      corrected_model_data = ApplyScalarSolution(solutions, solution_index1,
                                                 solution_index2, model_data);
    }

    // Always update sum_over_directions since 'if (itsSettings.subtract)'
    // is probably more expensive than adding 8 aligned doubles.
    sum_over_directions += corrected_model_data;

    if (settings_.keep_model_data) {
      corrected_model_data.AssignTo(model_data);
    }
  }

  if (settings_.subtract) {
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

  const size_t n_solutions = settings_.GetNSolutions();
  if (n_solutions != settings_.sub_solutions_per_direction.size()) {
    throw std::runtime_error(
        "Model correction is not implemented for DDECal with direction "
        "dependent solution intervals");
  }

  const size_t solution_index = first_solution_index_ + buffer_index;
  assert(solution_index < solutions_.size());
  std::vector<std::vector<std::complex<double>>>& solutions =
      solutions_[solution_index];

  assert(buffer_index < input_buffers_.size());
  std::vector<std::unique_ptr<DPBuffer>>& solution_interval =
      input_buffers_[buffer_index];

  for (std::unique_ptr<DPBuffer>& data_buffer : solution_interval) {
    for (size_t bl = 0; bl < getInfoOut().nbaselines(); ++bl) {
      size_t chanblock = 0;

      for (size_t ch = 0; ch < getInfoOut().nchan(); ++ch) {
        if (ch == channel_block_start_[chanblock + 1]) {
          chanblock++;
        }

        ApplySolution(*data_buffer, bl, ch, solutions[chanblock]);
      }
    }
    if (!settings_.keep_model_data) {
      for (const std::string& name : direction_names_) {
        data_buffer->RemoveData(name);
      }
    }
  }
}

void DDECal::SumModels(size_t buffer_index) {
  assert(buffer_index < input_buffers_.size());

  std::vector<std::unique_ptr<DPBuffer>>& solution_interval =
      input_buffers_[buffer_index];

  for (std::unique_ptr<DPBuffer>& data_buffer : solution_interval) {
    for (auto name_iterator = direction_names_.begin();
         name_iterator != direction_names_.end(); ++name_iterator) {
      if (direction_names_.begin() == name_iterator) {
        data_buffer->GetData().assign(data_buffer->GetData(*name_iterator));
      } else {
        data_buffer->GetData() += data_buffer->GetData(*name_iterator);
      }
      data_buffer->RemoveData(*name_iterator);
    }
  }
}

}  // namespace steps
}  // namespace dp3
