// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaDdeCal.h"

#include <algorithm>

#include <dp3/base/DP3.h>

#include "../base/SourceDBUtil.h"
#include "../common/StreamUtil.h"
#include "../ddecal/gain_solvers/SolveData.h"
#include "../ddecal/SolverFactory.h"

#include "BdaGroupPredict.h"
#include "Predict.h"
#include "Version.h"

using dp3::base::BDABuffer;
using dp3::base::DPInfo;
using dp3::ddecal::BdaSolverBuffer;
using dp3::ddecal::SolveData;
using dp3::steps::BdaGroupPredict;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

BdaDdeCal::BdaDdeCal(const common::ParameterSet& parset,
                     const std::string& prefix)
    : settings_(parset, prefix),
      solution_writer_(),
      steps_(),
      result_steps_(),
      patches_(),
      input_buffers_(),
      model_buffers_(),
      solver_buffer_(),
      solver_(),
      solution_interval_(0.0),
      chan_block_start_freqs_(),
      antennas1_(),
      antennas2_(),
      solutions_(),
      iterations_(),
      approx_iterations_(),
      constraint_solutions_(),
      timer_(),
      predict_timer_(),
      solve_timer_(),
      write_timer_() {
  if (settings_.keep_model_data) {
    throw std::runtime_error("BdaDdeCal does not support keeping model data");
  }

  uvw_flagger_step_ =
      std::make_unique<UVWFlagger>(parset, prefix, Step::MsType::kBda);
  uvw_flagger_result_step_ = std::make_shared<BDAResultStep>();
  uvw_flagger_step_->setNextStep(uvw_flagger_result_step_);
  InitializePredictSteps(parset, prefix);

  if (!settings_.only_predict) {
    solver_ = ddecal::CreateSolver(settings_, parset, prefix);
    solution_writer_ =
        std::make_unique<ddecal::SolutionWriter>(settings_.h5parm_name);
  }
}

void BdaDdeCal::InitializePredictSteps(const common::ParameterSet& parset,
                                       const std::string& prefix) {
  std::vector<std::vector<std::string>> directions =
      base::MakeDirectionList(settings_.directions, settings_.source_db);

  if (directions.empty()) {
    throw std::invalid_argument(
        "Invalid input parset: direction(s) or sourcedb must be specified");
  }

  const bool bda_group_predict = parset.getBool(prefix + "grouppredict", false);
  for (std::vector<std::string>& source_patterns : directions) {
    patches_.push_back(std::move(source_patterns));

    if (bda_group_predict) {
      steps_.push_back(
          std::make_shared<BdaGroupPredict>(parset, prefix, patches_.back()));
    } else {
      steps_.push_back(std::make_shared<Predict>(
          parset, prefix, patches_.back(), Step::MsType::kBda));
    }
    result_steps_.push_back(std::make_shared<BDAResultStep>());
    steps_.back()->setNextStep(result_steps_.back());
  }
}

common::Fields BdaDdeCal::getRequiredFields() const {
  common::Fields fields = uvw_flagger_step_->getRequiredFields();
  for (std::shared_ptr<Step> direction_first_step : steps_) {
    fields |= base::GetChainRequiredFields(direction_first_step);
  }
  return fields;
}

void BdaDdeCal::updateInfo(const DPInfo& _info) {
  Step::updateInfo(_info);
  uvw_flagger_step_->updateInfo(_info);

  // Update info for substeps
  for (unsigned int i = 0; i < patches_.size(); i++) {
    steps_[i]->setInfo(_info);
  }

  // Convert antenna positions from Casacore to STL, if necessary.
  std::vector<std::array<double, 3>> antenna_pos;
  if (!settings_.only_predict || solution_writer_) {
    antenna_pos.resize(info().antennaPos().size());
    for (size_t i = 0; i < info().antennaPos().size(); ++i) {
      casacore::Quantum<casacore::Vector<double>> pos =
          info().antennaPos()[i].get("m");
      antenna_pos[i][0] = pos.getValue()[0];
      antenna_pos[i][1] = pos.getValue()[1];
      antenna_pos[i][2] = pos.getValue()[2];
    }
  }

  if (!settings_.only_predict) {
    solution_interval_ = _info.timeInterval() * info().ntime();
    size_t n_solution_intervals = 1;
    if (settings_.solution_interval > 0) {
      solution_interval_ = _info.timeInterval() * settings_.solution_interval;
      n_solution_intervals =
          (info().ntime() + settings_.solution_interval - 1) /
          settings_.solution_interval;
    }

    double max_bda_interval = *std::max_element(info().ntimeAvgs().begin(),
                                                info().ntimeAvgs().end()) *
                              info().timeInterval();
    if (solution_interval_ < max_bda_interval) {
      throw std::invalid_argument(
          "Using BDA rows that are longer than the solution interval is not "
          "supported. Use less BDA time averaging or a larger solution "
          "interval.");
    }

    double max_chan_avg = 1.0;
    for (size_t i = 0; i < info().nbaselines(); i++) {
      double chan_avg =
          static_cast<double>(info().origNChan()) / info().chanWidths(i).size();
      if (chan_avg > max_chan_avg) {
        max_chan_avg = chan_avg;
      }
    }

    int channels_per_chan_block =
        (settings_.n_channels == 0) ? info().origNChan() : settings_.n_channels;
    if (max_chan_avg > channels_per_chan_block) {
      throw std::invalid_argument(
          "BDA frequency averaging scheme gives a maximum amount of channel "
          "averaged of " +
          std::to_string(max_chan_avg) +
          ", which is higher than the number of channels in each channel block "
          "set for BdaDdeCal of " +
          std::to_string(channels_per_chan_block) +
          ". Please adjust the BDA settings to have a lower frequency "
          "averaging factor.");
    }

    solver_buffer_ = std::make_unique<ddecal::BdaSolverBuffer>(
        patches_.size(), _info.startTime(), solution_interval_,
        _info.nbaselines());

    DetermineChannelBlocks();

    // Create lists with used antenna indices, similarly to
    // DPInfo::removeUnusedAnt.
    antennas1_.resize(info().getAnt1().size());
    antennas2_.resize(info().getAnt2().size());
    for (size_t i = 0; i < antennas1_.size(); ++i) {
      antennas1_[i] = info().antennaMap()[info().getAnt1()[i]];
      antennas2_[i] = info().antennaMap()[info().getAnt2()[i]];
    }

    std::vector<std::string> used_antenna_names;
    std::vector<std::array<double, 3>> used_antenna_positions;
    used_antenna_names.reserve(info().antennaUsed().size());
    used_antenna_positions.reserve(info().antennaUsed().size());
    for (const int& ant : info().antennaUsed()) {
      used_antenna_names.push_back(info().antennaNames()[ant]);
      used_antenna_positions.push_back(antenna_pos[ant]);
    }

    const std::vector<base::Direction> source_directions =
        GetSourceDirections();
    const std::vector<double> channel_block_frequencies =
        GetChannelBlockFrequencies();
    for (ddecal::SolverBase* solver : solver_->ConstraintSolvers()) {
      InitializeSolverConstraints(
          *solver, settings_, used_antenna_positions, used_antenna_names,
          std::vector<size_t>(source_directions.size(), 1), source_directions,
          channel_block_frequencies);
    }

    solver_->Initialize(info().antennaUsed().size(),
                        std::vector<size_t>(patches_.size(), 1),
                        chan_block_start_freqs_.size() - 1);

    // SolveCurrentInterval will add solution intervals to the solutions.
    solutions_.reserve(n_solution_intervals);
    constraint_solutions_.reserve(n_solution_intervals);
  }

  if (solution_writer_) {
    // Convert Casacore types to STL types and pass antenna info to the
    // SolutionWriter.
    const std::vector<std::string> antenna_names(info().antennaNames().begin(),
                                                 info().antennaNames().end());
    solution_writer_->AddAntennas(antenna_names, antenna_pos);
  }
}

void BdaDdeCal::DetermineChannelBlocks() {
  // Combine channels into channel blocks. Since many baselines are not
  // averaged, use the same strategy as the regular DdeCal.
  size_t n_channel_blocks = 1;
  if (settings_.n_channels > 0) {
    n_channel_blocks =
        std::max(info().nchan() / settings_.n_channels, size_t(1));
  }

  // Although info().chanWidths(baseline) differs between baselines,
  // min_freq and max_freq should be equal for all baselines.
  const double min_freq =
      info().chanFreqs(0).front() - (info().chanWidths(0).front() / 2);
  const double max_freq =
      info().chanFreqs(0).back() + (info().chanWidths(0).back() / 2);
  const double chan_width = (max_freq - min_freq) / info().nchan();

  chan_block_start_freqs_.clear();
  chan_block_start_freqs_.reserve(n_channel_blocks + 1);
  chan_block_start_freqs_.push_back(min_freq);
  size_t start_index = 0;
  double start_freq = min_freq;

  for (size_t ch_block = 0; ch_block < n_channel_blocks; ++ch_block) {
    const size_t next_index =
        (ch_block + 1) * info().nchan() / n_channel_blocks;
    const size_t block_size = next_index - start_index;
    const double next_freq = start_freq + block_size * chan_width;
    chan_block_start_freqs_.push_back(next_freq);
    start_index = next_index;
    start_freq = next_freq;
  }
}

std::vector<base::Direction> BdaDdeCal::GetSourceDirections() const {
  std::vector<base::Direction> source_directions;
  source_directions.reserve(steps_.size());
  for (const std::shared_ptr<ModelDataStep>& s : steps_) {
    source_directions.push_back(s->GetFirstDirection());
  }
  return source_directions;
}

std::vector<double> BdaDdeCal::GetChannelBlockFrequencies() const {
  std::vector<double> frequencies;
  if (chan_block_start_freqs_.empty()) {
    assert(false);
  } else {
    frequencies.reserve(chan_block_start_freqs_.size() - 1);
    for (size_t i = 0; i < chan_block_start_freqs_.size() - 1; ++i) {
      frequencies.push_back(
          (chan_block_start_freqs_[i] + chan_block_start_freqs_[i + 1]) * 0.5);
    }
  }
  return frequencies;
}

bool BdaDdeCal::process(std::unique_ptr<base::BDABuffer> buffer) {
  timer_.start();

  // Feed metadata-only copies of the buffer to the steps. Only allocate room
  // for 'data' in the buffers, since the steps should only produce data.
  BDABuffer::Fields fields(false);
  fields.data = true;
  BDABuffer::Fields copyfields(false);

  if (!uvw_flagger_step_->isDegenerate()) {
    uvw_flagger_step_->process(std::move(buffer));
    std::vector<std::unique_ptr<base::BDABuffer>> uvw_flagged_buffer =
        uvw_flagger_result_step_->Extract();
    assert(1 == uvw_flagged_buffer.size());
    buffer = std::move(uvw_flagged_buffer.front());
  }

  predict_timer_.start();
  for (std::shared_ptr<ModelDataStep>& step : steps_) {
    step->process(std::make_unique<BDABuffer>(*buffer, fields, copyfields));
  }
  predict_timer_.stop();

  // Always store the input buffer, since BdaDdeCal should forward the
  // unmodified weights and flags to its next step.

  // When only_predict is false and all predict sub-steps have completed a model
  // buffer, solver_buffer_ will receive the input buffer and all model buffers.
  input_buffers_.push_back(std::move(buffer));

  ExtractResults();

  timer_.stop();
  ProcessCompleteDirections();

  return true;
}

void BdaDdeCal::ExtractResults() {
  // The BDABuffers from the sub-steps should have the same shape, however, a
  // step may delay outputting steps, e.g., due to internal buffering.
  for (size_t direction = 0; direction < result_steps_.size(); direction++) {
    std::vector<std::unique_ptr<BDABuffer>> results =
        result_steps_[direction]->Extract();
    if (!results.empty()) {
      // Find the queue index of the first free slot for this direction.
      size_t queue_index = 0;
      while (model_buffers_.size() > queue_index &&
             model_buffers_[queue_index][direction]) {
        ++queue_index;
      }

      for (std::unique_ptr<BDABuffer>& result : results) {
        // Extend the queue if it's not big enough.
        if (queue_index == model_buffers_.size()) {
          model_buffers_.emplace_back(result_steps_.size());
        }
        model_buffers_[queue_index][direction] = std::move(result);
        ++queue_index;
      }
    }
  }
}

void BdaDdeCal::ProcessCompleteDirections() {
  const auto pointer_is_set = [](const std::unique_ptr<BDABuffer>& pointer) {
    return !!pointer;  // Convert the pointer to a boolean.
  };
  while (!model_buffers_.empty() &&
         std::all_of(model_buffers_.front().begin(),
                     model_buffers_.front().end(), pointer_is_set)) {
    assert(!input_buffers_.empty() && input_buffers_.front());
    std::vector<std::unique_ptr<base::BDABuffer>>& direction_buffers =
        model_buffers_.front();
    if (settings_.only_predict) {
      // Sum all model buffers into the saved data buffer.
      std::complex<float> restrict* data = input_buffers_.front()->GetData();
      const size_t data_size = input_buffers_.front()->GetNumberOfElements();

      assert(
          direction_buffers.front()->IsMetadataEqual(*input_buffers_.front()));
      std::copy_n(direction_buffers.front()->GetData(), data_size, data);
      direction_buffers.front().reset();

      for (size_t dir = 1; dir < direction_buffers.size(); ++dir) {
        assert(direction_buffers[dir]->GetNumberOfElements() == data_size);
        const std::complex<float> restrict* other_data =
            direction_buffers[dir]->GetData();
        for (size_t j = 0; j < data_size; ++j) data[j] += other_data[j];
        direction_buffers[dir].reset();
      }

      getNextStep()->process(std::move(input_buffers_.front()));
    } else {
      // Send data buffer and model_buffers to solver_buffer_.
      solver_buffer_->AppendAndWeight(std::move(input_buffers_.front()),
                                      std::move(direction_buffers));
    }
    input_buffers_.pop_front();
    model_buffers_.pop_front();
  }

  if (!settings_.only_predict) {
    while (solver_buffer_->IntervalIsComplete()) {
      SolveCurrentInterval();
      solver_buffer_->AdvanceInterval();
    }

    for (std::unique_ptr<BDABuffer>& done_buffer : solver_buffer_->GetDone()) {
      getNextStep()->process(std::move(done_buffer));
    }
  }
}

void BdaDdeCal::SolveCurrentInterval() {
  timer_.start();
  solve_timer_.start();
  const size_t n_channel_blocks = chan_block_start_freqs_.size() - 1;
  const size_t n_antennas = info().antennaUsed().size();

  const bool linear_mode =
      settings_.solver_algorithm == ddecal::SolverAlgorithm::kLowRank;
  dp3::ddecal::SolveData data(*solver_buffer_, n_channel_blocks,
                              patches_.size(), n_antennas, antennas1_,
                              antennas2_, linear_mode);

  const int current_interval = solutions_.size();
  assert(current_interval == solver_buffer_->GetCurrentInterval());
  const double current_center =
      info().startTime() + (current_interval + 0.5) * solution_interval_;

  solutions_.emplace_back(n_channel_blocks);
  const size_t block_solution_size =
      patches_.size() * n_antennas * solver_->NSolutionPolarizations();
  for (std::vector<std::complex<double>>& block_solution : solutions_.back()) {
    block_solution.assign(block_solution_size, 1.0);
  }

  // The variables below hold the count of the number of time intervals present
  // in current solution interval per each antenna: this is needed for a correct
  // normalization, as different antennas may have different averaging schemes.
  std::vector<double> antenna1_interval_count(n_antennas, 0);
  std::vector<double> antenna2_interval_count(n_antennas, 0);

  // Get antenna weights
  std::vector<double> weights_per_antenna(n_channel_blocks * n_antennas, 0.0);
  for (const BDABuffer::Row* data_row :
       solver_buffer_->GetUnweightedDataRows()) {
    const size_t antenna1 = antennas1_[data_row->baseline_nr];
    const size_t antenna2 = antennas2_[data_row->baseline_nr];
    ++antenna1_interval_count[antenna1];
    ++antenna2_interval_count[antenna2];

    if (antenna1 != antenna2) {
      for (size_t ch = 0; ch < data_row->n_channels; ++ch) {
        const size_t index = ch * data_row->n_correlations;
        const bool* const flag_ptr = data_row->flags + index;
        const float* const weights_ptr = data_row->weights + index;

        const double normalization_factor =
            1.0 /
            (data_row->n_correlations * data_row->n_channels * n_antennas);

        // Add the weights for current antenna and channel block
        const size_t channel_block_index =
            GetChanBlockIndex(ch, data_row->n_channels, n_channel_blocks);

        for (size_t corr = 0; corr < data_row->n_correlations; ++corr) {
          if (!(flag_ptr[corr])) {
            weights_per_antenna[antenna1 * n_channel_blocks +
                                channel_block_index] +=
                weights_ptr[corr] * normalization_factor;
            weights_per_antenna[antenna2 * n_channel_blocks +
                                channel_block_index] +=
                weights_ptr[corr] * normalization_factor;
          }
        }
      }
    }
  }

  // Normalize the weigths to the number of intervals contained in current
  // solution interval, per antenna
  for (size_t k = 0; k < antenna1_interval_count.size(); ++k) {
    for (size_t i = 0; i < n_channel_blocks; ++i) {
      weights_per_antenna[k * n_channel_blocks + i] /=
          (antenna1_interval_count[k] + antenna2_interval_count[k]);
    }
  }

  std::vector<ddecal::SolverBase*> solvers = solver_->ConstraintSolvers();
  for (ddecal::SolverBase* solver : solvers) {
    for (const std::unique_ptr<ddecal::Constraint>& constraint :
         solver->GetConstraints()) {
      constraint->SetWeights(weights_per_antenna);
    }
  }

  InitializeCurrentSolutions();

  ddecal::SolverBase::SolveResult result =
      solver_->Solve(data, solutions_.back(), current_center, nullptr);

  assert(iterations_.size() == solutions_.size() - 1);
  assert(approx_iterations_.size() == solutions_.size() - 1);
  iterations_.push_back(result.iterations);
  approx_iterations_.push_back(result.constraint_iterations);

  if (settings_.subtract) {
    solver_buffer_->SubtractCorrectedModel(
        solutions_.back(), chan_block_start_freqs_,
        solver_->NSolutionPolarizations(), antennas1_, antennas2_,
        info().BdaChanFreqs());
  }

  // Check for nonconvergence and flag if desired. Unconverged solutions are
  // identified by the number of iterations being one more than the max
  // allowed number.
  if (settings_.flag_unconverged &&
      result.iterations > solver_->GetMaxIterations()) {
    if (settings_.flag_diverged_only) {
      // Set negative weights (indicating unconverged solutions that diverged)
      // to zero. All other unconverged solutions remain unflagged.
      for (auto& constraint_results : result.results) {
        for (auto& constraint_result : constraint_results) {
          for (double& weight : constraint_result.weights) {
            weight = std::max(weight, 0.0);
          }
        }
      }
    } else {
      // Set all weights to zero
      for (auto& constraint_results : result.results) {
        for (auto& constraint_result : constraint_results) {
          std::fill(constraint_result.weights.begin(),
                    constraint_result.weights.end(), 0.0);
        }
      }
    }
  } else {
    // Set negative weights (indicating unconverged solutions that diverged)
    // to one. All other unconverged solutions are unflagged already.
    for (auto& constraint_results : result.results) {
      for (auto& constraint_result : constraint_results) {
        for (double& weight : constraint_result.weights) {
          if (weight < 0.0) weight = 1.0;
        }
      }
    }
  }

  // Store constraint solutions if any constaint has a non-empty result.
  if (std::any_of(
          result.results.begin(), result.results.end(),
          [](const std::vector<ddecal::Constraint::Result>&
                 constraint_results) { return !constraint_results.empty(); })) {
    constraint_solutions_.push_back(std::move(result.results));
  } else {  // Add an empty constraint solution for the solution interval.
    constraint_solutions_.emplace_back();
  }
  assert(solutions_.size() == constraint_solutions_.size());

  solve_timer_.stop();
  timer_.stop();
}

void BdaDdeCal::InitializeCurrentSolutions() {
  std::vector<std::vector<std::complex<double>>>& current_solution =
      solutions_.back();

  bool propagate = solutions_.size() > 1 && settings_.propagate_solutions;
  if (propagate && settings_.propagate_converged_only &&
      iterations_[solutions_.size() - 2] > solver_->GetMaxIterations()) {
    propagate = false;
  }

  if (propagate) {
    // Copy the data from the previous solution, without reallocating
    // data structures.
    const std::vector<std::vector<std::complex<double>>>& prev_solution =
        *(solutions_.end() - 2);
    assert(current_solution.size() == prev_solution.size());

    for (size_t block = 0; block < current_solution.size(); ++block) {
      auto& current_block = current_solution[block];
      const auto& prev_block = prev_solution[block];
      assert(current_block.size() == prev_block.size());
      std::copy(prev_block.begin(), prev_block.end(), current_block.begin());
    }
  } else if (solver_->NSolutionPolarizations() == 4) {
    // Initialize full jones solutions with unity matrix.
    for (auto& current_block : current_solution) {
      for (size_t i = 0; i < current_block.size(); i += 4) {
        current_block[i + 0] = 1.0;
        current_block[i + 1] = 0.0;
        current_block[i + 2] = 0.0;
        current_block[i + 3] = 1.0;
      }
    }
  } else {
    // Initial scalar and diagonal solutions with 1.0.
    for (auto& block_solution : current_solution) {
      std::fill(block_solution.begin(), block_solution.end(), 1.0);
    }
  }
}

void BdaDdeCal::finish() {
  timer_.start();
  predict_timer_.start();
  for (std::shared_ptr<ModelDataStep>& step : steps_) {
    step->finish();
  }
  predict_timer_.stop();
  timer_.stop();

  ExtractResults();
  ProcessCompleteDirections();
  assert(input_buffers_.empty());

  if (!settings_.only_predict) {
    while (solver_buffer_->BufferCount() > 0) {
      SolveCurrentInterval();
      solver_buffer_->AdvanceInterval();
    }
    for (std::unique_ptr<BDABuffer>& done_buffer : solver_buffer_->GetDone()) {
      getNextStep()->process(std::move(done_buffer));
    }
    if (solution_writer_) WriteSolutions();
  }

  getNextStep()->finish();
}

void BdaDdeCal::WriteSolutions() {
  timer_.start();
  write_timer_.start();
  // Create antenna info for H5Parm, used antennas only.
  std::vector<std::string> used_antenna_names;
  used_antenna_names.reserve(info().antennaUsed().size());
  for (size_t used_antenna : info().antennaUsed()) {
    used_antenna_names.emplace_back(info().antennaNames()[used_antenna]);
  }

  const std::string history = "CREATE by " + DP3Version::AsString() + "\n" +
                              "step " + settings_.name + " in parset: \n" +
                              settings_.parset_string;

  solution_writer_->Write(
      solutions_, constraint_solutions_, info().startTime(), solution_interval_,
      settings_.mode, used_antenna_names, GetSourceDirections(), patches_,
      info().chanFreqs(), GetChannelBlockFrequencies(), history);

  write_timer_.stop();
  timer_.stop();
}

void BdaDdeCal::show(std::ostream& stream) const {
  stream << "BdaDdeCal " << settings_.name << '\n'
         << "  mode (constraints):  " << ToString(settings_.mode) << '\n'
         << "  directions:          " << patches_ << '\n';
  if (solver_) {
    const size_t nchan = settings_.n_channels == 0 ? size_t(getInfo().nchan())
                                                   : settings_.n_channels;
    stream << "  solver algorithm:    "
           << ddecal::ToString(settings_.solver_algorithm) << '\n'
           << "  H5Parm:              " << settings_.h5parm_name << '\n'
           << "  subtract model:      " << std::boolalpha << settings_.subtract
           << '\n'
           << "  solution interval:   " << solution_interval_ << " s\n"
           << "  #channels/block:     " << nchan << '\n'
           << "  #channel blocks:     " << chan_block_start_freqs_.size() - 1
           << '\n'
           << "  tolerance:           " << solver_->GetAccuracy() << '\n'
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
    ShowConstraintSettings(stream, settings_);
  }

  for (size_t dir = 0; dir < steps_.size(); ++dir) {
    stream << "Model steps for direction " << patches_[dir] << '\n';
    for (std::shared_ptr<Step> step = steps_[dir]; step;
         step = step->getNextStep()) {
      step->show(stream);
    }
    stream << '\n';
  }

  uvw_flagger_step_->show(stream);
}

void BdaDdeCal::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " BdaDdeCal \n";

  const double total_time = timer_.getElapsed();
  os << "          ";
  base::FlagCounter::showPerc1(os, predict_timer_.getElapsed(), total_time);
  os << " of it spent in predict\n";

  if (!settings_.only_predict) {
    os << "          ";
    base::FlagCounter::showPerc1(os, solve_timer_.getElapsed(), total_time);
    os << " of it spent in estimating gains and computing residuals\n";

    solver_->GetTimings(os, solve_timer_.getElapsed());

    os << "          ";
    base::FlagCounter::showPerc1(os, write_timer_.getElapsed(), total_time);
    os << " of it spent in writing gain solutions to disk\n";
  }
}

size_t BdaDdeCal::GetChanBlockIndex(const size_t channel,
                                    const size_t n_channels,
                                    size_t n_channel_blocks) const {
  // Assert that the channel averaging schema is compatible with the channel
  // block division. One channel should be mapped to one single channel block.
  if (n_channels < n_channel_blocks) {
    throw std::runtime_error(
        "Number of BDA channels smaller than numbner of channel blocks");
  }

  return std::floor(static_cast<double>(channel) / n_channels *
                    n_channel_blocks);
}

}  // namespace steps
}  // namespace dp3
