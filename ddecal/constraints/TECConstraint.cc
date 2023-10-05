// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "TECConstraint.h"

#include <aocommon/dynamicfor.h>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

namespace dp3 {
namespace ddecal {

TECConstraintBase::TECConstraintBase(Mode mode)
    : mode_(mode), do_phase_reference_(true), phase_fitters_() {}

void TECConstraintBase::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);

  for (int32_t v : solutions_per_direction) {
    if (v != 1)
      throw std::runtime_error(
          "The TEC constraints do not yet support direction-dependent "
          "intervals");
  }

  const size_t n_threads = aocommon::ThreadPool::GetInstance().NThreads();
  phase_fitters_.resize(n_threads);
  for (PhaseFitter& fitter : phase_fitters_) fitter.Initialize(frequencies);

  weights_.assign(NChannelBlocks() * NAntennas(), 1.0);
  initializeChild();
}

void TECConstraintBase::SetWeights(const std::vector<double>& weights) {
  weights_ = weights;
}

void ApproximateTECConstraint::initializeChild() {
  const size_t n_threads = aocommon::ThreadPool::GetInstance().NThreads();
  pw_fitters_.resize(n_threads);
  thread_data_.resize(pw_fitters_.size());
  thread_fitted_data_.resize(pw_fitters_.size());
  thread_weights_.resize(pw_fitters_.size());
  for (size_t threadId = 0; threadId != pw_fitters_.size(); ++threadId) {
    thread_data_[threadId].resize(NChannelBlocks());
    thread_fitted_data_[threadId].resize(NChannelBlocks());
    thread_weights_[threadId].resize(NChannelBlocks());
  }

  if (fitting_chunk_size_ == 0) {
    fitting_chunk_size_ = PieceWisePhaseFitter::CalculateChunkSize(
        phase_fitters_.front().GetFrequencies());
  }
  for (size_t i = 0; i != pw_fitters_.size(); ++i)
    pw_fitters_[i].SetChunkSize(fitting_chunk_size_);
}

void TECConstraintBase::applyReferenceAntenna(SolutionSpan& solutions) const {
  // Choose reference antenna that has at least 20% channels unflagged
  size_t ref_antenna = 0;
  for (; ref_antenna != NAntennas(); ++ref_antenna) {
    // Only check flagged state for first direction
    const size_t n_unflagged_channels = xt::sum(
        xt::isfinite(xt::view(solutions, xt::all(), ref_antenna, 0, 0)))();
    if (n_unflagged_channels > 0.2 * NChannelBlocks())
      // Choose this refAntenna;
      break;
  }
  // All antennas are flagged, use first one (will lead to NaNs for this solint)
  if (ref_antenna == NAntennas()) ref_antenna = 0;

  for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
    auto reference_view = xt::view(solutions, ch, ref_antenna, xt::all(), 0);
    for (size_t antenna_index = 0; antenna_index != NAntennas();
         ++antenna_index) {
      if (antenna_index != ref_antenna) {
        xt::view(solutions, ch, antenna_index, xt::all(), 0) /= reference_view;
      }
    }
    reference_view.fill(1.0);
  }
}

std::vector<Constraint::Result> TECConstraint::Apply(
    SolutionSpan& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NSolutions());
  assert(solutions.shape(3) ==
         1);  // single polarization (phase only) solutions

  size_t nRes = 3;
  if (mode_ == Mode::kTecOnly) {
    nRes = 2;  // TEC and error
  } else {
    nRes = 3;  // TEC, phase and error
  }

  std::vector<Constraint::Result> res(nRes);
  res[0].vals.resize(NAntennas() * NSolutions());
  res[0].weights.resize(NAntennas() * NSolutions());
  res[0].axes = "ant,dir,freq";
  res[0].name = "tec";
  res[0].dims.resize(3);
  res[0].dims[0] = NAntennas();
  res[0].dims[1] = NSolutions();
  res[0].dims[2] = 1;
  if (mode_ == Mode::kTecAndCommonScalar) {
    res[1] = res[0];
    res[1].name = "phase";
  }
  res.back() = res[0];
  res.back().name = "error";

  // Divide out the reference antenna
  if (do_phase_reference_) applyReferenceAntenna(solutions);

  aocommon::DynamicFor<size_t> loop;
  loop.Run(
      0, NAntennas() * NSolutions(),
      [&](size_t antenna_and_solution_index, size_t thread) {
        const size_t antenna_index = antenna_and_solution_index / NSolutions();
        const size_t solution_index = antenna_and_solution_index % NSolutions();

        // Flag channels where calibration yielded inf or nan
        double weight_sum = 0.0;
        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          const std::complex<double>& solution =
              solutions(ch, antenna_index, solution_index, 0);
          if (isfinite(solution)) {
            phase_fitters_[thread].PhaseData()[ch] = std::arg(solution);
            phase_fitters_[thread].WeightData()[ch] =
                weights_[antenna_index * NChannelBlocks() + ch];
            weight_sum += weights_[antenna_index * NChannelBlocks() + ch];
          } else {
            phase_fitters_[thread].PhaseData()[ch] = 0.0;
            phase_fitters_[thread].WeightData()[ch] = 0.0;
          }
        }

        double alpha;
        double beta = 0.0;
        if (mode_ == Mode::kTecOnly) {
          res.back().vals[antenna_and_solution_index] =
              phase_fitters_[thread].FitDataToTEC1Model(alpha);
        } else {
          res.back().vals[antenna_and_solution_index] =
              phase_fitters_[thread].FitDataToTEC2Model(alpha, beta);
        }
        res.back().weights[antenna_and_solution_index] = weight_sum;

        res[0].vals[antenna_and_solution_index] = alpha / -8.44797245e9;
        res[0].weights[antenna_and_solution_index] = weight_sum;
        if (mode_ == Mode::kTecAndCommonScalar) {
          res[1].vals[antenna_and_solution_index] = beta;
          res[1].weights[antenna_and_solution_index] = weight_sum;
        }

        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          solutions(ch, antenna_index, solution_index, 0) =
              std::polar<double>(1.0, phase_fitters_[thread].PhaseData()[ch]);
        }
      });

  return res;
}

std::vector<Constraint::Result> ApproximateTECConstraint::Apply(
    SolutionSpan& solutions, double time, std::ostream* stat_stream) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NSolutions());
  assert(solutions.shape(3) ==
         1);  // single polarization (phase only) solutions

  if (finished_approximate_stage_) {
    return TECConstraint::Apply(solutions, time, stat_stream);
  } else {
    if (do_phase_reference_) applyReferenceAntenna(solutions);

    aocommon::DynamicFor<size_t> loop;
    loop.Run(0, NAntennas() * NSolutions(),
             [&](size_t antenna_and_solution_index, size_t thread) {
               const size_t antenna_index =
                   antenna_and_solution_index / NSolutions();
               const size_t solution_index =
                   antenna_and_solution_index % NSolutions();
               std::vector<double>& data = thread_data_[thread];
               std::vector<double>& fitted_data = thread_fitted_data_[thread];
               std::vector<double>& weights = thread_weights_[thread];

               // Flag channels where calibration yielded inf or nan
               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 const std::complex<double>& solution =
                     solutions(ch, antenna_index, solution_index, 0);
                 if (isfinite(solution)) {
                   data[ch] = std::arg(solution);
                   weights[ch] =
                       weights_[antenna_index * NChannelBlocks() + ch];
                 } else {
                   data[ch] = 0.0;
                   weights[ch] = 0.0;
                 }
               }

               // TODO might be nice to make it a user option whether to break
               // or not
               pw_fitters_[thread].SlidingFitWithBreak(
                   phase_fitters_[thread].GetFrequencies().data(), data.data(),
                   weights.data(), fitted_data.data(), data.size());

               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 solutions(ch, antenna_index, solution_index, 0) =
                     std::polar<double>(1.0, fitted_data[ch]);
               }
             });

    return std::vector<Constraint::Result>();
  }
}

}  // namespace ddecal
}  // namespace dp3
