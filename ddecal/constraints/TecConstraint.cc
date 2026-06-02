// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "TecConstraint.h"

#include <schaapcommon/threading/dynamicfor.h>
#include <schaapcommon/threading/staticfor.h>

#include <xtensor/core/xmath.hpp>

namespace dp3::ddecal {

TecConstraint::TecConstraint(Mode mode)
    : mode_(mode),
      do_phase_reference_(true),
      phase_fitters_(),
      results_(mode_ == Mode::kTecOnly
                   ? 2
                   : 3)  // TEC and error, or TEC, phase and error
{}

void TecConstraint::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);

  const size_t n_threads = schaapcommon::ThreadPool::GetInstance().NThreads();
  phase_fitters_.resize(n_threads);
  for (PhaseFitter& fitter : phase_fitters_) fitter.Initialize(frequencies);

  weights_.assign(NChannelBlocks() * NAntennas(), 1.0);
  initializeChild();

  results_[0].vals.resize(NAntennas() * NSubSolutions());
  results_[0].weights.resize(NAntennas() * NSubSolutions());
  results_[0].axes = "ant,dir,freq";
  results_[0].name = "tec";
  results_[0].dims.resize(3);
  results_[0].dims[0] = NAntennas();
  results_[0].dims[1] = NSubSolutions();
  results_[0].dims[2] = 1;
  if (mode_ == Mode::kTecAndCommonScalar) {
    results_[1] = results_[0];
    results_[1].name = "phase";
  }
  results_.back() = results_[0];
  results_.back().name = "error";
}

void TecConstraint::SetWeights(const std::vector<double>& weights) {
  weights_ = weights;
}

void ApproximateTECConstraint::initializeChild() {
  const size_t n_threads = schaapcommon::ThreadPool::GetInstance().NThreads();
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

void TecConstraint::Apply(SolutionSpan& solutions, double time) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) ==
         1);  // single polarization (phase only) solutions

  // Divide out the reference antenna
  if (do_phase_reference_) ApplyReferenceAntenna(solutions);

  // Use a DynamicFor, since the ternarySearch functions in PhaseFitter.cc run
  // a variable number of iterations.
  schaapcommon::DynamicFor<size_t> loop;
  loop.Run(
      0, NAntennas() * NSubSolutions(),
      [&](size_t antenna_and_solution_index, size_t thread) {
        const size_t antenna_index =
            antenna_and_solution_index / NSubSolutions();
        const size_t sub_solution_index =
            antenna_and_solution_index % NSubSolutions();

        // Flag channels where calibration yielded inf or nan
        double weight_sum = 0.0;
        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          const std::complex<double>& solution =
              solutions(ch, antenna_index, sub_solution_index, 0);
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
          results_.back().vals[antenna_and_solution_index] =
              phase_fitters_[thread].FitDataToTEC1Model(alpha);
        } else {
          results_.back().vals[antenna_and_solution_index] =
              phase_fitters_[thread].FitDataToTEC2Model(alpha, beta);
        }
        results_.back().weights[antenna_and_solution_index] = weight_sum;

        results_[0].vals[antenna_and_solution_index] = alpha / -8.44797245e9;
        results_[0].weights[antenna_and_solution_index] = weight_sum;
        if (mode_ == Mode::kTecAndCommonScalar) {
          results_[1].vals[antenna_and_solution_index] = beta;
          results_[1].weights[antenna_and_solution_index] = weight_sum;
        }

        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          solutions(ch, antenna_index, sub_solution_index, 0) =
              std::polar<double>(1.0, phase_fitters_[thread].PhaseData()[ch]);
        }
      });
}

void ApproximateTECConstraint::Apply(SolutionSpan& solutions, double time) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) ==
         1);  // single polarization (phase only) solutions

  if (finished_approximate_stage_) {
    TecConstraint::Apply(solutions, time);
  } else {
    if (do_phase_reference_) ApplyReferenceAntenna(solutions);

    // Use a StaticFor, since PieceWisePhaseFitter.cc has no dynamic code.
    schaapcommon::StaticFor<size_t> loop;
    loop.Run(
        0, NAntennas() * NSubSolutions(),
        [&](size_t antenna_and_solution_index, size_t end_index,
            size_t thread) {
          for (; antenna_and_solution_index < end_index;
               ++antenna_and_solution_index) {
            const size_t antenna_index =
                antenna_and_solution_index / NSubSolutions();
            const size_t sub_solution_index =
                antenna_and_solution_index % NSubSolutions();
            std::vector<double>& data = thread_data_[thread];
            std::vector<double>& fitted_data = thread_fitted_data_[thread];
            std::vector<double>& weights = thread_weights_[thread];

            // Flag channels where calibration yielded inf or nan
            for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
              const std::complex<double>& solution =
                  solutions(ch, antenna_index, sub_solution_index, 0);
              if (isfinite(solution)) {
                data[ch] = std::arg(solution);
                weights[ch] = weights_[antenna_index * NChannelBlocks() + ch];
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
              solutions(ch, antenna_index, sub_solution_index, 0) =
                  std::polar<double>(1.0, fitted_data[ch]);
            }
          }
        });
  }
}

}  // namespace dp3::ddecal
