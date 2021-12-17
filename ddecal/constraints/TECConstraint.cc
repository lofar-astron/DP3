// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "TECConstraint.h"

#include <aocommon/parallelfor.h>

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

  phase_fitters_.resize(NThreads());
  for (PhaseFitter& fitter : phase_fitters_) fitter.Initialize(frequencies);

  weights_.assign(NChannelBlocks() * NAntennas(), 1.0);
  initializeChild();
}

void TECConstraintBase::SetWeights(const std::vector<double>& weights) {
  weights_ = weights;
}

void ApproximateTECConstraint::initializeChild() {
  pw_fitters_.resize(NThreads());
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

void TECConstraintBase::applyReferenceAntenna(
    std::vector<std::vector<dcomplex>>& solutions) const {
  // Choose reference antenna that has at least 20% channels unflagged
  size_t ref_antenna = 0;
  for (; ref_antenna != NAntennas(); ++ref_antenna) {
    size_t n_unflagged_channels = 0;
    // Only check flagged state for first direction
    for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
      if (isfinite(solutions[ch][ref_antenna * NDirections()]))
        n_unflagged_channels++;
    }
    if (n_unflagged_channels * 1.0 / NChannelBlocks() > 0.2)
      // Choose this refAntenna;
      break;
  }
  // All antennas are flagged, use first one (will lead to NaNs for this solint)
  if (ref_antenna == NAntennas()) ref_antenna = 0;

  for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
    for (size_t antenna_index = 0; antenna_index != NAntennas();
         ++antenna_index) {
      for (size_t d = 0; d != NDirections(); ++d) {
        size_t solution_index = antenna_index * NDirections() + d;
        size_t ref_antenna_index = d + ref_antenna * NDirections();
        if (antenna_index != ref_antenna) {
          solutions[ch][solution_index] =
              solutions[ch][solution_index] / solutions[ch][ref_antenna_index];
        }
      }
    }
    for (size_t d = 0; d != NDirections(); ++d)
      solutions[ch][ref_antenna * NDirections() + d] = 1.0;
  }
}

std::vector<Constraint::Result> TECConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  size_t nRes = 3;
  if (mode_ == TECOnlyMode) {
    nRes = 2;  // TEC and error
  } else {
    nRes = 3;  // TEC, phase and error
  }

  std::vector<Constraint::Result> res(nRes);
  res[0].vals.resize(NAntennas() * NDirections());
  res[0].weights.resize(NAntennas() * NDirections());
  res[0].axes = "ant,dir,freq";
  res[0].name = "tec";
  res[0].dims.resize(3);
  res[0].dims[0] = NAntennas();
  res[0].dims[1] = NDirections();
  res[0].dims[2] = 1;
  if (mode_ == TECAndCommonScalarMode) {
    res[1] = res[0];
    res[1].name = "phase";
  }
  res.back() = res[0];
  res.back().name = "error";

  // Divide out the reference antenna
  if (do_phase_reference_) applyReferenceAntenna(solutions);

  aocommon::ParallelFor<size_t> loop(NThreads());
  loop.Run(0, NAntennas() * NDirections(),
           [&](size_t solution_index, size_t thread) {
             size_t antenna_index = solution_index / NDirections();

             // Flag channels where calibration yielded inf or nan
             double weight_sum = 0.0;
             for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
               if (isfinite(solutions[ch][solution_index])) {
                 phase_fitters_[thread].PhaseData()[ch] =
                     std::arg(solutions[ch][solution_index]);
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
             if (mode_ == TECOnlyMode) {
               res.back().vals[solution_index] =
                   phase_fitters_[thread].FitDataToTEC1Model(alpha);
             } else {
               res.back().vals[solution_index] =
                   phase_fitters_[thread].FitDataToTEC2Model(alpha, beta);
             }
             res.back().weights[solution_index] = weight_sum;

             res[0].vals[solution_index] = alpha / -8.44797245e9;
             res[0].weights[solution_index] = weight_sum;
             if (mode_ == TECAndCommonScalarMode) {
               res[1].vals[solution_index] = beta;
               res[1].weights[solution_index] = weight_sum;
             }

             for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
               solutions[ch][solution_index] = std::polar<double>(
                   1.0, phase_fitters_[thread].PhaseData()[ch]);
             }
           });

  return res;
}

std::vector<Constraint::Result> ApproximateTECConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double time,
    std::ostream* stat_stream) {
  if (finished_approximate_stage_)
    return TECConstraint::Apply(solutions, time, stat_stream);
  else {
    if (do_phase_reference_) applyReferenceAntenna(solutions);

    aocommon::ParallelFor<size_t> loop(NThreads());
    loop.Run(0, NAntennas() * NDirections(),
             [&](size_t solutionIndex, size_t thread) {
               size_t antennaIndex = solutionIndex / NDirections();
               std::vector<double>& data = thread_data_[thread];
               std::vector<double>& fittedData = thread_fitted_data_[thread];
               std::vector<double>& weights = thread_weights_[thread];

               // Flag channels where calibration yielded inf or nan
               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 if (isfinite(solutions[ch][solutionIndex])) {
                   data[ch] = std::arg(solutions[ch][solutionIndex]);
                   weights[ch] = weights_[antennaIndex * NChannelBlocks() + ch];
                 } else {
                   data[ch] = 0.0;
                   weights[ch] = 0.0;
                 }
               }

               // TODO might be nice to make it a user option whether to break
               // or not
               pw_fitters_[thread].SlidingFitWithBreak(
                   phase_fitters_[thread].GetFrequencies().data(), data.data(),
                   weights.data(), fittedData.data(), data.size());

               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 solutions[ch][solutionIndex] =
                     std::polar<double>(1.0, fittedData[ch]);
               }
             });

    return std::vector<Constraint::Result>();
  }
}

}  // namespace ddecal
}  // namespace dp3
