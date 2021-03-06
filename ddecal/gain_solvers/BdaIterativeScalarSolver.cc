// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaIterativeScalarSolver.h"

#include "SolverBuffer.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/parallelfor.h>

#include <algorithm>

namespace dp3 {
namespace ddecal {

BdaIterativeScalarSolver::SolveResult BdaIterativeScalarSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(n_channel_blocks_);

  SolveResult result;

  // Visibility vector ch_block x Nvis
  constexpr size_t n_solution_pols = 1;
  std::vector<std::vector<aocommon::MC2x2F>> v_residual(n_channel_blocks_);
  // The following loop allocates all structures
  for (size_t ch_block = 0; ch_block != n_channel_blocks_; ++ch_block) {
    next_solutions[ch_block].resize(n_directions_ * n_antennas_);
    v_residual[ch_block].resize(data.ChannelBlock(ch_block).NVisibilities());
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(max_iterations_);

  do {
    MakeSolutionsFinite1Pol(solutions);

    aocommon::ParallelFor<size_t> loop(n_threads_);
    loop.Run(0, n_channel_blocks_, [&](size_t ch_block, size_t /*thread*/) {
      PerformIteration(data.ChannelBlock(ch_block), v_residual[ch_block],
                       solutions[ch_block], next_solutions[ch_block]);
    });

    Step(solutions, next_solutions);

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes, n_solution_pols);
    iteration++;

    has_previously_converged = has_converged || has_previously_converged;

  } while (!ReachedStoppingCriterion(iteration, has_converged,
                                     constraints_satisfied, step_magnitudes));

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged iterations can be distinguished from converged ones.
  if (has_converged && constraints_satisfied)
    result.iterations = iteration;
  else
    result.iterations = iteration + 1;
  return result;
}

void BdaIterativeScalarSolver::PerformIteration(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != n_directions_; ++direction)
    AddOrSubtractDirection<false>(cb_data, v_residual, direction, solutions);

  const std::vector<aocommon::MC2x2F> v_copy = v_residual;

  for (size_t direction = 0; direction != n_directions_; ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    AddOrSubtractDirection<true>(cb_data, v_residual, direction, solutions);

    SolveDirection(cb_data, v_residual, direction, solutions, next_solutions);
  }
}

void BdaIterativeScalarSolver::SolveDirection(
    const SolveData::ChannelBlockData& cb_data,
    const std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  std::vector<std::complex<double>> numerator(n_antennas_, 0.0);
  std::vector<double> denominator(n_antennas_, 0.0);

  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  const std::vector<aocommon::MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    const Complex solution_ant1 =
        Complex(solutions[antenna1 * n_directions_ + direction]);
    const Complex solution_ant2 =
        Complex(solutions[antenna2 * n_directions_ + direction]);
    const aocommon::MC2x2F& data = v_residual[vis_index];
    const aocommon::MC2x2F& model = model_vector[vis_index];

    // Calculate the contribution of this baseline for antenna1
    const aocommon::MC2x2F cor_model_herm1(HermTranspose(model) *
                                           solution_ant2);
    numerator[antenna1] += Trace(data * cor_model_herm1);
    denominator[antenna1] += Norm(cor_model_herm1);

    // Calculate the contribution of this baseline for antenna2
    const aocommon::MC2x2F cor_model2(model * solution_ant1);
    numerator[antenna2] += Trace(HermTranspose(data) * cor_model2);
    denominator[antenna2] += Norm(cor_model2);
  }

  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    DComplex& destination = next_solutions[ant * n_directions_ + direction];

    if (denominator[ant] == 0.0)
      destination = std::numeric_limits<float>::quiet_NaN();
    else
      destination = numerator[ant] / denominator[ant];
  }
}

template <bool Add>
void BdaIterativeScalarSolver::AddOrSubtractDirection(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions) {
  const std::vector<aocommon::MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  const size_t n_visibilities = cb_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    const Complex solution1(solutions[antenna1 * n_directions_ + direction]);
    const Complex solution2Conj =
        std::conj(Complex(solutions[antenna2 * n_directions_ + direction]));
    aocommon::MC2x2F& data = v_residual[vis_index];
    const aocommon::MC2x2F& model = model_vector[vis_index];
    if (Add) {
      data[0] += solution1 * model[0] * solution2Conj;
      data[1] += solution1 * model[1] * solution2Conj;
      data[2] += solution1 * model[2] * solution2Conj;
      data[3] += solution1 * model[3] * solution2Conj;
    } else {
      data[0] -= solution1 * model[0] * solution2Conj;
      data[1] -= solution1 * model[1] * solution2Conj;
      data[2] -= solution1 * model[2] * solution2Conj;
      data[3] -= solution1 * model[3] * solution2Conj;
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
