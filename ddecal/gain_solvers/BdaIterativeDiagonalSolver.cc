// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BdaIterativeDiagonalSolver.h"

#include "SolverBuffer.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/parallelfor.h>

#include <algorithm>

namespace dp3 {
namespace ddecal {

BdaIterativeDiagonalSolver::SolveResult BdaIterativeDiagonalSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  PrepareConstraints();

  std::vector<std::vector<DComplex>> next_solutions(n_channel_blocks_);

  SolveResult result;

  // Visibility vector ch_block x Nvis
  std::vector<std::vector<aocommon::MC2x2F>> v_residual(n_channel_blocks_);
  // The following loop allocates all structures
  for (size_t ch_block = 0; ch_block != n_channel_blocks_; ++ch_block) {
    next_solutions[ch_block].resize(NDirections() * NAntennas() * 2);
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
    MakeSolutionsFinite2Pol(solutions);

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
    constexpr size_t n_solution_pols = 2;
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

void BdaIterativeDiagonalSolver::PerformIteration(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction)
    AddOrSubtractDirection<false>(cb_data, v_residual, direction, solutions);

  const std::vector<aocommon::MC2x2F> v_copy = v_residual;

  for (size_t direction = 0; direction != NDirections(); ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    AddOrSubtractDirection<true>(cb_data, v_residual, direction, solutions);

    SolveDirection(cb_data, v_residual, direction, solutions, next_solutions);
  }
}

void BdaIterativeDiagonalSolver::SolveDirection(
    const SolveData::ChannelBlockData& cb_data,
    const std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions,
    std::vector<DComplex>& next_solutions) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  std::vector<aocommon::MC2x2FDiag> numerator(NAntennas(),
                                              aocommon::MC2x2FDiag::Zero());
  std::vector<float> denominator(NAntennas() * 2, 0.0);

  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  const std::vector<aocommon::MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    const DComplex* solution_ant1 =
        &solutions[(antenna1 * NDirections() + direction) * 2];
    const DComplex* solution_ant2 =
        &solutions[(antenna2 * NDirections() + direction) * 2];
    const aocommon::MC2x2F& data = v_residual[vis_index];
    const aocommon::MC2x2F& model = model_vector[vis_index];

    // Calculate the contribution of this baseline for antenna1
    const aocommon::MC2x2FDiag solution1{Complex(solution_ant2[0]),
                                         Complex(solution_ant2[1])};
    const aocommon::MC2x2F cor_model_transp1(solution1 * HermTranspose(model));
    numerator[antenna1] += Diagonal(data * cor_model_transp1);
    // The indices (0, 2 / 1, 3) are following from the fact that we want
    // the contribution of antenna2's "X" polarization, and the matrix is
    // ordered [ XX XY / YX YY ].
    denominator[antenna1 * 2] +=
        std::norm(cor_model_transp1[0]) + std::norm(cor_model_transp1[2]);
    denominator[antenna1 * 2 + 1] +=
        std::norm(cor_model_transp1[1]) + std::norm(cor_model_transp1[3]);

    // Calculate the contribution of this baseline for antenna2
    // data_ba = data_ab^H, etc., therefore, numerator and denominator
    // become:
    // - num = data_ab^H * solutions_a * model_ab
    // - den = norm(model_ab^H * solutions_a)
    const aocommon::MC2x2FDiag solution2{Complex(solution_ant1[0]),
                                         Complex(solution_ant1[1])};
    const aocommon::MC2x2F cor_model2(solution2 * model);

    numerator[antenna2] += Diagonal(HermTranspose(data) * cor_model2);
    denominator[antenna2 * 2] +=
        std::norm(cor_model2[0]) + std::norm(cor_model2[2]);
    denominator[antenna2 * 2 + 1] +=
        std::norm(cor_model2[1]) + std::norm(cor_model2[3]);
  }

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    DComplex* destination =
        &next_solutions[(ant * NDirections() + direction) * 2];

    for (size_t pol = 0; pol != 2; ++pol) {
      if (denominator[ant * 2 + pol] == 0.0)
        destination[pol] = std::numeric_limits<double>::quiet_NaN();
      else
        destination[pol] =
            DComplex(numerator[ant][pol]) / double(denominator[ant * 2 + pol]);
    }
  }
}

template <bool Add>
void BdaIterativeDiagonalSolver::AddOrSubtractDirection(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions) {
  const std::vector<aocommon::MC2x2F>& model_vector =
      cb_data.ModelVisibilityVector(direction);
  const size_t n_visibilities = cb_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const size_t antenna1 = cb_data.Antenna1Index(vis_index);
    const size_t antenna2 = cb_data.Antenna2Index(vis_index);
    const DComplex* solution1 =
        &solutions[(antenna1 * NDirections() + direction) * 2];
    const DComplex* solution2 =
        &solutions[(antenna2 * NDirections() + direction) * 2];
    const Complex solution1_0(solution1[0]);
    const Complex solution1_1(solution1[1]);
    const Complex solution2_0_conj(std::conj(solution2[0]));
    const Complex solution2_1_conj(std::conj(solution2[1]));

    aocommon::MC2x2F& data = v_residual[vis_index];
    const aocommon::MC2x2F& model = model_vector[vis_index];
    const aocommon::MC2x2F contribution(
        solution1_0 * model[0] * solution2_0_conj,
        solution1_0 * model[1] * solution2_1_conj,
        solution1_1 * model[2] * solution2_0_conj,
        solution1_1 * model[3] * solution2_1_conj);
    if (Add)
      data += contribution;
    else
      data -= contribution;
  }
}

}  // namespace ddecal
}  // namespace dp3
