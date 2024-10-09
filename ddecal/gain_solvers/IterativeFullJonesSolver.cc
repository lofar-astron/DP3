// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeFullJonesSolver.h"

#include <algorithm>

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/staticfor.h>

using aocommon::MC2x2F;

namespace dp3 {
namespace ddecal {

IterativeFullJonesSolver::SolveResult IterativeFullJonesSolver::Solve(
    const FullSolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  PrepareConstraints();

  SolutionTensor next_solutions(
      {NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()});

  SolveResult result;

  // Visibility vector v_residual[cb][vis] of size NChannelBlocks() x
  // n_visibilities
  std::vector<std::vector<MC2x2F>> v_residual(NChannelBlocks());
  // The following loop allocates all structures
  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
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
  step_magnitudes.reserve(GetMaxIterations());

  aocommon::RecursiveFor recursive_for;
  do {
    MakeSolutionsFinite4Pol(solutions);

    aocommon::StaticFor<size_t> loop;
    loop.Run(0, NChannelBlocks(), [&](size_t ch_block, size_t end_index) {
      for (; ch_block < end_index; ++ch_block) {
        PerformIteration(ch_block, data.ChannelBlock(ch_block),
                         v_residual[ch_block], solutions[ch_block],
                         next_solutions);
      }
    });

    Step(solutions, next_solutions);

    constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

    double avg_squared_diff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes);
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

void IterativeFullJonesSolver::PerformIteration(
    size_t ch_block, const FullSolveData::ChannelBlockData& cb_data,
    std::vector<MC2x2F>& v_residual, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction)
    AddOrSubtractDirection<false>(cb_data, v_residual, direction, solutions);

  const std::vector<MC2x2F> v_copy = v_residual;

  for (size_t direction = 0; direction != NDirections(); ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    AddOrSubtractDirection<true>(cb_data, v_residual, direction, solutions);

    SolveDirection(ch_block, cb_data, v_residual, direction, solutions,
                   next_solutions);
  }
}

void IterativeFullJonesSolver::SolveDirection(
    size_t ch_block, const FullSolveData::ChannelBlockData& cb_data,
    const std::vector<MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions, SolutionTensor& next_solutions) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  constexpr size_t n_solution_pols = 4;
  const uint32_t n_dir_solutions = cb_data.NSolutionsForDirection(direction);
  std::vector<MC2x2F> numerator(NAntennas() * n_dir_solutions, MC2x2F::Zero());
  std::vector<MC2x2F> denominator(NAntennas() * n_dir_solutions,
                                  MC2x2F::Zero());

  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  const uint32_t solution_index0 = cb_data.SolutionIndex(direction, 0);

  std::mutex mutex;
  aocommon::StaticFor<size_t> loop;
  loop.Run(
      0, n_visibilities, [&](size_t start_vis_index, size_t end_vis_index) {
        std::vector<MC2x2F> local_numerator(NAntennas() * n_dir_solutions,
                                            MC2x2F::Zero());
        std::vector<MC2x2F> local_denominator(NAntennas() * n_dir_solutions,
                                              MC2x2F::Zero());
        for (size_t vis_index = start_vis_index; vis_index != end_vis_index;
             ++vis_index) {
          const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
          const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
          const uint32_t solution_index =
              cb_data.SolutionIndex(direction, vis_index);

          const MC2x2F solution_ant_1(
              &solutions[(antenna_1 * NSolutions() + solution_index) *
                         n_solution_pols]);

          const MC2x2F solution_ant_2(
              &solutions[(antenna_2 * NSolutions() + solution_index) *
                         n_solution_pols]);

          const MC2x2F data(v_residual[vis_index]);
          const MC2x2F model(cb_data.ModelVisibility(direction, vis_index));

          const uint32_t rel_solution_index = solution_index - solution_index0;
          // Calculate the contribution of this baseline for antenna_1
          const MC2x2F cor_model_herm_1(solution_ant_2.MultiplyHerm(model));
          const uint32_t full_solution_1_index =
              antenna_1 * n_dir_solutions + rel_solution_index;

          // sum(D^H J M) [ sum(M^H J^H J M) ]^-1
          local_numerator[full_solution_1_index] +=
              static_cast<MC2x2F>(data * cor_model_herm_1);
          local_denominator[full_solution_1_index] += static_cast<MC2x2F>(
              HermTranspose(cor_model_herm_1) * cor_model_herm_1);

          // Calculate the contribution of this baseline for antenna_2
          const MC2x2F cor_model_2(solution_ant_1 * model);
          const uint32_t full_solution_2_index =
              antenna_2 * n_dir_solutions + rel_solution_index;
          // sum(D^H J M) [ sum(M^H J^H J M) ]^-1
          local_numerator[full_solution_2_index] +=
              static_cast<MC2x2F>(HermTranspose(data) * cor_model_2);
          local_denominator[full_solution_2_index] +=
              static_cast<MC2x2F>(HermTranspose(cor_model_2) * cor_model_2);
        }
        std::scoped_lock lock(mutex);
        for (size_t i = 0; i != numerator.size(); ++i)
          numerator[i] += local_numerator[i];
        for (size_t i = 0; i != denominator.size(); ++i)
          denominator[i] += local_denominator[i];
      });

  for (size_t ant = 0; ant != NAntennas(); ++ant) {
    for (uint32_t rel_sol = 0; rel_sol != n_dir_solutions; ++rel_sol) {
      const uint32_t solution_index = rel_sol + solution_index0;
      const uint32_t index = ant * n_dir_solutions + rel_sol;
      MC2x2F result;
      if (denominator[index].Invert())
        result = numerator[index] * denominator[index];
      else
        result = MC2x2F::NaN();
      result.AssignTo(&next_solutions(ch_block, ant, solution_index, 0));
    }
  }
}

template <bool Add>
void IterativeFullJonesSolver::AddOrSubtractDirection(
    const FullSolveData::ChannelBlockData& cb_data,
    std::vector<MC2x2F>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions) {
  constexpr size_t n_solution_polarizations = 4;
  const size_t n_visibilities = cb_data.NVisibilities();
  aocommon::StaticFor<size_t> loop;
  loop.Run(
      0, n_visibilities, [&](size_t start_vis_index, size_t end_vis_index) {
        for (size_t vis_index = start_vis_index; vis_index != end_vis_index;
             ++vis_index) {
          const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
          const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
          const uint32_t solution_index =
              cb_data.SolutionIndex(direction, vis_index);
          const MC2x2F solution_1(
              &solutions[(antenna_1 * NSolutions() + solution_index) *
                         n_solution_polarizations]);
          const MC2x2F solution_2(
              &solutions[(antenna_2 * NSolutions() + solution_index) *
                         n_solution_polarizations]);
          const MC2x2F model(cb_data.ModelVisibility(direction, vis_index));
          const MC2x2F term =
              static_cast<MC2x2F>(solution_1 * model.MultiplyHerm(solution_2));
          if (Add) {
            v_residual[vis_index] += term;
          } else {
            v_residual[vis_index] -= term;
          }
        }
      });
}

}  // namespace ddecal
}  // namespace dp3
