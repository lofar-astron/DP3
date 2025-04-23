// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeDiagonalSolver.h"

#include <algorithm>

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>
#include <aocommon/recursivefor.h>
#include <aocommon/staticfor.h>

using aocommon::MC2x2F;
using aocommon::MC2x2FDiag;

namespace dp3::ddecal {

namespace {
void AddNormToDenominator(float* denominator, const MC2x2F& corrected_model) {
  // The indices (0, 2 ; 1, 3) are following from the fact that we want
  // the contribution of antenna2's "X" polarization, and the matrix is
  // ordered [ XX XY ; YX YY ].
  denominator[0] +=
      std::norm(corrected_model.Get(0)) + std::norm(corrected_model.Get(2));
  denominator[1] +=
      std::norm(corrected_model.Get(1)) + std::norm(corrected_model.Get(3));
}

void AddNormToDenominator(float* denominator,
                          const MC2x2FDiag& corrected_model) {
  denominator[0] += std::norm(corrected_model.Get(0));
  denominator[1] += std::norm(corrected_model.Get(1));
}
}  // namespace

template <typename VisMatrix>
SolverBase::SolveResult IterativeDiagonalSolver<VisMatrix>::Solve(
    const SolveData<VisMatrix>& data,
    std::vector<std::vector<DComplex>>& solutions, double time,
    std::ostream* stat_stream) {
  PrepareConstraints();

  SolutionTensor next_solutions({NChannelBlocks(), NAntennas(), NSubSolutions(),
                                 NSolutionPolarizations()});

  SolveResult result;

  // Visibility vector v_residual[cb][vis] of size NChannelBlocks() x
  // n_visibilities
  std::vector<std::vector<VisMatrix>> v_residual(NChannelBlocks());
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

  std::unique_ptr<aocommon::RecursiveFor> recursive_for =
      MakeOptionalRecursiveFor();
  do {
    MakeSolutionsFinite2Pol(solutions);

    aocommon::RunStaticFor<size_t>(
        0, NChannelBlocks(), [&](size_t start_block, size_t end_block) {
          for (size_t ch_block = start_block; ch_block < end_block;
               ++ch_block) {
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

template <typename VisMatrix>
void IterativeDiagonalSolver<VisMatrix>::PerformIteration(
    size_t ch_block, const ChannelBlockData& cb_data,
    std::vector<VisMatrix>& v_residual, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction)
    DiagonalAddOrSubtractDirection<false, VisMatrix>(cb_data, v_residual,
                                                     direction, NSubSolutions(),
                                                     solutions, NSubThreads());

  const std::vector<VisMatrix> v_copy = v_residual;

  for (size_t direction = 0; direction != NDirections(); ++direction) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving
    if (direction != 0) v_residual = v_copy;
    DiagonalAddOrSubtractDirection<true>(cb_data, v_residual, direction,
                                         NSubSolutions(), solutions,
                                         NSubThreads());

    SolveDirection(ch_block, cb_data, v_residual, direction, solutions,
                   next_solutions);
  }
}

template <typename VisMatrix>
void IterativeDiagonalSolver<VisMatrix>::SolveDirection(
    size_t ch_block, const ChannelBlockData& cb_data,
    const std::vector<VisMatrix>& v_residual, size_t direction,
    const std::vector<DComplex>& solutions, SolutionTensor& next_solutions) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)

  const uint32_t n_dir_solutions = cb_data.NSolutionsForDirection(direction);
  std::vector<MC2x2FDiag> numerator(NAntennas() * n_dir_solutions,
                                    MC2x2FDiag::Zero());
  std::vector<float> denominator(NAntennas() * n_dir_solutions * 2, 0.0);

  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  const uint32_t solution_index0 = cb_data.SolutionIndex(direction, 0);

  std::mutex mutex;
  aocommon::RunConstrainedStaticFor<size_t>(
      0u, n_visibilities, NSubThreads(),
      [&](size_t start_vis_index, size_t end_vis_index) {
        std::vector<MC2x2FDiag> local_numerator(NAntennas() * n_dir_solutions,
                                                MC2x2FDiag::Zero());
        std::vector<float> local_denominator(NAntennas() * n_dir_solutions * 2,
                                             0.0);
        for (size_t vis_index = start_vis_index; vis_index != end_vis_index;
             ++vis_index) {
          const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
          const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
          const uint32_t solution_index =
              cb_data.SolutionIndex(direction, vis_index);
          const DComplex* solution_ant_1 =
              &solutions[(antenna_1 * NSubSolutions() + solution_index) * 2];
          const DComplex* solution_ant_2 =
              &solutions[(antenna_2 * NSubSolutions() + solution_index) * 2];
          const VisMatrix& data = v_residual[vis_index];
          const VisMatrix& model =
              cb_data.ModelVisibility(direction, vis_index);

          const uint32_t rel_solution_index = solution_index - solution_index0;
          // Calculate the contribution of this baseline for antenna_1
          const MC2x2FDiag solution_1{Complex(solution_ant_2[0]),
                                      Complex(solution_ant_2[1])};
          const VisMatrix cor_model_transp_1(solution_1 * HermTranspose(model));
          const uint32_t full_solution_1_index =
              antenna_1 * n_dir_solutions + rel_solution_index;
          local_numerator[full_solution_1_index] +=
              Diagonal(data * cor_model_transp_1);
          AddNormToDenominator(&local_denominator[full_solution_1_index * 2],
                               cor_model_transp_1);

          // Calculate the contribution of this baseline for antenna_2
          // data_ba = data_ab^H, etc., therefore, numerator and denominator
          // become:
          // - num = data_ab^H * solutions_a * model_ab
          // - den = norm(model_ab^H * solutions_a)
          const MC2x2FDiag solution_2{Complex(solution_ant_1[0]),
                                      Complex(solution_ant_1[1])};
          const VisMatrix cor_model_2(solution_2 * model);

          const uint32_t full_solution_2_index =
              antenna_2 * n_dir_solutions + rel_solution_index;
          local_numerator[full_solution_2_index] +=
              Diagonal(HermTranspose(data) * cor_model_2);
          AddNormToDenominator(&local_denominator[full_solution_2_index * 2],
                               cor_model_2);
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

      for (size_t pol = 0; pol != 2; ++pol) {
        if (denominator[index * 2 + pol] == 0.0)
          next_solutions(ch_block, ant, solution_index, pol) =
              std::numeric_limits<double>::quiet_NaN();
        else
          next_solutions(ch_block, ant, solution_index, pol) =
              DComplex(numerator[index].Get(pol)) /
              double(denominator[index * 2 + pol]);
      }
    }
  }
}

template class IterativeDiagonalSolver<aocommon::MC2x2F>;
template class IterativeDiagonalSolver<aocommon::MC2x2FDiag>;

}  // namespace dp3::ddecal
