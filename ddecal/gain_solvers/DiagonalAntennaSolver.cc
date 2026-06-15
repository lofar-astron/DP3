#include "DiagonalAntennaSolver.h"

#include <cassert>

#include <schaapcommon/threading/dynamicfor.h>

#include <gsl/gsl_multilarge_nlinear.h>

using aocommon::MC2x2;
using aocommon::MC2x2Diag;
using aocommon::MC2x2F;

namespace dp3::ddecal {

namespace {
const MC2x2 GetVisibility(const gsl_vector* x, size_t visibility_index) {
  return MC2x2{{gsl_vector_get(x, visibility_index * 8 + 0),
                gsl_vector_get(x, visibility_index * 8 + 1)},
               {gsl_vector_get(x, visibility_index * 8 + 2),
                gsl_vector_get(x, visibility_index * 8 + 3)},
               {gsl_vector_get(x, visibility_index * 8 + 4),
                gsl_vector_get(x, visibility_index * 8 + 5)},
               {gsl_vector_get(x, visibility_index * 8 + 6),
                gsl_vector_get(x, visibility_index * 8 + 7)}};
}

void SetVisibility(gsl_vector* x, size_t visibility_index, const MC2x2& m) {
  for (size_t i = 0; i != 4; ++i) {
    gsl_vector_set(x, visibility_index * 8 + i * 2, m.Get(i).real());
    gsl_vector_set(x, visibility_index * 8 + i * 2 + 1, m.Get(i).imag());
  }
}

MC2x2Diag GetSolution(const gsl_vector* solutions, size_t solution_index) {
  return MC2x2Diag{{gsl_vector_get(solutions, solution_index * 4),
                    gsl_vector_get(solutions, solution_index * 4 + 1)},
                   {gsl_vector_get(solutions, solution_index * 4 + 2),
                    gsl_vector_get(solutions, solution_index * 4 + 3)}};
}

void AddSolution(gsl_vector* solutions, size_t solution_index,
                 const MC2x2Diag& added_term) {
  const size_t index = solution_index * 4;
  gsl_vector_set(solutions, index,
                 added_term.Get(0).real() + gsl_vector_get(solutions, index));
  gsl_vector_set(
      solutions, index + 1,
      added_term.Get(0).imag() + gsl_vector_get(solutions, index + 1));
  gsl_vector_set(
      solutions, index + 2,
      added_term.Get(1).real() + gsl_vector_get(solutions, index + 2));
  gsl_vector_set(
      solutions, index + 3,
      added_term.Get(1).imag() + gsl_vector_get(solutions, index + 3));
}

void AddSolutionToMatrix(gsl_matrix* matrix, size_t row, size_t col,
                         const MC2x2& added_term) {
  const size_t i = row * 4;
  const size_t j = col * 4;
  const auto real_part = [&](size_t index) {
    return added_term.Get(index).real();
  };
  const auto imag_part = [&](size_t index) {
    return added_term.Get(index).imag();
  };
  // XX block (0-1 <-> 0-1)
  gsl_matrix_set(matrix, i + 0, j + 0,
                 gsl_matrix_get(matrix, i + 0, j + 0) + real_part(0));
  gsl_matrix_set(matrix, i + 0, j + 1,
                 gsl_matrix_get(matrix, i + 0, j + 1) + imag_part(0));
  gsl_matrix_set(
      matrix, i + 1, j + 0,
      gsl_matrix_get(matrix, i + 1, j + 0) - imag_part(0));  // Hermitian
  gsl_matrix_set(matrix, i + 1, j + 1,
                 gsl_matrix_get(matrix, i + 1, j + 1) + real_part(0));

  // YY block (2-3 <-> 2-3)
  gsl_matrix_set(matrix, i + 2, j + 2,
                 gsl_matrix_get(matrix, i + 2, j + 2) + real_part(3));
  gsl_matrix_set(matrix, i + 2, j + 3,
                 gsl_matrix_get(matrix, i + 2, j + 3) + imag_part(3));
  gsl_matrix_set(matrix, i + 3, j + 2,
                 gsl_matrix_get(matrix, i + 3, j + 2) - imag_part(3));
  gsl_matrix_set(matrix, i + 3, j + 3,
                 gsl_matrix_get(matrix, i + 3, j + 3) + real_part(3));

  // XY block (rows 0-1, cols 2-3)
  gsl_matrix_set(matrix, i + 0, j + 2,
                 gsl_matrix_get(matrix, i + 0, j + 2) + real_part(1));
  gsl_matrix_set(matrix, i + 0, j + 3,
                 gsl_matrix_get(matrix, i + 0, j + 3) + imag_part(1));
  gsl_matrix_set(matrix, i + 1, j + 2,
                 gsl_matrix_get(matrix, i + 1, j + 2) - imag_part(1));
  gsl_matrix_set(matrix, i + 1, j + 3,
                 gsl_matrix_get(matrix, i + 1, j + 3) + real_part(1));

  // YX block (rows 2-3, cols 0-1)
  gsl_matrix_set(matrix, i + 2, j + 0,
                 gsl_matrix_get(matrix, i + 2, j + 0) + real_part(2));
  gsl_matrix_set(matrix, i + 2, j + 1,
                 gsl_matrix_get(matrix, i + 2, j + 1) + imag_part(2));
  gsl_matrix_set(matrix, i + 3, j + 0,
                 gsl_matrix_get(matrix, i + 3, j + 0) - imag_part(2));
  gsl_matrix_set(matrix, i + 3, j + 1,
                 gsl_matrix_get(matrix, i + 3, j + 1) + real_part(2));
}

}  // namespace

int DiagonalAntennaSolver::PenaltyDerivativeForward(
    CBLAS_TRANSPOSE_t transpose_j, const gsl_vector* x, const gsl_vector* u,
    void* solver_info, gsl_vector* v, gsl_matrix* jtj) {
  SolverInfo& si = *static_cast<SolverInfo*>(solver_info);
  if (jtj != nullptr) {
    si.solver.PenaltyJTJ(x, si, jtj);
  }
  if (transpose_j == CblasTrans)
    return si.solver.PenaltyTransposedDerivative(x, u, si, v);
  else
    return si.solver.PenaltyDerivative(x, u, si, v);
}

int DiagonalAntennaSolver::PenaltyVVForward(const gsl_vector* x,
                                            const gsl_vector* v,
                                            void* solver_info,
                                            gsl_vector* fvv) {
  SolverInfo& si = *static_cast<SolverInfo*>(solver_info);
  return si.solver.PenaltyVV(x, v, si, fvv);
}

int DiagonalAntennaSolver::PenaltyFunctionForward(const gsl_vector* x,
                                                  void* solver_info,
                                                  gsl_vector* f) {
  SolverInfo& si = *static_cast<SolverInfo*>(solver_info);
  return si.solver.PenaltyFunction(x, si, f);
}

DiagonalAntennaSolver::SolveResult DiagonalAntennaSolver::Solve(
    const FullSolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time) {
  PrepareConstraints();

  const bool subtract_immediately = GetStepSize() > 0.99;
  if (subtract_immediately) {
    CalculateNormPerDirection(data);
  } else {
    direction_ordering_.reserve(NDirections());
    for (size_t i = 0; i != NDirections(); ++i)
      direction_ordering_.emplace_back(i);
  }

  SolutionTensor next_solutions({NChannelBlocks(), NAntennas(), NSubSolutions(),
                                 NSolutionPolarizations()});

  SolveResult result;

  // Visibility vector residual_data[cb][vis] is of size NChannelBlocks() x
  // n_visibilities
  std::vector<std::vector<MC2x2F>> residual_data(NChannelBlocks());
  // Allocates all structures
  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    residual_data[ch_block].resize(data.ChannelBlock(ch_block).NVisibilities());
  }

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(GetMaxIterations());
  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;
  do {
    MakeSolutionsFinite2Pol(solutions);

    schaapcommon::DynamicFor<size_t> loop;
    loop.Run(0, NChannelBlocks(), [&](size_t ch_block) {
      PerformIteration(ch_block, data.ChannelBlock(ch_block),
                       residual_data[ch_block], solutions[ch_block],
                       next_solutions, iteration);
    });

    Step(solutions, next_solutions);

    constraints_satisfied = ApplyConstraints(
        iteration, time, has_previously_converged, next_solutions);

    double avg_squared_diff;
    bool has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avg_squared_diff, step_magnitudes);
    iteration++;

    has_previously_converged = has_converged || has_previously_converged;

  } while (!ReachedStoppingCriterion(iteration, has_converged,
                                     constraints_satisfied, step_magnitudes));

  if (has_converged && constraints_satisfied)
    result.iterations = iteration;
  else
    result.iterations = iteration + 1;
  return result;
}

void DiagonalAntennaSolver::CalculateNormPerDirection(
    const FullSolveData& data) {
  std::vector<std::pair<float, size_t>> norm_sum_per_direction(NDirections());
  for (size_t channel_block_index = 0;
       channel_block_index != data.NChannelBlocks(); ++channel_block_index) {
    const FullSolveData::ChannelBlockData& ch_block_data =
        data.ChannelBlock(channel_block_index);
    const size_t n_visibilities = ch_block_data.NVisibilities();
    for (size_t direction = 0; direction != ch_block_data.NDirections();
         ++direction) {
      norm_sum_per_direction[direction].second = direction;
      for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
        const aocommon::MC2x2F& model =
            ch_block_data.ModelVisibility(direction, vis_index);
        if (ch_block_data.Weight(vis_index) != 0.0)
          norm_sum_per_direction[direction].first += Norm(model);
      }
    }
  }
  std::sort(norm_sum_per_direction.begin(), norm_sum_per_direction.end(),
            std::greater<std::pair<float, size_t>>());
  direction_ordering_.reserve(NDirections());
  for (const std::pair<float, size_t>& direction : norm_sum_per_direction) {
    direction_ordering_.emplace_back(direction.second);
  }
}

double DiagonalAntennaSolver::ChiSquared(
    const FullSolveData::ChannelBlockData& ch_block_data,
    std::vector<aocommon::MC2x2F>& residual_data, size_t direction,
    const SolutionSpan& solutions) const {
  double chi_squared = 0.0;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t antenna_1 = ch_block_data.Antenna1Index(vis_index);
    const uint32_t antenna_2 = ch_block_data.Antenna2Index(vis_index);
    const uint32_t solution_index =
        ch_block_data.SolutionIndex(direction, vis_index);
    const DComplex* solution_1 =
        &solutions[(antenna_1 * NSubSolutions() + solution_index) * 2];
    const DComplex* solution_2 =
        &solutions[(antenna_2 * NSubSolutions() + solution_index) * 2];
    const Complex solution_1_0(solution_1[0]);
    const Complex solution_1_1(solution_1[1]);
    const Complex solution_2_0_conj(std::conj(solution_2[0]));
    const Complex solution_2_1_conj(std::conj(solution_2[1]));

    const MC2x2F& model = ch_block_data.ModelVisibility(direction, vis_index);
    const MC2x2F contribution(solution_1_0 * model.Get(0) * solution_2_0_conj,
                              solution_1_0 * model.Get(1) * solution_2_1_conj,
                              solution_1_1 * model.Get(2) * solution_2_0_conj,
                              solution_1_1 * model.Get(3) * solution_2_1_conj);
    MC2x2F data = residual_data[vis_index];
    data -= contribution;
    const float weight = ch_block_data.Weight(vis_index);
    chi_squared += Norm(data) * weight;
  }
  return chi_squared;
}

void DiagonalAntennaSolver::PerformIteration(
    size_t ch_block, const FullSolveData::ChannelBlockData& ch_block_data,
    std::vector<MC2x2F>& residual_data, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions, size_t iteration) {
  // Fill residual_data
  std::copy(ch_block_data.DataBegin(), ch_block_data.DataEnd(),
            residual_data.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction) {
    DiagonalAddOrSubtractDirection<false>(ch_block_data, residual_data,
                                          direction, NSubSolutions(), solutions,
                                          NSubThreads());
  }

  const bool subtract_immediately = GetStepSize() > 0.99;

  const std::vector<MC2x2F> v_copy = residual_data;

  const size_t n_solve_directions = subtract_immediately
                                        ? std::min(NDirections(), iteration + 1)
                                        : NDirections();
  for (size_t direction_counter = 0; direction_counter != n_solve_directions;
       ++direction_counter) {
    const size_t direction = direction_ordering_[direction_counter];
    const uint32_t n_direction_solutions =
        ch_block_data.NSolutionsForDirection(direction);
    if (direction != 0 && !subtract_immediately) residual_data = v_copy;
    DiagonalAddOrSubtractDirection<true>(ch_block_data, residual_data,
                                         direction, NSubSolutions(), solutions,
                                         NSubThreads());

    const uint32_t solution_index0 = ch_block_data.SolutionIndex(direction, 0);
    for (uint32_t direction_solution = 0;
         direction_solution != n_direction_solutions; ++direction_solution) {
      SolveDirectionSolution(ch_block, ch_block_data, residual_data, direction,
                             solution_index0 + direction_solution, solutions,
                             next_solutions);
    }
    if (subtract_immediately) {
      std::vector<DComplex> new_solutions(next_solutions.begin(),
                                          next_solutions.end());
      DiagonalAddOrSubtractDirection<false>(ch_block_data, residual_data,
                                            direction, NSubSolutions(),
                                            new_solutions, NSubThreads());
    }
  }
}

void DiagonalAntennaSolver::SolveDirectionSolution(
    size_t ch_block, const FullSolveData::ChannelBlockData& ch_block_data,
    const std::vector<aocommon::MC2x2F>& residual_data, size_t direction_index,
    size_t solution_index, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions) {
  const size_t n_visibilities =
      ch_block_data.NSolutionVisibilities(direction_index, solution_index);
  SolverInfo solver_info{*this,         ch_block,        ch_block_data,
                         residual_data, direction_index, solution_index};

  gsl_multilarge_nlinear_fdf fdf;
  fdf.f = &PenaltyFunctionForward;
  fdf.df = &PenaltyDerivativeForward;
  fdf.fvv = &PenaltyVVForward;
  // The number of data values: 4 elements per visibility, with real and
  // imaginary values
  fdf.n = n_visibilities * 4 * 2;
  // The number of parameters: 2 elements (diagonal) with real and imaginary
  fdf.p = NAntennas() * 4;
  fdf.params = &solver_info;

  gsl_multilarge_nlinear_parameters parameters =
      gsl_multilarge_nlinear_default_parameters();

  // Enable geodesic acceleration
  parameters.trs = gsl_multilarge_nlinear_trs_lmaccel;
  parameters.avmax = 0.5;  // stability factor (0.0-1.0); lower = more stable

  gsl_multilarge_nlinear_workspace* work = gsl_multilarge_nlinear_alloc(
      gsl_multilarge_nlinear_trust, &parameters, fdf.n, fdf.p);
  gsl_vector* initial_value = gsl_vector_alloc(fdf.p);
  for (size_t antenna_index = 0; antenna_index < NAntennas(); ++antenna_index) {
    const size_t index = antenna_index * NSubSolutions() + solution_index;
    const size_t base = antenna_index * 4;

    gsl_vector_set(initial_value, base + 0, solutions[index * 2 + 0].real());
    gsl_vector_set(initial_value, base + 1, solutions[index * 2 + 0].imag());
    gsl_vector_set(initial_value, base + 2, solutions[index * 2 + 1].real());
    gsl_vector_set(initial_value, base + 3, solutions[index * 2 + 1].imag());
  }
  gsl_multilarge_nlinear_init(initial_value, &fdf, work);

  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  const size_t max_iter = 200;
  int info = 0;
  gsl_multilarge_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr,
                                &info, work);

  const gsl_vector* final_x = gsl_multilarge_nlinear_position(
      work);  // or gsl_multilarge_nlinear_x(work) in some versions

  for (size_t ant = 0; ant < NAntennas(); ++ant) {
    const size_t base = ant * 4;
    const DComplex gxx(gsl_vector_get(final_x, base + 0),
                       gsl_vector_get(final_x, base + 1));
    const DComplex gyy(gsl_vector_get(final_x, base + 2),
                       gsl_vector_get(final_x, base + 3));

    next_solutions(ch_block, ant, solution_index, 0) = gxx;
    next_solutions(ch_block, ant, solution_index, 1) = gyy;
  }

  gsl_vector_free(initial_value);
  gsl_multilarge_nlinear_free(work);
}

int DiagonalAntennaSolver::PenaltyFunction(const gsl_vector* x,
                                           SolverInfo& solver_info,
                                           gsl_vector* f) {
  // Sets f so that it contains the residual between the data and the model
  // corrected by the current parameter estimate x. f is of size
  // n_visibilities*8.
  const FullSolveData::ChannelBlockData& ch_block_data =
      solver_info.ch_block_data;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t vis_solution_index =
        ch_block_data.SolutionIndex(solver_info.direction_index, vis_index);
    if (vis_solution_index == solver_info.solution_index) {
      const MC2x2Diag solution_1 =
          GetSolution(x, ch_block_data.Antenna1Index(vis_index));
      const MC2x2Diag solution_2 =
          GetSolution(x, ch_block_data.Antenna2Index(vis_index));
      const MC2x2 data(solver_info.v_residual[vis_index]);
      const MC2x2 model(ch_block_data.ModelVisibility(
          solver_info.direction_index, vis_index));
      MC2x2 residual = solution_1 * model * solution_2.HermTranspose();
      residual -= data;
      SetVisibility(f, vis_index, residual);
    }
  }
  return 0;
}

int DiagonalAntennaSolver::PenaltyDerivative(const gsl_vector* x,
                                             const gsl_vector* u,
                                             SolverInfo& solver_info,
                                             gsl_vector* v) {
  // v should be set to J times u, where J the jacobian: J_ij = df_i / dS_j.
  // Therefore: v_i = u_{a_i} df_{ab_i} / dS_{a_i} + u_{b_i} df_{ab_i} /
  // dS_{b_i}. Here, a and b represent antenna_1 and antenna_2. On input:
  // - x has n_antennas*4 elements, giving the current solutions.
  // - u has n_antennas*4 elements.
  // On output:
  // - v has n_visibilities*8 elements.
  gsl_vector_set_all(v, 0.0);
  const FullSolveData::ChannelBlockData& ch_block_data =
      solver_info.ch_block_data;
  const size_t direction = solver_info.direction_index;
  const size_t solution_index = solver_info.solution_index;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t vis_solution_index =
        ch_block_data.SolutionIndex(direction, vis_index);
    if (vis_solution_index == solution_index) {
      const size_t antenna_1 = ch_block_data.Antenna1Index(vis_index);
      const size_t antenna_2 = ch_block_data.Antenna2Index(vis_index);
      const MC2x2Diag solution_1 = GetSolution(x, antenna_1);
      const MC2x2Diag solution_2 = GetSolution(x, antenna_2);
      const MC2x2Diag u1 = GetSolution(u, antenna_1);
      const MC2x2Diag u2 = GetSolution(u, antenna_2);
      const MC2x2 model(ch_block_data.ModelVisibility(direction, vis_index));

      // Derivative of (g1 * M * g2^H):
      //   dr/du = u1 * M * g2^H   +   g1 * M * u2^H
      const MC2x2 result = u1 * model * solution_2.HermTranspose() +
                           solution_1 * model * u2.HermTranspose();
      SetVisibility(v, vis_index, result);
    }
  }
  return 0;
}

int DiagonalAntennaSolver::PenaltyTransposedDerivative(const gsl_vector* x,
                                                       const gsl_vector* u,
                                                       SolverInfo& solver_info,
                                                       gsl_vector* v) {
  // v should be set to J^T times u, where J the jacobian: J_ij = df_i / dS_j.
  // Therefore: v_c = sum {a in ant}: u_{i_ac} df_{ab_i} / dS_{a_i} + u_{i_ca}
  // df_{ca} / dS_c, where a and c represent antenna indices and i is the
  // visibility index. On input:
  // - x has n_antennas*4 elements, giving the current solutions.
  // - u has n_visibilities*8 elements.
  // On output:
  // - v has n_antennas*4 elements.
  gsl_vector_set_all(v, 0.0);
  const FullSolveData::ChannelBlockData& ch_block_data =
      solver_info.ch_block_data;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  const size_t direction_index = solver_info.direction_index;
  const size_t solution_index = solver_info.solution_index;
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t vis_solution_index =
        ch_block_data.SolutionIndex(direction_index, vis_index);
    if (vis_solution_index != solution_index) continue;
    const size_t antenna_1 = ch_block_data.Antenna1Index(vis_index);
    const size_t antenna_2 = ch_block_data.Antenna2Index(vis_index);
    const MC2x2Diag solution_1 = GetSolution(x, antenna_1);
    const MC2x2Diag solution_2 = GetSolution(x, antenna_2);
    const MC2x2 u_value = GetVisibility(u, vis_index);
    const MC2x2 model(
        ch_block_data.ModelVisibility(direction_index, vis_index));
    const MC2x2Diag contrib1 =
        Diagonal(u_value.HermTranspose() * model * solution_2.HermTranspose());
    const MC2x2Diag contrib2 =
        Diagonal((solution_1 * model).HermTranspose() * u_value);

    AddSolution(v, antenna_1, contrib1);
    AddSolution(v, antenna_2, contrib2.Conjugate());
  }

  return 0;
}

void DiagonalAntennaSolver::PenaltyJTJ(const gsl_vector* x,
                                       SolverInfo& solver_info,
                                       gsl_matrix* jtj) {
  // - jtj has (n_antennas*4)^2 elements.
  gsl_matrix_set_all(jtj, 0.0);
  const FullSolveData::ChannelBlockData& ch_block_data =
      solver_info.ch_block_data;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  const size_t direction_index = solver_info.direction_index;
  const size_t solution_index = solver_info.solution_index;

  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t vis_solution_index =
        ch_block_data.SolutionIndex(direction_index, vis_index);
    if (vis_solution_index != solution_index) continue;

    const size_t antenna_1 = ch_block_data.Antenna1Index(vis_index);
    const size_t antenna_2 = ch_block_data.Antenna2Index(vis_index);
    const MC2x2Diag g1 = GetSolution(x, antenna_1);
    const MC2x2Diag g2 = GetSolution(x, antenna_2);
    const MC2x2 model(
        ch_block_data.ModelVisibility(direction_index, vis_index));

    const MC2x2 a = model * g2.HermTranspose();
    const MC2x2 b = g1 * model;

    // Because gains are diagonal, each visibility contributes a rank-2 update
    // to four 4×4 blocks of J^T J:
    //   - block(a1, a1)
    //   - block(a2, a2)
    //   - block(a1, a2) and block(a2, a1)  (the cross terms)

    AddSolutionToMatrix(jtj, antenna_1, antenna_1, a.HermTranspose() * a);
    AddSolutionToMatrix(jtj, antenna_2, antenna_2, b.HermTranspose() * b);

    // Cross terms:  A^H * B  and its Hermitian (B^H * A)
    const MC2x2 cross = a.HermTranspose() * b;
    AddSolutionToMatrix(jtj, antenna_1, antenna_2, cross);
    AddSolutionToMatrix(jtj, antenna_2, antenna_1, cross.HermTranspose());
  }
}

int DiagonalAntennaSolver::PenaltyVV(const gsl_vector* x, const gsl_vector* v,
                                     SolverInfo& solver_info, gsl_vector* fvv) {
  // fvv[i] = second directional derivative of the i-th residual
  //          in the direction of the velocity vector v.
  // That is:  fvv = D²f / Dv²   (where f is the residual vector)

  gsl_vector_set_all(fvv, 0.0);  // important: start from zero

  const FullSolveData::ChannelBlockData& ch_block_data =
      solver_info.ch_block_data;
  const size_t n_visibilities = ch_block_data.NVisibilities();
  const size_t direction = solver_info.direction_index;
  const size_t target_solution = solver_info.solution_index;

  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    if (ch_block_data.SolutionIndex(direction, vis_index) != target_solution)
      continue;

    const size_t a1 = ch_block_data.Antenna1Index(vis_index);
    const size_t a2 = ch_block_data.Antenna2Index(vis_index);

    const MC2x2Diag v1 = GetSolution(v, a1);  // velocity for antenna 1
    const MC2x2Diag v2 = GetSolution(v, a2);  // velocity for antenna 2

    const MC2x2 model(ch_block_data.ModelVisibility(direction, vis_index));

    // The model prediction is:   residual = g1 * model * g2^H - data
    // Because data is constant w.r.t. the parameters, its second derivative is
    // zero. So we only need the second directional derivative of (g1 * M *
    // g2^H).

    // For diagonal gains the second directional derivative simplifies:
    // Let A = g1 * M * g2^H
    // Then D_v A = v1 * M * g2^H + g1 * M * v2^H
    // D_v² A = 2 * (v1 * M * v2^H)     (the cross term; all other terms vanish)

    const MC2x2 second_deriv = (v1 * model * v2.HermTranspose()) * 2.0;

    SetVisibility(fvv, vis_index, second_deriv);
  }

  return GSL_SUCCESS;
}

}  // namespace dp3::ddecal
