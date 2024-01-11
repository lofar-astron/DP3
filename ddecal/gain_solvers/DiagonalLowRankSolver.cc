#include "DiagonalLowRankSolver.h"

#include <cassert>
#include <random>

#include <aocommon/dynamicfor.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor-blas/xlinalg.hpp>

using aocommon::MC2x2F;

namespace dp3::ddecal {

float DominantEigenPair(const xt::xtensor<std::complex<float>, 2>& matrix,
                        xt::xtensor<std::complex<float>, 1>& eigen_vector,
                        size_t n_iterations) {
  const size_t width = matrix.shape()[0];
  assert(width == matrix.shape()[1]);
  const float max_value = std::max<float>(xt::amax(xt::real(matrix))(),
                                          xt::amax(xt::imag(matrix))());
  std::mt19937 gen;
  std::uniform_real_distribution distribution(0.0, 1.0);
  eigen_vector.resize({width});
  for (std::complex<float>& v : eigen_vector) {
    v = std::complex<float>(distribution(gen) * max_value,
                            distribution(gen) * max_value);
  }
  return DominantEigenPairNear(matrix, eigen_vector, n_iterations);
}

float DominantEigenPairNear(const xt::xtensor<std::complex<float>, 2>& matrix,
                            xt::xtensor<std::complex<float>, 1>& eigen_vector,
                            size_t n_iterations) {
  for (size_t iteration = 0; iteration != n_iterations; ++iteration) {
    eigen_vector = xt::linalg::dot(matrix, xt::transpose(eigen_vector));
    eigen_vector = eigen_vector / xt::linalg::norm(eigen_vector, 2);

    std::real(xt::linalg::dot(xt::linalg::dot(xt::conj(eigen_vector), matrix),
                              eigen_vector)[0]);
  }

  // Calculate eigen value using Rayleigh's quotient. Note that the denominator
  // is already normalized to 1, so the eigen value is just the numerator of the
  // quotient.
  const float eigen_value = std::real(xt::linalg::dot(
      xt::linalg::dot(xt::conj(eigen_vector), matrix), eigen_vector)[0]);
  return eigen_value;
}

DiagonalLowRankSolver::SolveResult DiagonalLowRankSolver::Solve(
    const SolveData& data, std::vector<std::vector<DComplex>>& solutions,
    double time, std::ostream* stat_stream) {
  PrepareConstraints();

  const bool subtract_immediately = GetStepSize() > 0.99;
  if (subtract_immediately) {
    CalculateNormPerDirection(data);
  } else {
    direction_ordering_.reserve(NDirections());
    for (size_t i = 0; i != NDirections(); ++i)
      direction_ordering_.emplace_back(i);
  }

  SolutionTensor next_solutions(
      {NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()});

  SolveResult result;

  // Visibility vector v_residual[cb][vis] of size NChannelBlocks() x
  // n_visibilities
  std::vector<std::vector<MC2x2F>> v_residual(NChannelBlocks());
  // Allocates all structures
  for (size_t ch_block = 0; ch_block != NChannelBlocks(); ++ch_block) {
    v_residual[ch_block].resize(data.ChannelBlock(ch_block).NVisibilities());
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

    aocommon::DynamicFor<size_t> loop;
    loop.Run(0, NChannelBlocks(), [&](size_t ch_block) {
      PerformIteration(ch_block, data.ChannelBlock(ch_block),
                       v_residual[ch_block], solutions[ch_block],
                       next_solutions, iteration);
    });

    Step(solutions, next_solutions);

    bool constraints_satisfied =
        ApplyConstraints(iteration, time, has_previously_converged, result,
                         next_solutions, stat_stream);

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

void DiagonalLowRankSolver::CalculateNormPerDirection(const SolveData& data) {
  std::vector<std::pair<float, size_t>> norm_sum_per_direction(NDirections());
  for (size_t channel_block_index = 0;
       channel_block_index != data.NChannelBlocks(); ++channel_block_index) {
    const SolveData::ChannelBlockData& cb_data =
        data.ChannelBlock(channel_block_index);
    const size_t n_visibilities = cb_data.NVisibilities();
    for (size_t direction = 0; direction != cb_data.NDirections();
         ++direction) {
      norm_sum_per_direction[direction].second = direction;
      for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
        const aocommon::MC2x2F& model =
            cb_data.ModelVisibility(direction, vis_index);
        if (cb_data.Weight(vis_index) != 0.0)
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

double DiagonalLowRankSolver::ChiSquared(
    const SolveData::ChannelBlockData& cb_data,
    std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
    const SolutionSpan& solutions) const {
  using DComplex = std::complex<double>;
  using Complex = std::complex<float>;
  using aocommon::MC2x2F;
  double chi_squared = 0.0;
  const size_t n_visibilities = cb_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
    const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
    const uint32_t solution_index = cb_data.SolutionIndex(direction, vis_index);
    const DComplex* solution_1 =
        &solutions[(antenna_1 * NSolutions() + solution_index) * 2];
    const DComplex* solution_2 =
        &solutions[(antenna_2 * NSolutions() + solution_index) * 2];
    const Complex solution_1_0(solution_1[0]);
    const Complex solution_1_1(solution_1[1]);
    const Complex solution_2_0_conj(std::conj(solution_2[0]));
    const Complex solution_2_1_conj(std::conj(solution_2[1]));

    const MC2x2F& model = cb_data.ModelVisibility(direction, vis_index);
    const MC2x2F contribution(solution_1_0 * model[0] * solution_2_0_conj,
                              solution_1_0 * model[1] * solution_2_1_conj,
                              solution_1_1 * model[2] * solution_2_0_conj,
                              solution_1_1 * model[3] * solution_2_1_conj);
    MC2x2F data = v_residual[vis_index];
    data -= contribution;
    const float weight = cb_data.Weight(vis_index);
    chi_squared += Norm(data) * weight;
  }
  return chi_squared;
}

void DiagonalLowRankSolver::PerformIteration(
    size_t ch_block, const SolveData::ChannelBlockData& cb_data,
    std::vector<MC2x2F>& v_residual, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions, size_t iteration) {
  // Fill v_residual
  std::copy(cb_data.DataBegin(), cb_data.DataEnd(), v_residual.begin());

  // Subtract all directions with their current solutions
  for (size_t direction = 0; direction != NDirections(); ++direction) {
    DiagonalAddOrSubtractDirection<false>(cb_data, v_residual, direction,
                                          NSolutions(), solutions);
  }

  const bool subtract_immediately = GetStepSize() > 0.99;

  const std::vector<MC2x2F> v_copy = v_residual;

  const size_t n_solve_directions = subtract_immediately
                                        ? std::min(NDirections(), iteration + 1)
                                        : NDirections();
  for (size_t direction_counter = 0; direction_counter != n_solve_directions;
       ++direction_counter) {
    const size_t direction = direction_ordering_[direction_counter];
    const uint32_t n_direction_solutions =
        cb_data.NSolutionsForDirection(direction);
    if (direction != 0 && !subtract_immediately) v_residual = v_copy;
    DiagonalAddOrSubtractDirection<true>(cb_data, v_residual, direction,
                                         NSolutions(), solutions);

    const uint32_t solution_index0 = cb_data.SolutionIndex(direction, 0);
    for (uint32_t direction_solution = 0;
         direction_solution != n_direction_solutions; ++direction_solution) {
      SolveDirectionSolution(ch_block, cb_data, v_residual, direction,
                             solution_index0 + direction_solution, solutions,
                             next_solutions);
    }
    if (subtract_immediately) {
      std::vector<DComplex> new_solutions(next_solutions.begin(),
                                          next_solutions.end());
      DiagonalAddOrSubtractDirection<false>(cb_data, v_residual, direction,
                                            NSolutions(), new_solutions);
    }
  }
}

void AddToCorrelation(std::complex<float>& correlation_element,
                      float& variance_element, float& divisor_element,
                      const std::complex<float>& data, float weight,
                      const std::complex<float>& model) {
  if (std::isfinite(data.real()) && std::isfinite(data.imag())) {
    const std::complex<float> division = data / model;
    // VAR(division) = VAR(data / model) = VAR(data) / model^2
    const float inv_divisor_variance = std::norm(model);
    // Perform inverse-variance weighting, hence multiply with 1/model^2
    correlation_element += weight * inv_divisor_variance * division;
    // VAR(weight * division / model^2) = VAR(data) * weight^2 / norm(model)^2 *
    // norm(model)^2
    variance_element +=
        weight * weight * inv_divisor_variance * inv_divisor_variance;
    divisor_element += weight * inv_divisor_variance;
  }
}

void DiagonalLowRankSolver::SolveDirectionSolution(
    size_t ch_block, const SolveData::ChannelBlockData& cb_data,
    const std::vector<aocommon::MC2x2F>& v_residual, size_t direction_index,
    size_t solution_index, const std::vector<DComplex>& solutions,
    SolutionTensor& next_solutions) {
  const std::array<size_t, 2> shape{NAntennas() * 2, NAntennas() * 2};
  xt::xtensor<std::complex<float>, 2> correlation_matrix(shape);
  correlation_matrix.fill(0.0);
  xt::xtensor<float, 2> variance_matrix(shape);
  variance_matrix.fill(0.0);
  xt::xtensor<float, 2> divisor_sum(shape);
  divisor_sum.fill(0.0);
  // Iterate over all data
  const size_t n_visibilities = cb_data.NVisibilities();
  for (size_t vis_index = 0; vis_index != n_visibilities; ++vis_index) {
    const uint32_t vis_solution_index =
        cb_data.SolutionIndex(direction_index, vis_index);
    if (vis_solution_index == solution_index) {
      const uint32_t antenna_1 = cb_data.Antenna1Index(vis_index);
      const uint32_t antenna_2 = cb_data.Antenna2Index(vis_index);
      const MC2x2F& data = v_residual[vis_index];
      const float* weight = &cb_data.Weight(vis_index);
      const MC2x2F& model = cb_data.ModelVisibility(direction_index, vis_index);

      AddToCorrelation(correlation_matrix(antenna_1 * 2, antenna_2 * 2),
                       variance_matrix(antenna_1 * 2, antenna_2 * 2),
                       divisor_sum(antenna_1 * 2, antenna_2 * 2), data[0],
                       weight[0], model[0]);
      AddToCorrelation(correlation_matrix(antenna_1 * 2, antenna_2 * 2 + 1),
                       variance_matrix(antenna_1 * 2, antenna_2 * 2 + 1),
                       divisor_sum(antenna_1 * 2, antenna_2 * 2 + 1), data[1],
                       weight[1], model[1]);
      AddToCorrelation(correlation_matrix(antenna_1 * 2 + 1, antenna_2 * 2),
                       variance_matrix(antenna_1 * 2 + 1, antenna_2 * 2),
                       divisor_sum(antenna_1 * 2 + 1, antenna_2 * 2), data[2],
                       weight[2], model[2]);
      AddToCorrelation(correlation_matrix(antenna_1 * 2 + 1, antenna_2 * 2 + 1),
                       variance_matrix(antenna_1 * 2 + 1, antenna_2 * 2 + 1),
                       divisor_sum(antenna_1 * 2 + 1, antenna_2 * 2 + 1),
                       data[3], weight[3], model[3]);

      AddToCorrelation(correlation_matrix(antenna_2 * 2, antenna_1 * 2),
                       variance_matrix(antenna_2 * 2, antenna_1 * 2),
                       divisor_sum(antenna_2 * 2, antenna_1 * 2),
                       std::conj(data[0]), weight[0], std::conj(model[0]));
      AddToCorrelation(correlation_matrix(antenna_2 * 2, antenna_1 * 2 + 1),
                       variance_matrix(antenna_2 * 2, antenna_1 * 2 + 1),
                       divisor_sum(antenna_2 * 2, antenna_1 * 2 + 1),
                       std::conj(data[2]), weight[2], std::conj(model[2]));
      AddToCorrelation(correlation_matrix(antenna_2 * 2 + 1, antenna_1 * 2),
                       variance_matrix(antenna_2 * 2 + 1, antenna_1 * 2),
                       divisor_sum(antenna_2 * 2 + 1, antenna_1 * 2),
                       std::conj(data[1]), weight[1], std::conj(model[1]));
      AddToCorrelation(correlation_matrix(antenna_2 * 2 + 1, antenna_1 * 2 + 1),
                       variance_matrix(antenna_2 * 2 + 1, antenna_1 * 2 + 1),
                       divisor_sum(antenna_2 * 2 + 1, antenna_1 * 2 + 1),
                       std::conj(data[3]), weight[3], std::conj(model[3]));
    }
  }

  for (size_t i = 0; i != correlation_matrix.size(); ++i) {
    if (divisor_sum[i] == 0.0) {
      correlation_matrix[i] = 0.0;
      variance_matrix[i] = 0.0;
    } else {
      const float divisor = divisor_sum[i];
      correlation_matrix[i] /= divisor;
      variance_matrix[i] /= divisor * divisor;
    }
  }
  const float max_weight = xt::amax(variance_matrix)[0];
  if (max_weight != 0.0) variance_matrix = variance_matrix / max_weight;

  xt::xtensor<std::complex<float>, 1> x_t;
  x_t.resize({NAntennas() * 2});
  for (size_t antenna = 0; antenna != NAntennas(); ++antenna) {
    x_t[antenna * 2] = solutions[(antenna * NSolutions() + solution_index) * 2];
    x_t[antenna * 2 + 1] =
        solutions[(antenna * NSolutions() + solution_index) * 2 + 1];
  }
  // Compute:  (from: Srebro & Jaakkola, 2003)
  // X_{t+1} = LRA_k ( weights * A + (1 - weights) * X_t )
  // Perform the initial iteration of above formula. This iteration is
  // different, as we do not use an estimate for the eigen value yet.
  float eigen_value = 1.0;
  xt::xtensor<std::complex<float>, 2> term =
      variance_matrix * correlation_matrix +
      (1.0f - variance_matrix) *
          (eigen_value * xt::linalg::outer(x_t, xt::conj(x_t)));
  // Because we start with no estimate, we iterate 2x more during the
  // first power method.
  eigen_value = DominantEigenPair(term, x_t, n_power_iterations_ * 2);

  // Continue iterating the same formula
  for (size_t i = 1; i < n_low_rank_approximation_iterations_; ++i) {
    xt::noalias(term) =
        variance_matrix * correlation_matrix +
        (1.0f - variance_matrix) *
            (eigen_value * xt::linalg::outer(x_t, xt::conj(x_t)));
    eigen_value = DominantEigenPairNear(term, x_t, n_power_iterations_);
  }
  // A complex sqrt is performed to allow for negative eigen values
  // TODO find way to do this without extra copy
  xt::xarray<std::complex<float>> copy =
      x_t * std::sqrt(std::complex<float>(eigen_value));
  std::complex<float> phase_reference = copy[0] / std::abs(copy[0]);
  copy = copy / phase_reference;
  copy.reshape({NAntennas(), 2});
  xt::view(next_solutions, ch_block, xt::all(), solution_index, xt::all()) =
      xt::cast<std::complex<double>>(copy);
}

}  // namespace dp3::ddecal
