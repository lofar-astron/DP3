// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_DIAGONAL_LOW_RANK_SOLVER_H_
#define DDECAL_DIAGONAL_LOW_RANK_SOLVER_H_

#include "SolverBase.h"
#include "SolveData.h"

namespace dp3::ddecal {

/**
 * Calculates the largest eigenvalue and corresponding eigen vector using the
 * power iteration method. The eigen vector is calculated from the eigen vector
 * using the Rayleigh quotient.
 * @param [out] eigen_vector Used to store the calculated eigen vector.
 * @param n_iterations Number of power iterations. Typical values are 10-20.
 * @returns Largest eigen value of the given matrix.
 */
float DominantEigenPair(const xt::xtensor<std::complex<float>, 2>& matrix,
                        xt::xtensor<std::complex<float>, 1>& eigen_vector,
                        size_t n_iterations);

/**
 * Like @ref DominantEigenPair, but starts from an estimate of the eigen vector
 * to speed up the calculation.
 * @param [in,out] eigen_vector on input, an estimate of the eigen vector
 * corresponding to the dominant eigen value. On output, the calculated eigen
 * vector.
 * @param n_iterations Number of power iterations. Typical values are 10-20.
 * @returns Largest eigen value of the given matrix.
 */
float DominantEigenPairNear(const xt::xtensor<std::complex<float>, 2>& matrix,
                            xt::xtensor<std::complex<float>, 1>& eigen_vector,
                            size_t n_iterations);

/**
 * Experimental solver that make use of the property that the antenna x antenna
 * visibility matrix is a low-rank (rank 1 to be precise) matrix. A normal
 * low-rank approximation can be solved exactly and without iterations, but
 * because there are weights and missing values in the matrix, an iterative
 * procedure is required, and the result is not exact.
 *
 * When the step size is set to 1, this solver makes an immediate step. Because
 * a solving iteration uses all data for one direction at once, this may be a
 * good approach (in contrast to the iterative diagonal solver). To help
 * convergence in this case, the directions are solved in order of largest power
 * first.
 */
class DiagonalLowRankSolver final : public SolverBase {
 public:
  DiagonalLowRankSolver(size_t n_low_rank_approximation_iterations,
                        size_t n_power_iterations)
      : SolverBase(),
        n_low_rank_approximation_iterations_(
            n_low_rank_approximation_iterations),
        n_power_iterations_(n_power_iterations) {}

  SolveResult Solve(const SolveData& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 2; }

  /**
   * Number of power iterations performed to calculate the eigen value.
   * See also @ref DominantEigenPair().
   */
  void SetNPowerIterations(size_t n_power_iterations) {
    n_power_iterations_ = n_power_iterations;
  }

  /**
   * Number of full low-rank approximation iterations.
   */
  void SetNLowRankApproximationIterations(size_t n_iterations) {
    n_low_rank_approximation_iterations_ = n_iterations;
  }

 private:
  void PerformIteration(size_t ch_block,
                        const SolveData::ChannelBlockData& cb_data,
                        std::vector<aocommon::MC2x2F>& v_residual,
                        const std::vector<DComplex>& solutions,
                        SolutionTensor& next_solutions, size_t iteration);

  /**
   * Calculates the chi-squared value, i.e. the weighted squared sum of the
   * difference between the data and the corrected residual.
   */
  double ChiSquared(const SolveData::ChannelBlockData& cb_data,
                    std::vector<aocommon::MC2x2F>& v_residual, size_t direction,
                    const SolutionSpan& solutions) const;

  void SolveDirectionSolution(size_t ch_block,
                              const SolveData::ChannelBlockData& cb_data,
                              const std::vector<aocommon::MC2x2F>& v_residual,
                              size_t direction_index, size_t solution_index,
                              const std::vector<DComplex>& solutions,
                              SolutionTensor& next_solutions);

  void CalculateNormPerDirection(const SolveData& data);

  std::vector<size_t> direction_ordering_;
  size_t n_low_rank_approximation_iterations_ = 25;
  size_t n_power_iterations_ = 10;
};

}  // namespace dp3::ddecal

#endif  // DDECAL_DIAGONAL_LOW_RANK_SOLVER_H_
