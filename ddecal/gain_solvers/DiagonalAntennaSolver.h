// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_DIAGONAL_ANTENNA_SOLVER_H_
#define DDECAL_DIAGONAL_ANTENNA_SOLVER_H_

#include "SolverBase.h"
#include "SolveData.h"

#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace dp3::ddecal {

class DiagonalAntennaSolver final : public SolverBase {
 public:
  SolveResult Solve(const FullSolveData& data,
                    std::vector<std::vector<DComplex>>& solutions,
                    double time) override;

  size_t NSolutionPolarizations() const override { return 2; }

 private:
  void PerformIteration(size_t ch_block,
                        const FullSolveData::ChannelBlockData& ch_block_data,
                        std::vector<aocommon::MC2x2F>& residual_data,
                        const std::vector<DComplex>& solutions,
                        SolutionTensor& next_solutions, size_t iteration);

  /**
   * Calculates the chi-squared value, i.e. the weighted squared sum of the
   * difference between the data and the corrected residual.
   */
  double ChiSquared(const FullSolveData::ChannelBlockData& ch_block_data,
                    std::vector<aocommon::MC2x2F>& residual_data,
                    size_t direction, const SolutionSpan& solutions) const;

  void SolveDirectionSolution(
      size_t ch_block, const FullSolveData::ChannelBlockData& ch_block_data,
      const std::vector<aocommon::MC2x2F>& residual_data,
      size_t direction_index, size_t solution_index,
      const std::vector<DComplex>& solutions, SolutionTensor& next_solutions);

  void CalculateNormPerDirection(const FullSolveData& data);

  struct SolverInfo {
    DiagonalAntennaSolver& solver;
    size_t ch_block;
    const FullSolveData::ChannelBlockData& ch_block_data;
    const std::vector<aocommon::MC2x2F>& v_residual;
    size_t direction_index;
    size_t solution_index;
  };

  /**
   * @param [in] x The solutions
   * @param [out] f The differences between corrected model data and observed
   * data.
   */
  int PenaltyFunction(const gsl_vector* x, SolverInfo& solver_info,
                      gsl_vector* f);
  static int PenaltyFunctionForward(const gsl_vector* x, void* solver_info,
                                    gsl_vector* f);

  /**
   * Calculates the Jacobian matrix times a vector u.
   * @param [in] x The solutions
   * @param [in] u The values to be multiplied with the Jacobian
   * @param [out] v The result
   */
  int PenaltyDerivative(const gsl_vector* x, const gsl_vector* u,
                        SolverInfo& solver_info, gsl_vector* v);
  /**
   * Calculates the transposed Jacobian times a vector u.
   * @param [in] x The solutions
   * @param [in] u The values to be multiplied with the Jacobian
   * @param [out] v The result
   */
  int PenaltyTransposedDerivative(const gsl_vector* x, const gsl_vector* u,
                                  SolverInfo& solver_info, gsl_vector* v);
  /**
   * Calculates J^T J.
   * @param [in] x The solutions
   * @param [out] jtj The result
   */
  void PenaltyJTJ(const gsl_vector* x, SolverInfo& solver_info,
                  gsl_matrix* jtj);

  static int PenaltyDerivativeForward(CBLAS_TRANSPOSE_t transpose_j,
                                      const gsl_vector* x, const gsl_vector* u,
                                      void* solver_info, gsl_vector* v,
                                      gsl_matrix* jtj);

  /**
   * Calculate second directional derivative in direction v.
   * @param [in] x The solutions
   * @param [out] fvv The result
   */
  int PenaltyVV(const gsl_vector* x, const gsl_vector* v,
                SolverInfo& solver_info, gsl_vector* fvv);
  static int PenaltyVVForward(const gsl_vector* x, const gsl_vector* v,
                              void* solver_info, gsl_vector* fvv);

  std::vector<size_t> direction_ordering_;
};

}  // namespace dp3::ddecal

#endif  // DDECAL_DIAGONAL_ANTENNA_SOLVER_H_
