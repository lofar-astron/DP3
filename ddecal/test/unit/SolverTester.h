// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLVER_TESTER_H
#define DP3_DDECAL_SOLVER_TESTER_H

#include "../../../base/DPBuffer.h"
#include "../../../ddecal/gain_solvers/BDASolverBuffer.h"
#include "../../../ddecal/gain_solvers/SolverBuffer.h"

#include <aocommon/uvector.h>

#include <complex>
#include <vector>

namespace dp3 {
namespace base {
class RegularSolverBase;
class BdaSolverBase;

namespace test {

/**
 * The SolverTester is a fixture for testing solvers. It supports both
 * regular data and BDA data.
 */
class SolverTester {
 public:
  SolverTester();

  /**
   * Creates regular data, for testing regular solvers.
   * @return The internal solver buffer that contains the data.
   */
  const SolverBuffer& FillData();

  /**
   * Creates BDA data, for testing BDA solvers.
   * @return The internal solver buffer that contains the data.
   */
  const BDASolverBuffer& FillBDAData();

  /**
   * Initializes a solver using default values. After using this function, a
   * test can adjust the default values.
   */
  void InitializeSolver(dp3::base::RegularSolverBase& solver) const;
  void InitializeSolver(dp3::base::BdaSolverBase& solver) const;

  /**
   * Initializes input solutions and solver solutions.
   * The input solutions are for generating data.
   * Sets the solver solutions to unit values with the appropriate dimensions.
   * The solver can then use these values as initial values.
   * @{
   */
  void SetScalarSolutions();
  void SetDiagonalSolutions();
  /** @} */

  /**
   * Get the solution data that can be passed to the solver as initial values.
   * and that will contain the solve result after running the solver.
   * The various Check functions like @ref CheckScalarResults() can verify
   * the solve result.
   * @see SetScalarSolutions() et al.
   */
  std::vector<std::vector<std::complex<double>>>& GetSolverSolutions() {
    return solver_solutions_;
  }

  /**
   * Checks if the solver solutions match the generated input solutions.
   * @param tolerance Tolerance value for BOOST_CHECK_CLOSE.
   */
  void CheckScalarResults(double tolerance);
  void CheckDiagonalResults(double tolerance);

  const size_t kNPolarizations = 4;
  const size_t kNAntennas = 50;
  const size_t kNDirections = 3;
  const size_t kNChannels = 10;
  const size_t kNChannelBlocks = 4;
  const size_t kNRegularTimes = 50;
  // Use more times with BDA so the number of visibilities is similar to the
  // regular data.
  const size_t kNBDATimes = 128;
  const size_t kNBaselines = kNAntennas * (kNAntennas - 1) / 2;
  // For each block of 8 baselines, the averaged data size corresponds to
  // 1 + 1/2 + 1/2 + 1/3 + 1/3 + 1/6 + 1/6 + 1/16 = 3 + 1/16
  // non-averaged baselines. -> Divide the non-averaged data by two.
  const size_t kBDABufferSize =
      kNBDATimes * kNBaselines * kNChannels * kNPolarizations / 2;

  // Default solver settings:
  const size_t kMaxIterations = 100;
  const double kAccuracy = 1e-8;
  const double kStepSize = 0.2;
  const size_t kNThreads = 4;
  const bool kPhaseOnly = false;

  const std::vector<int>& Antennas1() const { return antennas1_; }
  const std::vector<int>& Antennas2() const { return antennas2_; }

 private:
  std::vector<int> antennas1_;
  std::vector<int> antennas2_;
  std::vector<std::complex<float>> input_solutions_;
  std::vector<std::vector<std::complex<double>>> solver_solutions_;

  std::vector<DPBuffer> data_buffers_;
  std::vector<std::vector<DPBuffer>> model_buffer_store_;
  std::vector<aocommon::UVector<std::complex<float>>> data_store_;
  std::vector<aocommon::UVector<float>> weight_store_;
  SolverBuffer solver_buffer_;

  std::unique_ptr<BDABuffer> bda_data_buffer_;
  std::vector<std::unique_ptr<BDABuffer>> bda_model_buffers_;
  BDASolverBuffer bda_solver_buffer_;
};

}  // namespace test
}  // namespace base
}  // namespace dp3

#endif  // DP3_DDECAL_SOLVER_TESTER_H
