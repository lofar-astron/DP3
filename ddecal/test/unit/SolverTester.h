// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLVER_TESTER_H
#define DP3_DDECAL_SOLVER_TESTER_H

#include "../../../base/DPBuffer.h"
#include "../../../ddecal/gain_solvers/BdaSolverBuffer.h"

#include <aocommon/uvector.h>

#include <complex>
#include <vector>

namespace dp3 {
namespace ddecal {

class BdaSolverBase;
class RegularSolverBase;
class SolverBase;

namespace test {

/**
 * The SolverTester is a fixture for testing solvers. It supports both
 * regular data and BDA data.
 */
class SolverTester {
 public:
  SolverTester();

  /**
   * Creates BDA data, for testing BDA solvers.
   * @return The internal solver buffer that contains the data.
   */
  const ddecal::BdaSolverBuffer& FillBDAData();

  /**
   * Initializes a solver using default values. After using this function, a
   * test can adjust the default values.
   */
  void InitializeSolver(dp3::ddecal::RegularSolverBase& solver) const;
  void InitializeSolver(dp3::ddecal::BdaSolverBase& solver) const;

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

  /**
   * @return The first antennas for the generated baselines.
   */
  const std::vector<int>& Antennas1() const { return antennas1_; }

  /**
   * @return The second antennas for the generated baselines.
   */
  const std::vector<int>& Antennas2() const { return antennas2_; }

  static constexpr size_t kNPolarizations = 4;
  static constexpr size_t kNAntennas = 50;
  static constexpr size_t kNDirections = 3;
  static constexpr size_t kNChannels = 10;
  static constexpr size_t kNChannelBlocks = 4;
  static constexpr size_t kNBDATimes = 128;
  static constexpr size_t kNBaselines = kNAntennas * (kNAntennas - 1) / 2;

  // Default solver settings:
  static constexpr size_t kMaxIterations = 100;
  static constexpr double kAccuracy = 1e-8;
  static constexpr double kStepSize = 0.2;
  static constexpr size_t kNThreads = 4;
  static constexpr bool kPhaseOnly = false;

 private:
  void InitializeSolverSettings(dp3::ddecal::SolverBase& solver) const;

  std::vector<int> antennas1_;
  std::vector<int> antennas2_;
  std::vector<std::complex<float>> input_solutions_;
  std::vector<std::vector<std::complex<double>>> solver_solutions_;

  BdaSolverBuffer bda_solver_buffer_;
};

}  // namespace test
}  // namespace ddecal
}  // namespace dp3

#endif  // DP3_DDECAL_SOLVER_TESTER_H
