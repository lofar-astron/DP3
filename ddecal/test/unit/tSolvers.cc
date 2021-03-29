// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <vector>
#include <complex>

#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

#include "../../gain_solvers/DiagonalSolver.h"
#include "../../gain_solvers/FullJonesSolver.h"
#include "../../gain_solvers/IterativeDiagonalSolver.h"
#include "../../gain_solvers/IterativeScalarSolver.h"
#include "../../gain_solvers/ScalarSolver.h"

#include "../../linear_solvers/LSMRSolver.h"
#include "../../linear_solvers/NormalEquationsSolver.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using aocommon::MC2x2;
using namespace dp3::base;

class SolverTester {
 public:
  SolverTester()
      : input_solutions(n_ant * n_dir * 2),
        model_data_store(n_times),
        data_store(n_times),
        weight_store(n_times) {
    for (size_t a1 = 0; a1 != n_ant; ++a1) {
      for (size_t a2 = a1 + 1; a2 != n_ant; ++a2) {
        ant1s.push_back(a1);
        ant2s.push_back(a2);
      }
    }
  }

  void FillData() {
    std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
    std::mt19937 mt(0);
    for (size_t timestep = 0; timestep != n_times; ++timestep) {
      aocommon::UVector<cf>& time_data = data_store[timestep];
      time_data.assign(n_pol * n_chan * nBl, 0);
      aocommon::UVector<float>& time_weights = weight_store[timestep];
      time_weights.assign(n_pol * n_chan * nBl, 0);
      // modelTimeData has dimensions [direction][pol x ch x bl]
      std::vector<aocommon::UVector<cf>>& model_time_data =
          model_data_store[timestep];

      for (size_t d = 0; d != n_dir; ++d) {
        model_time_data.emplace_back(n_pol * n_chan * nBl);
        aocommon::UVector<cf>& this_direction = model_time_data[d];

        for (size_t bl = 0; bl != nBl; ++bl) {
          for (size_t ch = 0; ch != n_chan; ++ch) {
            const size_t matrix_index = (bl * n_chan + ch) * 4;
            this_direction[matrix_index + 0] =
                cf(uniform_data(mt), uniform_data(mt));
            this_direction[matrix_index + 1] =
                cf(uniform_data(mt) * 0.1, uniform_data(mt) * 0.1);
            this_direction[matrix_index + 2] =
                cf(uniform_data(mt) * 0.1, uniform_data(mt) * 0.1);
            this_direction[matrix_index + 3] =
                cf(1.5 * uniform_data(mt), 1.5 * uniform_data(mt));
          }
        }
      }
      size_t baseline_index = 0;
      for (size_t a1 = 0; a1 != n_ant; ++a1) {
        for (size_t a2 = a1 + 1; a2 != n_ant; ++a2) {
          for (size_t ch = 0; ch != n_chan; ++ch) {
            MC2x2 perturbed_model = MC2x2::Zero();
            for (size_t d = 0; d != n_dir; ++d) {
              MC2x2 val(
                  &model_time_data[d][(baseline_index * n_chan + ch) * 4]);
              MC2x2 left(input_solutions[(a1 * n_dir + d) * 2 + 0], 0.0, 0.0,
                         input_solutions[(a1 * n_dir + d) * 2 + 1]);
              MC2x2 right(input_solutions[(a2 * n_dir + d) * 2 + 0], 0.0, 0.0,
                          input_solutions[(a2 * n_dir + d) * 2 + 1]);
              MC2x2 res;
              MC2x2::ATimesB(res, left, val);
              MC2x2::ATimesHermB(left, res,
                                 right);  // left is used as scratch for result
              perturbed_model += left;
            }
            for (size_t p = 0; p != 4; ++p) {
              time_data[(baseline_index * n_chan + ch) * 4 + p] =
                  perturbed_model[p];
              time_weights[(baseline_index * n_chan + ch) * 4 + p] = 1.0;
            }
          }
          ++baseline_index;
        }
      }

      data.push_back(time_data.data());
      weights.push_back(time_weights.data());
      std::vector<cf*> model_ptrs;
      model_ptrs.reserve(model_data.size());
      for (aocommon::UVector<cf>& model_dir : model_time_data)
        model_ptrs.push_back(model_dir.data());
      model_data.push_back(model_ptrs);
    }
  }

  void CheckScalarResults(
      const std::vector<std::vector<std::complex<double>>>& solutions,
      double tolerance) {
    for (size_t ch = 0; ch != n_chan_blocks; ++ch) {
      for (size_t ant = 0; ant != n_ant; ++ant) {
        for (size_t d = 0; d != n_dir; ++d) {
          const std::complex<double> sol0 = solutions[ch][d];
          const std::complex<double> inp0 = input_solutions[d * 2];

          const std::complex<double> sol = solutions[ch][d + ant * n_dir];
          const std::complex<double> inp =
              input_solutions[(d + ant * n_dir) * 2];

          // Compare the squared quantities, because the phase has an ambiguity
          BOOST_CHECK_CLOSE(std::norm(sol), std::norm(inp), tolerance);

          // Reference to antenna0 to check if relative phase is correct
          BOOST_CHECK_LT(std::abs((sol * std::conj(sol0)).real() -
                                  (inp * std::conj(inp0)).real()),
                         tolerance);
          BOOST_CHECK_LT(std::abs((sol * std::conj(sol0)).imag() -
                                  (inp * std::conj(inp0)).imag()),
                         tolerance);
        }
      }
    }
  }

  void CheckDiagonalResults(
      const std::vector<std::vector<std::complex<double>>>& solutions,
      double tolerance) {
    for (size_t ch = 0; ch != n_chan_blocks; ++ch) {
      for (size_t ant = 0; ant != n_ant; ++ant) {
        for (size_t d = 0; d != n_dir; ++d) {
          std::complex<double> solX0 = solutions[ch][d * 2];
          std::complex<double> solY0 = solutions[ch][d * 2 + 1];
          std::complex<double> inpX0 = input_solutions[d * 2];
          std::complex<double> inpY0 = input_solutions[d * 2 + 1];

          std::complex<double> solX = solutions[ch][(d + ant * n_dir) * 2];
          std::complex<double> solY = solutions[ch][(d + ant * n_dir) * 2 + 1];
          std::complex<double> inpX = input_solutions[(d + ant * n_dir) * 2];
          std::complex<double> inpY =
              input_solutions[(d + ant * n_dir) * 2 + 1];

          // Compare the squared quantities, because the phase has an ambiguity
          BOOST_CHECK_CLOSE(std::norm(solX), std::norm(inpX), tolerance);
          BOOST_CHECK_CLOSE(std::norm(solY), std::norm(inpY), tolerance);

          // Reference to antenna0 to check if relative phase is correct
          BOOST_CHECK_CLOSE((solX * std::conj(solX0)).real(),
                            (inpX * std::conj(inpX0)).real(), tolerance);
          BOOST_CHECK_CLOSE((solY * std::conj(solY0)).real(),
                            (inpY * std::conj(inpY0)).real(), tolerance);
        }
      }
    }
  }

  void SetScalarSolutions() {
    std::mt19937 mt;
    std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
    for (size_t a = 0; a != n_ant; ++a) {
      for (size_t d = 0; d != n_dir; ++d) {
        if (d == 0)
          input_solutions[(a * n_dir + d) * 2] =
              cf(uniform_sols(mt), uniform_sols(mt));
        else
          input_solutions[(a * n_dir + d) * 2] =
              cf(uniform_sols(mt) * 0.5, uniform_sols(mt) * 0.5);
        input_solutions[(a * n_dir + d) * 2 + 1] =
            input_solutions[(a * n_dir + d) * 2];
      }
    }
  }

  typedef std::complex<float> cf;
  size_t n_pol = 4, n_ant = 50, n_dir = 3, n_chan = 10, n_chan_blocks = 4,
         n_times = 50, nBl = n_ant * (n_ant - 1) / 2, max_iter = 100;
  std::vector<int> ant1s, ant2s;
  std::vector<cf> input_solutions;
  std::vector<cf*> data;
  std::vector<float*> weights;
  std::vector<std::vector<cf*>> model_data;
  std::vector<std::vector<aocommon::UVector<cf>>> model_data_store;
  std::vector<aocommon::UVector<cf>> data_store;
  std::vector<aocommon::UVector<float>> weight_store;
};

BOOST_AUTO_TEST_SUITE(solvers)

#ifdef USE_LSMR
BOOST_AUTO_TEST_CASE(lsmr_solver) {
  int n = 2;
  int m = 4;
  std::vector<std::complex<float>> A{1.0, 1.5, 3.5, 2.0,
                                     0.5, 0.7, 0.8, 0.4};  // column major
  std::vector<std::complex<float>> b{1.0, 2.0, 1.5, 1.2};

  LSMRSolver solver(m, n, 1);
  solver.Solve(A.data(), b.data());

  BOOST_CHECK_CLOSE(b[0].real(), -0.14141126, 1.0E-3);
  BOOST_CHECK_CLOSE(b[1].real(), 2.79757662, 1.0E-3);
}

BOOST_FIXTURE_TEST_CASE(scalar_solver_lsmr, SolverTester) {
  ScalarSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-9);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.SetLLSSolverType(LLSSolverType::LSMR,
                          std::pair<double, double>(1.0E-2, 1.0E-2));
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  SetScalarSolutions();

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-values as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckScalarResults(solutions, 1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}
#endif

BOOST_AUTO_TEST_CASE(normaleq_solver) {
  int n = 2;
  int m = 4;
  std::vector<std::complex<float>> A{1.0, 1.5, 3.5, 2.0,
                                     0.5, 0.7, 0.8, 0.4};  // column major
  std::vector<std::complex<float>> b{1.0, 2.0, 1.5, 1.2};

  NormalEquationsSolver solver(m, n, 1);
  solver.Solve(A.data(), b.data());

  BOOST_CHECK_CLOSE(b[0].real(), -0.14141126, 1.0E-3);
  BOOST_CHECK_CLOSE(b[1].real(), 2.79757662, 1.0E-3);
}

BOOST_FIXTURE_TEST_CASE(scalar_solver, SolverTester) {
  ScalarSolver solver;
  solver.SetLLSSolverType(LLSSolverType::QR,
                          std::pair<double, double>(0.0, 0.0));
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  SetScalarSolutions();

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-values as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckScalarResults(solutions, 1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_scalar_solver, SolverTester) {
  IterativeScalarSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(1);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  SetScalarSolutions();

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-values as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckScalarResults(solutions, 1.0e-3);
}

BOOST_FIXTURE_TEST_CASE(scalar_solver_normaleq, SolverTester) {
  ScalarSolver solver;
  solver.SetLLSSolverType(LLSSolverType::NORMAL_EQUATIONS,
                          std::pair<double, double>(0.0, 0.0));
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  SetScalarSolutions();

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-values as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckScalarResults(solutions, 1.0E-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}

#ifdef USE_LSMR
BOOST_FIXTURE_TEST_CASE(diagonal_solver_lsmr, SolverTester) {
  typedef std::complex<float> cf;
  DiagonalSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.SetLLSSolverType(LLSSolverType::LSMR,
                          std::pair<double, double>(1.0E-7, 1.0E-2));
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t a = 0; a != n_ant; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != n_dir; ++d) {
        if (d == 0)
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt), uniform_sols(mt));
        else
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt) * 0.5, uniform_sols(mt) * 0.5);
      }
    }
  }

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant * 2, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckDiagonalResults(solutions, 2e-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}
#endif

BOOST_FIXTURE_TEST_CASE(diagonal_solver, SolverTester) {
  typedef std::complex<float> cf;
  DiagonalSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t a = 0; a != n_ant; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != n_dir; ++d) {
        if (d == 0)
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt), uniform_sols(mt));
        else
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt) * 0.5, uniform_sols(mt) * 0.5);
      }
    }
  }

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant * 2, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckDiagonalResults(solutions, 2e-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}

BOOST_FIXTURE_TEST_CASE(iterative_diagonal_solver, SolverTester) {
  typedef std::complex<float> cf;
  IterativeDiagonalSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);

  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t a = 0; a != n_ant; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != n_dir; ++d) {
        if (d == 0)
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt), uniform_sols(mt));
        else
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt) * 0.5, uniform_sols(mt) * 0.5);
      }
    }
  }

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant * 2, 1.0);
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  CheckDiagonalResults(solutions, 1e-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}

BOOST_FIXTURE_TEST_CASE(full_jones_solver, SolverTester) {
  typedef std::complex<float> cf;
  FullJonesSolver solver;
  solver.SetMaxIterations(max_iter);
  solver.SetAccuracy(1e-8);
  solver.SetStepSize(0.2);
  solver.SetNThreads(4);
  solver.SetPhaseOnly(false);
  solver.Initialize(n_ant, n_dir, n_chan, n_chan_blocks, ant1s, ant2s);
  DiagonalConstraint diagonal_constraint(4);
  solver.AddConstraint(diagonal_constraint);

  std::mt19937 mt(0);
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t a = 0; a != n_ant; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != n_dir; ++d) {
        if (d == 0)
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt), uniform_sols(mt));
        else
          input_solutions[(a * n_dir + d) * 2 + p] =
              cf(uniform_sols(mt) * 0.5, uniform_sols(mt) * 0.5);
      }
    }
  }

  FillData();

  SolverBase::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(n_chan_blocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(n_dir * n_ant * 4, 0.0);
    for (size_t i = 0; i != n_dir * n_ant * 4; i += 4) {
      vec[i] = 1.0;
      vec[i + 3] = 1.0;
    }
  }

  // Call the solver
  result = solver.Solve(data, weights, std::move(model_data), solutions, 0.0,
                        nullptr);

  // Convert full matrices to diagonals
  std::vector<std::vector<std::complex<double>>> diagonals(solutions);
  for (size_t chBlock = 0; chBlock != solutions.size(); ++chBlock) {
    for (size_t s = 0; s != solutions[chBlock].size() / 4; s++) {
      diagonals[chBlock][s * 2] = solutions[chBlock][s * 4];
      diagonals[chBlock][s * 2 + 1] = solutions[chBlock][s * 4 + 3];
    }
    diagonals[chBlock].resize(diagonals[chBlock].size() / 2);
  }

  CheckDiagonalResults(diagonals, 2e-2);
  BOOST_CHECK_EQUAL(result.iterations, max_iter + 1);
}

BOOST_AUTO_TEST_SUITE_END()
