// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <vector>
#include <complex>

#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

#include "../../DiagonalSolver.h"
#include "../../MultiDirSolver.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using aocommon::MC2x2;

class SolverTester {
 public:
  SolverTester()
      : inputSolutions(nAnt * nDir * 2),
        modelDataStore(nTimes),
        dataStore(nTimes),
        weightStore(nTimes) {
    for (size_t a1 = 0; a1 != nAnt; ++a1) {
      for (size_t a2 = a1 + 1; a2 != nAnt; ++a2) {
        ant1s.push_back(a1);
        ant2s.push_back(a2);
      }
    }
  }

  void FillData() {
    std::uniform_real_distribution<float> uniformData(-1.0, 1.0);
    std::mt19937 mt(0);
    for (size_t timestep = 0; timestep != nTimes; ++timestep) {
      aocommon::UVector<cf>& timeData = dataStore[timestep];
      timeData.assign(nPol * nChan * nBl, 0);
      aocommon::UVector<float>& timeWeights = weightStore[timestep];
      timeWeights.assign(nPol * nChan * nBl, 0);
      // modelTimeData has dimensions [direction][pol x ch x bl]
      std::vector<aocommon::UVector<cf>>& modelTimeData =
          modelDataStore[timestep];

      for (size_t d = 0; d != nDir; ++d) {
        modelTimeData.emplace_back(nPol * nChan * nBl);
        aocommon::UVector<cf>& thisDirection = modelTimeData[d];

        for (size_t bl = 0; bl != nBl; ++bl) {
          for (size_t ch = 0; ch != nChan; ++ch) {
            const size_t matrixIndex = (bl * nChan + ch) * 4;
            thisDirection[matrixIndex + 0] =
                cf(uniformData(mt), uniformData(mt));
            thisDirection[matrixIndex + 1] =
                cf(uniformData(mt) * 0.1, uniformData(mt) * 0.1);
            thisDirection[matrixIndex + 2] =
                cf(uniformData(mt) * 0.1, uniformData(mt) * 0.1);
            thisDirection[matrixIndex + 3] =
                cf(1.5 * uniformData(mt), 1.5 * uniformData(mt));
          }
        }
      }
      size_t baselineIndex = 0;
      for (size_t a1 = 0; a1 != nAnt; ++a1) {
        for (size_t a2 = a1 + 1; a2 != nAnt; ++a2) {
          for (size_t ch = 0; ch != nChan; ++ch) {
            MC2x2 perturbedModel = MC2x2::Zero();
            for (size_t d = 0; d != nDir; ++d) {
              MC2x2 val(&modelTimeData[d][(baselineIndex * nChan + ch) * 4]);
              MC2x2 left(inputSolutions[(a1 * nDir + d) * 2 + 0], 0.0, 0.0,
                         inputSolutions[(a1 * nDir + d) * 2 + 1]);
              MC2x2 right(inputSolutions[(a2 * nDir + d) * 2 + 0], 0.0, 0.0,
                          inputSolutions[(a2 * nDir + d) * 2 + 1]);
              MC2x2 res;
              MC2x2::ATimesB(res, left, val);
              MC2x2::ATimesHermB(left, res,
                                 right);  // left is used as scratch for result
              perturbedModel += left;
            }
            for (size_t p = 0; p != 4; ++p) {
              timeData[(baselineIndex * nChan + ch) * 4 + p] =
                  perturbedModel[p];
              timeWeights[(baselineIndex * nChan + ch) * 4 + p] = 1.0;
            }
          }
          ++baselineIndex;
        }
      }

      data.push_back(timeData.data());
      weights.push_back(timeWeights.data());
      std::vector<cf*> modelPtrs;
      modelPtrs.reserve(modelData.size());
      for (aocommon::UVector<cf>& modelDir : modelTimeData)
        modelPtrs.push_back(modelDir.data());
      modelData.push_back(modelPtrs);
    }
  }

  void CheckDiagonalResults(
      const std::vector<std::vector<std::complex<double>>>& solutions) {
    for (size_t ch = 0; ch != nChanBlocks; ++ch) {
      for (size_t ant = 0; ant != nAnt; ++ant) {
        for (size_t d = 0; d != nDir; ++d) {
          std::complex<double> solX0 = solutions[ch][d * 2];
          std::complex<double> solY0 = solutions[ch][d * 2 + 1];
          std::complex<double> inpX0 = inputSolutions[d * 2];
          std::complex<double> inpY0 = inputSolutions[d * 2 + 1];

          std::complex<double> solX = solutions[ch][(d + ant * nDir) * 2];
          std::complex<double> solY = solutions[ch][(d + ant * nDir) * 2 + 1];
          std::complex<double> inpX = inputSolutions[(d + ant * nDir) * 2];
          std::complex<double> inpY = inputSolutions[(d + ant * nDir) * 2 + 1];

          // Compare the squared quantities, because the phase has an ambiguity
          BOOST_CHECK_CLOSE(std::norm(solX), std::norm(inpX), 2e-2);
          BOOST_CHECK_CLOSE(std::norm(solY), std::norm(inpY), 2e-2);

          // Reference to antenna0 to check if relative phase is correct
          BOOST_CHECK_CLOSE((solX * std::conj(solX0)).real(),
                            (inpX * std::conj(inpX0)).real(), 2e-2);
          BOOST_CHECK_CLOSE((solY * std::conj(solY0)).real(),
                            (inpY * std::conj(inpY0)).real(), 2e-2);
        }
      }
    }
  }

  typedef std::complex<float> cf;
  size_t nPol = 4, nAnt = 50, nDir = 3, nChan = 10, nChanBlocks = 4,
         nTimes = 50, nBl = nAnt * (nAnt - 1) / 2, maxIter = 100;
  std::vector<int> ant1s, ant2s;
  std::vector<cf> inputSolutions;
  std::vector<cf*> data;
  std::vector<float*> weights;
  std::vector<std::vector<cf*>> modelData;
  std::vector<std::vector<aocommon::UVector<cf>>> modelDataStore;
  std::vector<aocommon::UVector<cf>> dataStore;
  std::vector<aocommon::UVector<float>> weightStore;
};

BOOST_AUTO_TEST_SUITE(solvers)

BOOST_FIXTURE_TEST_CASE(diagonal_solver, SolverTester) {
  typedef std::complex<float> cf;
  DiagonalSolver solver;
  solver.set_max_iterations(maxIter);
  solver.set_accuracy(1e-8);
  solver.set_step_size(0.2);
  solver.set_nthreads(4);
  solver.set_phase_only(false);
  solver.init(nAnt, nDir, nChan, ant1s, ant2s);
  solver.set_channel_blocks(nChanBlocks);

  std::mt19937 mt;
  std::uniform_real_distribution<float> uniformSols(1.0, 2.0);
  for (size_t a = 0; a != nAnt; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != nDir; ++d) {
        if (d == 0)
          inputSolutions[(a * nDir + d) * 2 + p] =
              cf(uniformSols(mt), uniformSols(mt));
        else
          inputSolutions[(a * nDir + d) * 2 + p] =
              cf(uniformSols(mt) * 0.5, uniformSols(mt) * 0.5);
      }
    }
  }

  FillData();

  DiagonalSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(nChanBlocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(nDir * nAnt * 2, 1.0);
  }

  // Call the solver
  result = solver.process(data, weights, modelData, solutions, 0.0, nullptr);

  CheckDiagonalResults(solutions);
  BOOST_CHECK_EQUAL(result.iterations, maxIter + 1);
}

BOOST_FIXTURE_TEST_CASE(full_jones_solver, SolverTester) {
  typedef std::complex<float> cf;
  MultiDirSolver solver;
  solver.set_max_iterations(maxIter);
  solver.set_accuracy(1e-8);
  solver.set_step_size(0.2);
  solver.set_nthreads(4);
  solver.set_phase_only(false);
  solver.init(nAnt, nDir, nChan, ant1s, ant2s);
  solver.set_channel_blocks(nChanBlocks);
  DiagonalConstraint diagonal_constraint(4);
  solver.add_constraint(&diagonal_constraint);

  std::mt19937 mt(0);
  std::uniform_real_distribution<float> uniformSols(1.0, 2.0);
  for (size_t a = 0; a != nAnt; ++a) {
    for (size_t p = 0; p != 2; ++p) {
      for (size_t d = 0; d != nDir; ++d) {
        if (d == 0)
          inputSolutions[(a * nDir + d) * 2 + p] =
              cf(uniformSols(mt), uniformSols(mt));
        else
          inputSolutions[(a * nDir + d) * 2 + p] =
              cf(uniformSols(mt) * 0.5, uniformSols(mt) * 0.5);
      }
    }
  }

  FillData();

  MultiDirSolver::SolveResult result;
  std::vector<std::vector<std::complex<double>>> solutions(nChanBlocks);

  // Initialize unit-matrices as initial values
  for (auto& vec : solutions) {
    vec.assign(nDir * nAnt * 4, 0.0);
    for (size_t i = 0; i != nDir * nAnt * 4; i += 4) {
      vec[i] = 1.0;
      vec[i + 3] = 1.0;
    }
  }

  // Call the solver
  result = solver.processFullMatrix(data, weights, modelData, solutions, 0.0,
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

  CheckDiagonalResults(diagonals);
  BOOST_CHECK_EQUAL(result.iterations, maxIter + 1);
}

BOOST_AUTO_TEST_SUITE_END()
