// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DiagonalSolver.h"
#include "QRSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using namespace aocommon;

#include <iomanip>
#include <iostream>

DiagonalSolver::DiagonalSolver()
    : _nAntennas(0),
      _nDirections(0),
      _nChannels(0),
      _nChannelBlocks(0),
      _maxIterations(100),
      _accuracy(1e-5),
      _constraintAccuracy(1e-4),
      _stepSize(0.2),
      _detectStalling(true),
      _phaseOnly(false) {}

void DiagonalSolver::init(size_t nAntennas, size_t nDirections,
                          size_t nChannels, const std::vector<int>& ant1,
                          const std::vector<int>& ant2) {
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannels = nChannels;
  _nChannelBlocks = nChannels;
  _ant1 = ant1;
  _ant2 = ant2;
  _buffer.SetDimensions(nDirections, nChannels, ant1.size());
}

void DiagonalSolver::makeStep(
    const std::vector<std::vector<DComplex>>& solutions,
    std::vector<std::vector<DComplex>>& nextSolutions) const {
  // Move the solutions towards nextSolutions
  // (the moved solutions are stored in 'nextSolutions')
  ParallelFor<size_t> loop(_nThreads);
  loop.Run(0, _nChannelBlocks, [&](size_t chBlock, size_t /*thread*/) {
    for (size_t i = 0; i != nextSolutions[chBlock].size(); ++i) {
      if (_phaseOnly) {
        // In phase only mode, a step is made along the complex circle,
        // towards the shortest direction.
        double phaseFrom = std::arg(solutions[chBlock][i]);
        double distance = std::arg(nextSolutions[chBlock][i]) - phaseFrom;
        if (distance > M_PI)
          distance = distance - 2.0 * M_PI;
        else if (distance < -M_PI)
          distance = distance + 2.0 * M_PI;
        nextSolutions[chBlock][i] =
            std::polar(1.0, phaseFrom + _stepSize * distance);
      } else {
        nextSolutions[chBlock][i] = solutions[chBlock][i] * (1.0 - _stepSize) +
                                    nextSolutions[chBlock][i] * _stepSize;
      }
    }
  });
}

void DiagonalSolver::makeSolutionsFinite2pol(
    std::vector<std::vector<DComplex>>& solutions) {
  for (std::vector<DComplex>& solVector : solutions) {
    // Find the average abs solution for this channel
    size_t count = 0;
    double average[2] = {0.0, 0.0};
    for (std::vector<DComplex>::iterator iter = solVector.begin();
         iter != solVector.end(); iter += 2) {
      if (isfinite(*iter) && isfinite(*(iter + 1))) {
        for (size_t p = 0; p != 2; ++p) average[p] += std::abs(iter[0]);
        ++count;
      }
    }
    if (count == 0) {
      average[0] = 1.0;
      average[1] = 1.0;
    } else {
      for (size_t p = 0; p != 2; ++p) average[p] /= count;
    }
    for (std::vector<DComplex>::iterator iter = solVector.begin();
         iter != solVector.end(); iter += 2) {
      if (!isfinite(*iter) || !isfinite(*(iter + 1))) {
        for (size_t p = 0; p != 2; ++p) iter[p] = average[p];
      }
    }
  }
}

template <size_t NPol>
bool DiagonalSolver::assignSolutions(
    std::vector<std::vector<DComplex>>& solutions,
    std::vector<std::vector<DComplex>>& nextSolutions,
    bool useConstraintAccuracy, double& avgAbsDiff,
    std::vector<double>& stepMagnitudes) const {
  avgAbsDiff = 0.0;
  //  Calculate the norm of the difference between the old and new solutions
  size_t n = 0;
  for (size_t chBlock = 0; chBlock < _nChannelBlocks; ++chBlock) {
    for (size_t i = 0; i != solutions[chBlock].size(); i += NPol) {
      // A normalized squared difference is calculated between the solutions of
      // this and the previous step:
      //   sqrt { 1/n sum over | (t1 - t0) t0^(-1) |_2 }
      // This criterion is scale independent: all solutions can be scaled
      // without affecting the number of iterations. Also, when the polarized
      // version is given scalar matrices, it will use the same number of
      // iterations as the scalar version.
      if (NPol == 2) {
        DComplex s[2] = {solutions[chBlock][i], solutions[chBlock][i + 1]};
        DComplex sInv[2] = {1.0 / s[0], 1.0 / s[1]};
        if (isfinite(sInv[0]) && isfinite(sInv[1])) {
          DComplex ns[2] = {nextSolutions[chBlock][i],
                            nextSolutions[chBlock][i + 1]};
          ns[0] = (ns[0] - s[0]) * sInv[0];
          ns[1] = (ns[1] - s[1]) * sInv[1];
          double sumabs = 0.0;
          for (size_t p = 0; p != NPol; ++p) {
            sumabs += std::abs(ns[p]);
          }
          if (std::isfinite(sumabs)) {
            avgAbsDiff += sumabs;
            n += 2;
          }
        }
        for (size_t p = 0; p != NPol; ++p) {
          solutions[chBlock][i + p] = nextSolutions[chBlock][i + p];
        }
      }
    }
  }

  // The stepsize is taken out, so that a small stepsize won't cause
  // a premature stopping criterion.
  double stepMagnitude = (n == 0 ? 0 : avgAbsDiff / _stepSize / n);
  stepMagnitudes.emplace_back(stepMagnitude);

  if (useConstraintAccuracy)
    return stepMagnitude <= _constraintAccuracy;
  else {
    return stepMagnitude <= _accuracy;
  }
}

DiagonalSolver::SolveResult DiagonalSolver::process(
    const std::vector<Complex*>& dataNoW, const std::vector<float*>& weights,
    const std::vector<std::vector<Complex*>>& modelDataNoW,
    std::vector<std::vector<DComplex>>& solutions, double time,
    std::ostream* statStream) {
  const size_t nTimes = dataNoW.size();

  _buffer.CopyAndWeight(dataNoW, weights, modelDataNoW);

  for (size_t i = 0; i != _constraints.size(); ++i)
    _constraints[i]->PrepareIteration(false, 0, false);

  std::vector<std::vector<DComplex>> nextSolutions(_nChannelBlocks);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != _nChannelBlocks) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  result._results.resize(_constraints.size());

  // Model matrix ant x [N x D] and visibility vector ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<Matrix>> gTimesCs(_nChannelBlocks);
  std::vector<std::vector<std::vector<Complex>>> vs(_nChannelBlocks);
  for (size_t chBlock = 0; chBlock != _nChannelBlocks; ++chBlock) {
    nextSolutions[chBlock].resize(_nDirections * _nAntennas * 2);
    const size_t channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
                 channelIndexEnd = (chBlock + 1) * _nChannels / _nChannelBlocks,
                 curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas * 2);
    vs[chBlock].resize(_nAntennas * 2);

    for (size_t ant = 0; ant != _nAntennas * 2; ++ant) {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set
      // to 0.
      // X and Y polarizations are treated as two different antennas.
      size_t m = (_nAntennas * 2) * nTimes * curChannelBlockSize,
             n = _nDirections;
      gTimesCs[chBlock][ant] = Matrix(m, n);
      vs[chBlock][ant].resize(std::max(m, n));
    }
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0, constrainedIterations = 0;
  bool hasConverged = false, hasPreviouslyConverged = false,
       constraintsSatisfied = false, hasStalled = false;

  std::vector<double> stepMagnitudes;
  stepMagnitudes.reserve(_maxIterations);

  do {
    makeSolutionsFinite2pol(solutions);

    ParallelFor<size_t> loop(_nThreads);
    loop.Run(0, _nChannelBlocks, [&](size_t chBlock, size_t /*thread*/) {
      performIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                       solutions[chBlock], nextSolutions[chBlock]);
    });

    makeStep(solutions, nextSolutions);

    if (statStream) {
      (*statStream) << iteration << '\t';
    }

    constraintsSatisfied = true;

    for (size_t i = 0; i != _constraints.size(); ++i) {
      // PrepareIteration() might change Satisfied(), and since we always want
      // to iterate at least once more when a constraint is not yet satisfied,
      // we evaluate Satisfied() before preparing.
      constraintsSatisfied =
          _constraints[i]->Satisfied() && constraintsSatisfied;
      _constraints[i]->PrepareIteration(hasPreviouslyConverged, iteration,
                                        iteration + 1 >= _maxIterations);
      result._results[i] =
          _constraints[i]->Apply(nextSolutions, time, statStream);
    }

    if (!constraintsSatisfied) constrainedIterations = iteration + 1;

    double avgSquaredDiff;
    hasConverged =
        assignSolutions<2>(solutions, nextSolutions, !constraintsSatisfied,
                           avgSquaredDiff, stepMagnitudes);
    if (statStream) {
      (*statStream) << stepMagnitudes.back() << '\t' << avgSquaredDiff << '\n';
    }
    iteration++;

    hasPreviouslyConverged = hasConverged || hasPreviouslyConverged;

    if (_detectStalling && constraintsSatisfied)
      hasStalled = detectStall(iteration, stepMagnitudes);

  } while (iteration < _maxIterations &&
           (!hasConverged || !constraintsSatisfied) && !hasStalled);

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged iterations can be distinguished from converged ones.
  if ((!hasConverged || !constraintsSatisfied) && !hasStalled)
    result.iterations = iteration + 1;
  else
    result.iterations = iteration;
  result.constraintIterations = constrainedIterations;
  return result;
}

void DiagonalSolver::performIteration(size_t channelBlockIndex,
                                      std::vector<Matrix>& gTimesCs,
                                      std::vector<std::vector<Complex>>& vs,
                                      const std::vector<DComplex>& solutions,
                                      std::vector<DComplex>& nextSolutions) {
  for (size_t ant = 0; ant != _nAntennas * 2; ++ant) {
    gTimesCs[ant].zeros();
    std::fill(vs[ant].begin(), vs[ant].end(), 0.0);
  }

  const size_t channelIndexStart =
      channelBlockIndex * _nChannels / _nChannelBlocks;
  const size_t channelIndexEnd =
      (channelBlockIndex + 1) * _nChannels / _nChannelBlocks;
  const size_t curChannelBlockSize = channelIndexEnd - channelIndexStart,
               nTimes = _buffer.Data().size();

  // The following loop fills the matrices for all antennas
  std::vector<const Complex*> modelPtrs(_nDirections);
  for (size_t timeIndex = 0; timeIndex != nTimes; ++timeIndex) {
    for (size_t baseline = 0; baseline != _ant1.size(); ++baseline) {
      size_t antenna1 = _ant1[baseline];
      size_t antenna2 = _ant2[baseline];
      if (antenna1 != antenna2) {
        for (size_t d = 0; d != _nDirections; ++d) {
          modelPtrs[d] =
              &_buffer.ModelData()[timeIndex][d]
                                  [(channelIndexStart + baseline * _nChannels) *
                                   4];
        }
        const Complex* dataPtr =
            &_buffer.Data()[timeIndex]
                           [(channelIndexStart + baseline * _nChannels) * 4];
        for (size_t ch = channelIndexStart; ch != channelIndexEnd; ++ch) {
          for (size_t p = 0; p != 4; ++p) {
            size_t p1 = p / 2;
            size_t p2 = p % 2;
            const size_t dataIndex1 =
                             ch - channelIndexStart +
                             (timeIndex + (antenna1 * 2 + p1) * nTimes) *
                                 curChannelBlockSize,
                         dataIndex2 =
                             ch - channelIndexStart +
                             (timeIndex + (antenna2 * 2 + p2) * nTimes) *
                                 curChannelBlockSize;
            Matrix& gTimesC1 = gTimesCs[antenna1 * 2 + p1];
            std::vector<Complex>& v1 = vs[antenna1 * 2 + p1];
            Matrix& gTimesC2 = gTimesCs[antenna2 * 2 + p2];
            std::vector<Complex>& v2 = vs[antenna2 * 2 + p2];

            for (size_t d = 0; d != _nDirections; ++d) {
              std::complex<double> predicted = *modelPtrs[d];

              size_t solIndex1 = (antenna1 * _nDirections + d) * 2 + p1;
              size_t solIndex2 = (antenna2 * _nDirections + d) * 2 + p2;
              gTimesC2(dataIndex1, d) = std::conj(
                  solutions[solIndex1] * predicted);  // using a* b* = (ab)*
              gTimesC1(dataIndex2, d) =
                  std::conj(solutions[solIndex2]) * predicted;

              ++modelPtrs[d];  // Goto the next polarization of this 2x2 matrix.
            }
            v1[dataIndex2] = *dataPtr;
            v2[dataIndex1] = std::conj(*dataPtr);
            ++dataPtr;  // Goto the next polarization of this 2x2 matrix.
          }
        }
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.
  size_t m = _nAntennas * 2 * nTimes * curChannelBlockSize;
  size_t n = _nDirections, nrhs = 1;
  QRSolver solver(m, n, nrhs);
  for (size_t ant = 0; ant != _nAntennas; ++ant) {
    for (size_t pol = 0; pol != 2; ++pol) {
      // solve x^H in [g C] x^H  = v
      bool success = solver.Solve(gTimesCs[ant * 2 + pol].data(),
                                  vs[ant * 2 + pol].data());
      std::vector<Complex>& x = vs[ant * 2 + pol];
      if (success && x[0] != Complex(0.0, 0.0)) {
        for (size_t d = 0; d != _nDirections; ++d)
          nextSolutions[(ant * _nDirections + d) * 2 + pol] = x[d];
      } else {
        for (size_t d = 0; d != _nDirections; ++d)
          nextSolutions[(ant * _nDirections + d) * 2 + pol] =
              std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

bool DiagonalSolver::detectStall(
    size_t iteration, const std::vector<double>& step_magnitudes) const {
  if (iteration < 30) {
    return false;
  } else {
    return std::abs(step_magnitudes[iteration - 1] /
                        step_magnitudes[iteration - 2] -
                    1) < 1.e-4 &&
           std::abs(step_magnitudes[iteration - 2] /
                        step_magnitudes[iteration - 3] -
                    1) < 1.e-4;
  }
}
