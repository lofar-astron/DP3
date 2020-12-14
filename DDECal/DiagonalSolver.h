// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDE_DIAGONAL_SOLVER_H
#define DDE_DIAGONAL_SOLVER_H

#include "Constraint.h"
#include "MultiDirBuffer.h"
#include "Stopwatch.h"

#include <complex>
#include <vector>
#include <memory>

class DiagonalSolver {
 public:
  typedef std::complex<double> DComplex;
  typedef std::complex<float> Complex;

  class Matrix : public std::vector<Complex> {
   public:
    Matrix() : _m(0) {}
    Matrix(size_t m, size_t n) : std::vector<Complex>(n * m, 0.0), _m(m) {}
    void zeros() { assign(size(), Complex(0.0, 0.0)); }
    Complex& operator()(size_t i, size_t j) { return (*this)[i + j * _m]; }

   private:
    size_t _m;
  };

  struct SolveResult {
    size_t iterations, constraintIterations;
    std::vector<std::vector<Constraint::Result>> _results;
  };

  DiagonalSolver();

  void init(size_t nAntennas, size_t nDirections, size_t nChannels,
            const std::vector<int>& ant1, const std::vector<int>& ant2);

  /// data[i] is a pointer to the data for timestep i. Those data are then
  /// ordered as they are in the MS (bl, chan, pol). mdata[i] is a pointer for
  /// timestep i to arrays of ndir model data pointers Each of these pointers is
  /// in the same order as the data. solutions[ch] is a pointer for channelblock
  /// ch to antenna x directions x pol(=2) solutions.
  SolveResult process(const std::vector<Complex*>& data,
                      const std::vector<float*>& weights,
                      const std::vector<std::vector<Complex*>>& modelData,
                      std::vector<std::vector<DComplex>>& solutions,
                      double time, std::ostream* statStream);

  void set_phase_only(bool phaseOnly) { _phaseOnly = phaseOnly; }

  void set_channel_blocks(size_t nChannelBlocks) {
    _nChannelBlocks = nChannelBlocks;
  }

  size_t max_iterations() const { return _maxIterations; }
  void set_max_iterations(size_t maxIterations) {
    _maxIterations = maxIterations;
  }

  void set_accuracy(double accuracy) { _accuracy = accuracy; }
  double get_accuracy() const { return _accuracy; }
  void set_constraint_accuracy(double constraintAccuracy) {
    _constraintAccuracy = constraintAccuracy;
  }
  void set_step_size(double stepSize) { _stepSize = stepSize; }
  double get_step_size() const { return _stepSize; }

  void set_detect_stalling(bool detectStalling) {
    _detectStalling = detectStalling;
  }

  bool get_detect_stalling() const { return _detectStalling; }

  void add_constraint(Constraint* constraint) {
    _constraints.push_back(constraint);
  }

  void set_nthreads(size_t nthreads) { _nThreads = nthreads; }

 private:
  void performIteration(size_t channelBlockIndex, std::vector<Matrix>& gTimesCs,
                        std::vector<std::vector<Complex>>& vs,
                        const std::vector<DComplex>& solutions,
                        std::vector<DComplex>& nextSolutions);

  void makeStep(const std::vector<std::vector<DComplex>>& solutions,
                std::vector<std::vector<DComplex>>& nextSolutions) const;

  bool detectStall(size_t iteration,
                   const std::vector<double>& stepMagnitudes) const;

  static void makeSolutionsFinite2pol(
      std::vector<std::vector<DComplex>>& solutions);

  template <typename T>
  static bool isfinite(const std::complex<T>& val) {
    return std::isfinite(val.real()) && std::isfinite(val.imag());
  }

  /**
   * Assign the solutions in nextSolutions to the solutions.
   * @returns whether the solutions have converged. Appends the current step
   * magnitude to step_magnitudes
   */
  template <size_t NPol>
  bool assignSolutions(std::vector<std::vector<DComplex>>& solutions,
                       std::vector<std::vector<DComplex>>& nextSolutions,
                       bool useConstraintAccuracy, double& sum,
                       std::vector<double>& step_magnitudes) const;

  size_t _nAntennas, _nDirections, _nChannels, _nChannelBlocks;
  std::vector<int> _ant1, _ant2;
  MultiDirBuffer _buffer;

  /// @name Calibration setup
  /// @{
  size_t _maxIterations, _nThreads;
  double _accuracy, _constraintAccuracy;
  double _stepSize;
  bool _detectStalling;
  bool _phaseOnly;
  std::vector<Constraint*> _constraints;
  /// @}
};

#endif
