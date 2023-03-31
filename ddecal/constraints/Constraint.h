// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_CONSTRAINT_H_
#define DP3_DDECAL_CONSTRAINT_H_

#include <complex>
#include <numeric>
#include <ostream>
#include <vector>

#include <xtensor/xtensor.hpp>

#include "../Solutions.h"

namespace dp3 {
namespace ddecal {

/**
 * \brief This class is the base class for classes that implement a constraint
 * on calibration solutions. Constraints are used to increase the converge of
 * calibration by applying these inside the solving step.
 *
 * The MultiDirSolver class uses this class for constrained calibration.
 */
class Constraint {
 public:
  typedef std::complex<double> dcomplex;
  struct Result {
   public:
    /// Both vals and weights have the dimensions described in dims and axes.
    std::vector<double> vals;
    std::vector<double> weights;
    /// Comma-separated string with axis names, fastest varying last.
    std::string axes;
    std::vector<size_t> dims;
    std::string name;
  };

  virtual ~Constraint() {}

  /**
   * Function that initializes the constraint for the next calibration
   * iteration. It should be called each time all antenna solutions have been
   * calculated, but before the constraint has been applied to all those antenna
   * solutions.
   *
   * Unlike Apply(), this method is not thread safe.
   *
   * @param bool This can be used to specify whether the previous solution
   * "step" is smaller than the requested precision, i.e. calibration with the
   * constrained has converged. This allows a constraint to apply its constraint
   * in steps: apply a better-converging constraint as long as the solutions are
   * far from the correct answer, then switch to a different constraint when
   * hasReachedPrecision=true.
   */
  virtual void PrepareIteration([[maybe_unused]] bool hasReachedPrecision,
                                [[maybe_unused]] size_t iteration,
                                [[maybe_unused]] bool finalIter) {}

  /**
   * Whether the constraint has been satisfied. The calibration process will
   * continue at least as long as Satisfied()=false, and performs at least one
   * more iteration after Satisfied()=true. Together with SetPrecisionReached(),
   * this can make the algorithm change the constraining method based on amount
   * of convergence.
   */
  virtual bool Satisfied() const { return true; }

  /**
   * This method applies the constraints to the solutions.
   * @param solutions A 4D array with dimensions n_channel_blocks x n_antennas
   * x n_solutions x n_pol solutions.
   * n_pol is the dimension with the fastest changing index.
   * @param time Central time of interval.
   */
  virtual std::vector<Result> Apply(SolutionSpan& solutions, double time,
                                    std::ostream* statStream) = 0;

  /**
   * Perform common constraint initialization. Should be overridden when
   * something more than assigning dimensions is needed (e.g. resizing vectors).
   * @param frequencies For each channel block, the mean frequency.
   */
  virtual void Initialize(size_t n_antennas,
                          const std::vector<uint32_t>& solutions_per_direction,
                          const std::vector<double>& frequencies) {
    n_antennas_ = n_antennas;
    solutions_per_direction_ = solutions_per_direction;
    n_channel_blocks_ = frequencies.size();
    n_solutions_ = std::accumulate(solutions_per_direction.begin(),
                                   solutions_per_direction.end(), 0u);
  }

  /**
   * Set weights. The vector should contain an array of size nAntennas *
   * nChannelBlocks, where the channel index varies fastest.
   */
  virtual void SetWeights([[maybe_unused]] const std::vector<double>& weights) {
  }

  /**
   * Set the number of threads for parallel loops etc.
   * @param n_threads Desired number of threads. If it is zero, it becomes one.
   */
  void SetNThreads(size_t n_threads) {
    n_threads_ = std::min(n_threads, size_t(1));
  }

  virtual void GetTimings([[maybe_unused]] std::ostream& os,
                          [[maybe_unused]] double duration) const {}

  size_t NAntennas() const { return n_antennas_; }
  size_t NDirections() const { return solutions_per_direction_.size(); }
  /**
   * Number of solutions over all directions, i.e. the sum over
   * solutions_per_direction_.
   */
  size_t NSolutions() const { return n_solutions_; }
  size_t NChannelBlocks() const { return n_channel_blocks_; }
  size_t NThreads() const { return n_threads_; }

  static bool isfinite(const dcomplex& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  }

 private:
  size_t n_antennas_ = 0;
  size_t n_channel_blocks_ = 0;
  size_t n_threads_ = 1;
  size_t n_solutions_ = 0;
  std::vector<uint32_t> solutions_per_direction_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
