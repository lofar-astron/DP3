// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <complex>
#include <memory>
#include <set>
#include <vector>
#include <ostream>

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
    std::vector<double> vals;
    std::vector<double> weights;
    std::string
        axes;  ///< Comma-separated string with axis names, fastest varying last
    std::vector<size_t> dims;
    std::string name;
  };

  Constraint()
      : n_antennas_(0), n_directions_(0), n_channel_blocks_(0), n_threads_(0) {}

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
   * @param solutions is an array of array, such that:
   * - solutions[ch] is a pointer for channelblock ch to antenna x directions x
   * pol solutions.
   * - pol is the dimension with the fastest changing index.
   * @param time Central time of interval.
   */
  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream) = 0;

  /**
   * Perform common constraint initialization. Should be overridden when
   * something more than assigning dimensions is needed (e.g. resizing vectors).
   * @param frequencies For each channel block, the mean frequency.
   */
  virtual void Initialize(size_t nAntennas, size_t nDirections,
                          const std::vector<double>& frequencies) {
    n_antennas_ = nAntennas;
    n_directions_ = nDirections;
    n_channel_blocks_ = frequencies.size();
  }

  /**
   * Set weights. The vector should contain an array of size nAntennas *
   * nChannelBlocks, where the channel index varies fastest.
   */
  virtual void SetWeights([[maybe_unused]] const std::vector<double>& weights) {
  }

  void SetNThreads(size_t n_threads) { n_threads_ = n_threads; }

  virtual void GetTimings([[maybe_unused]] std::ostream& os,
                          [[maybe_unused]] double duration) const {}

 protected:
  size_t NAntennas() const { return n_antennas_; }
  size_t NDirections() const { return n_directions_; }
  size_t NChannelBlocks() const { return n_channel_blocks_; }
  size_t NThreads() const { return n_threads_; }

  static bool isfinite(const dcomplex& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  }

 private:
  size_t n_antennas_, n_directions_, n_channel_blocks_, n_threads_;
};

/**
 * @brief This class constraints the amplitudes of the solution to be unity, but
 * keeps the phase.
 */
class PhaseOnlyConstraint : public Constraint {
 public:
  PhaseOnlyConstraint(){};

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream);
};

/**
 * @brief This class constraints the phases of the solution to be zero, but
 * keeps the amplitude information.
 */
class AmplitudeOnlyConstraint : public Constraint {
 public:
  AmplitudeOnlyConstraint(){};

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream) final override;
};

class DiagonalConstraint : public Constraint {
 public:
  DiagonalConstraint(size_t polsPerSolution)
      : _polsPerSolution(polsPerSolution){};

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream) final override;

 private:
  const size_t _polsPerSolution;
};

/**
 * @brief This constraint averages the solutions of several groups of antennas,
 * so that antennas within the same group have equal solutions.
 *
 * The DDE solver can use this constraint e.g. to average the solutions of
 * the core antennas. Core antennas are determined by a given maximum distance
 * from a reference antenna. The reference antenna is by default the first
 * antenna. This constraint is meant to force all core stations to
 * have the same solution, thereby decreasing the noise in their solutions.
 */
class AntennaConstraint : public Constraint {
 public:
  AntennaConstraint() {}

  void SetAntennaSets(std::vector<std::set<size_t>>&& antenna_sets) {
    antenna_sets_ = std::move(antenna_sets);
  }

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream) final override;

 private:
  std::vector<std::set<size_t>> antenna_sets_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
