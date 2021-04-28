// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDE_SOLVER_BASE_H
#define DDE_SOLVER_BASE_H

#include "../constraints/Constraint.h"

#include "../linear_solvers/LLSSolver.h"

#include <boost/algorithm/string.hpp>

#include <complex>
#include <vector>
#include <memory>

namespace dp3 {
namespace base {

class DPBuffer;
class SolverBuffer;

class SolverBase {
 public:
  typedef std::complex<double> DComplex;
  typedef std::complex<float> Complex;

  class Matrix final {
   public:
    Matrix() : data_(), columns_(0) {}
    Matrix(size_t columns, size_t rows)
        : data_(columns * rows, 0.0), columns_(columns) {}
    void SetZero() { data_.assign(data_.size(), Complex(0.0, 0.0)); }
    Complex& operator()(size_t column, size_t row) {
      return data_[column + row * columns_];
    }
    Complex* data() { return data_.data(); }

   private:
    std::vector<Complex> data_;
    size_t columns_;
  };

  struct SolveResult {
    size_t iterations = 0;
    size_t constraint_iterations = 0;
    std::vector<std::vector<Constraint::Result>> results;
  };

  SolverBase();

  virtual ~SolverBase() {}

  /**
   * Prepares the solver with the given dimensionality info
   * and antenna mapping.
   * The antenna arrays map the data provided in @Solve to the antennas.
   */
  virtual void Initialize(size_t nAntennas, size_t nDirections,
                          size_t nChannels, size_t nChannelBlocks,
                          const std::vector<int>& ant1,
                          const std::vector<int>& ant2);

  /**
   * Solves multi-directional Jones matrices. Takes the (single) measured data
   * and the (multi-directional) model data, and solves the optimization
   * problem that minimizes the norm of the differences.
   *
   * @param solver_buffer Buffer with unweighted data, weights and model data.
   * @param solutions The per-channel and per-antenna solutions.
   * solutions[ch] is a pointer for channelblock ch to antenna x directions x
   * pol solutions.
   * @param statStream Optional pointer to a stream for displaying statistics.
   */
  virtual SolveResult Solve(const SolverBuffer& solver_buffer,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) = 0;

  /**
   * @return The number of polarizations in the solution.
   */
  virtual size_t NSolutionPolarizations() const = 0;

  /**
   * Add a constraint to the solver.
   * @param constraint A valid constraint pointer.
   * @throw std::runtime_error If the constraint is invalid.
   */
  void AddConstraint(std::unique_ptr<Constraint> constraint) {
    if (!constraint) {
      throw std::runtime_error("Solver constraint is empty.");
    }
    constraints_.push_back(std::move(constraint));
  }

  /**
   * Get the constraints for the solver.
   * @return A non-modifyable list of modifyable Constraints. All pointers in
   *         the list are valid.
   */
  const std::vector<std::unique_ptr<Constraint>>& GetConstraints() {
    return constraints_;
  }

  /**
   * If enabled, the solver will perform steps along the complex
   * circle, instead of moving freely through complex space.
   * See the implementation of @ref MakeStep().
   */
  void SetPhaseOnly(bool phase_only) { phase_only_ = phase_only; }

  /**
   * Max nr of iterations (stopping criterion).
   * @{
   */
  size_t GetMaxIterations() const { return max_iterations_; }
  void SetMaxIterations(size_t max_iterations) {
    max_iterations_ = max_iterations;
  }
  /** @} */

  /**
   * Min nr of iterations before stopping.
   * @{
   */
  size_t GetMinIterations() const { return min_iterations_; }
  void SetMinIterations(size_t min_iterations) {
    min_iterations_ = min_iterations;
  }
  /** @} */

  /**
   * Required relative accuracy.
   * @{
   */
  void SetAccuracy(double accuracy) { accuracy_ = accuracy; }
  double GetAccuracy() const { return accuracy_; }
  /** @} */

  /**
   * Required relative accuracy for the constraints to finish.
   */
  void SetConstraintAccuracy(double constraint_accuracy) {
    constraint_accuracy_ = constraint_accuracy;
  }

  /**
   * The step size taken each iteration. Higher values might
   * make convergence faster, but may cause instability.
   * @{
   */
  void SetStepSize(double step_size) { step_size_ = step_size; }
  double GetStepSize() const { return step_size_; }
  /** @} */

  /**
   * Whether stalling of the solutions should abort the solving.
   * @{
   */
  void SetDetectStalling(bool detect_stalling) {
    detect_stalling_ = detect_stalling;
  }
  bool GetDetectStalling() const { return detect_stalling_; }
  /** @} */

  /**
   * Number of threads to use in parts that can be parallelized.
   * The solving is parallelized over channel blocks.
   */
  virtual void SetNThreads(size_t n_threads) {
    n_threads_ = n_threads;
    for (std::unique_ptr<Constraint>& constraint : constraints_) {
      constraint->SetNThreads(n_threads);
    }
  }
  size_t GetNThreads() const { return n_threads_; }

  /**
   * Output timing information to a stream.
   */
  void GetTimings(std::ostream& os, double duration) const;

  void SetLLSSolverType(LLSSolverType solver_type, double min_tolerance,
                        double max_tolerance);

 protected:
  void Step(const std::vector<std::vector<DComplex>>& solutions,
            std::vector<std::vector<DComplex>>& next_solutions) const;

  bool DetectStall(size_t iteration,
                   const std::vector<double>& step_magnitudes) const;

  static void MakeSolutionsFinite1Pol(
      std::vector<std::vector<DComplex>>& solutions);

  static void MakeSolutionsFinite2Pol(
      std::vector<std::vector<DComplex>>& solutions);

  static void MakeSolutionsFinite4Pol(
      std::vector<std::vector<DComplex>>& solutions);

  void PrepareConstraints();

  bool ApplyConstraints(size_t iteration, double time,
                        bool has_previously_converged, SolveResult& result,
                        std::vector<std::vector<DComplex>>& next_solutions,
                        std::ostream* stat_stream) const;

  /**
   * Assign the solutions in nextSolutions to the solutions.
   * @returns whether the solutions have converged. Appends the current step
   * magnitude to step_magnitudes
   */
  bool AssignSolutions(std::vector<std::vector<DComplex>>& solutions,
                       const std::vector<std::vector<DComplex>>& new_solutions,
                       bool use_constraint_accuracy, double& avg_abs_diff,
                       std::vector<double>& step_magnitudes, size_t nPol) const;

  bool ReachedStoppingCriterion(
      size_t iteration, bool has_converged, bool constraints_satisfied,
      const std::vector<double>& step_magnitudes) const {
    bool has_stalled = false;
    if (detect_stalling_ && constraints_satisfied)
      has_stalled = DetectStall(iteration, step_magnitudes);

    const bool is_ready = iteration >= max_iterations_ ||
                          (has_converged && constraints_satisfied) ||
                          has_stalled;
    return iteration >= min_iterations_ && is_ready;
  }

  size_t NAntennas() const { return n_antennas_; }
  size_t NDirections() const { return n_directions_; }
  size_t NChannels() const { return n_channels_; }
  size_t NChannelBlocks() const { return n_channel_blocks_; }
  size_t NBaselines() const { return ant1_.size(); }
  int AntennaIndex1(size_t baseline) const { return ant1_[baseline]; }
  int AntennaIndex2(size_t baseline) const { return ant2_[baseline]; }

  /**
   * @param block A channel block index, less than NChannelBlocks().
   * @return The index of the first channel in the given channel block.
   */
  size_t FirstChannel(size_t block) const {
    return block * n_channels_ / n_channel_blocks_;
  }

  /**
   * Create an LLSSolver with the given matrix dimensions.
   * Set the tolerance using 'iteration_fraction' and 'solver_precision'.
   */
  std::unique_ptr<LLSSolver> CreateLLSSolver(size_t m, size_t n, size_t nrhs,
                                             double iteration_fraction,
                                             double solver_precision) const;

  size_t n_antennas_;
  size_t n_directions_;
  size_t n_channels_;
  size_t n_channel_blocks_;
  std::vector<int> ant1_, ant2_;

  /**
   * Calibration setup
   * @{
   */
  size_t min_iterations_;
  size_t max_iterations_;
  size_t n_threads_;
  double accuracy_;
  double constraint_accuracy_;
  double step_size_;
  bool detect_stalling_;

 private:
  bool phase_only_;
  std::vector<std::unique_ptr<Constraint>> constraints_;

  LLSSolverType lls_solver_type_;
  double lls_min_tolerance_;
  double lls_max_tolerance_;
  /** @} */
};

}  // namespace base
}  // namespace dp3

#endif
