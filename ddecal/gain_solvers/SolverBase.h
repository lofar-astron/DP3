// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_BASE_H
#define DDECAL_SOLVER_BASE_H

#include <algorithm>
#include <cassert>
#include <complex>
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <vector>

#include <aocommon/recursivefor.h>

#include "../constraints/Constraint.h"
#include "../linear_solvers/LLSSolver.h"
#include "SolveData.h"

namespace dp3 {
namespace ddecal {

class SolverBase {
 public:
  typedef std::complex<double> DComplex;
  typedef std::complex<float> Complex;

  class Matrix final {
   public:
    Matrix() : data_(), columns_(0) {}
    Matrix(size_t columns, size_t rows)
        : data_(columns * rows, Complex(0.0, 0.0)), columns_(columns) {}
    void SetZero() { std::fill(data_.begin(), data_.end(), Complex(0.0, 0.0)); }
    Complex& operator()(size_t column, size_t row) {
      return data_[column + row * columns_];
    }
    Complex* data() { return data_.data(); }

    // Re-initialize the object to a new size.
    //
    // Like the constructor this sets all data element to 0+0i.
    void Reset(size_t columns, size_t rows) {
      size_t size = columns * rows;
      // Minimize the number of elements to modify.
      if (size < data_.size()) {
        data_.resize(size);
        SetZero();
      } else {
        SetZero();
        data_.resize(size, Complex(0.0, 0.0));
      }
      columns_ = columns;
    }

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

  virtual ~SolverBase() = default;

  /**
   * Prepares the solver with the given dimensionality info
   * and antenna mapping.
   * The antenna arrays map the data provided in @Solve to the antennas.
   */
  virtual void Initialize(size_t n_antennas,
                          const std::vector<size_t>& n_solutions_per_direction,
                          size_t n_channel_blocks);

  /**
   * @return The number of polarizations in the solution.
   */
  virtual size_t NSolutionPolarizations() const = 0;

  /**
   * Add a constraint to the solver.
   * @param constraint A valid constraint pointer, must not be nullptr.
   */
  void AddConstraint(std::unique_ptr<Constraint> constraint) {
    assert(constraint);
    constraints_.push_back(std::move(constraint));
  }

  /**
   * Get the constraints for the solver.
   * @return A non-modifiable list of modifiable Constraints. All pointers in
   *         the list are valid.
   */
  const std::vector<std::unique_ptr<Constraint>>& GetConstraints() {
    return constraints_;
  }

  /**
   * If enabled, the solver will perform steps along the complex
   * circle, instead of moving freely through complex space.
   * See the implementation of @ref MakeStep().
   *{
   */
  bool GetPhaseOnly() const { return phase_only_; }
  void SetPhaseOnly(bool phase_only) { phase_only_ = phase_only; }
  /** @} */

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
   * @{
   */
  void SetConstraintAccuracy(double constraint_accuracy) {
    constraint_accuracy_ = constraint_accuracy;
  }
  double GetConstraintAccuracy() const { return constraint_accuracy_; }
  /** @} */

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
  void SetDetectStalling(bool detect_stalling, double step_diff_sigma) {
    detect_stalling_ = detect_stalling;
    step_diff_sigma_ = step_diff_sigma;
  }
  bool GetDetectStalling() const { return detect_stalling_; }
  /** @} */

  /**
   * Output timing information to a stream.
   */
  void GetTimings(std::ostream& os, double duration) const;

  void SetLLSSolverType(LLSSolverType solver_type);
  LLSSolverType GetLLSSolverType() const { return lls_solver_type_; }

  /**
   * True if this solver supports direction-dependent solution intervals. In
   * that case, the solver can solve different directions with different
   * solution intervals. More specific, each timestep of each direction can
   * individually be set to belong to a specific solution.
   */
  virtual bool SupportsDdSolutionIntervals() const { return false; }

  /**
   * Returns a list of solvers that this solver uses and for which
   * constraints should be set up. The @ref HybridSolver overrides
   * this function to make it possible to initialize the solvers it
   * combines.
   */
  virtual std::vector<SolverBase*> ConstraintSolvers() { return {this}; }

  /**
   * Solves multi-directional Jones matrices. Takes the (single) measured data
   * and the (multi-directional) model data, and solves the optimization
   * problem that minimizes the norm of the differences.
   *
   * @param data Buffer with weighted data and model data.
   * @param solutions The per-channel and per-antenna solutions.
   * solutions[ch] is a pointer for channelblock ch to antenna x directions x
   * pol solutions.
   * @param statStream Optional pointer to a stream for displaying statistics.
   */
  virtual SolveResult Solve(const FullSolveData& data,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) {
    throw std::logic_error(
        "Full-visibility solver called for a solver that does not "
        "support full-visibility solving");
  }

  virtual SolveResult Solve(const DuoSolveData& data,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) {
    throw std::logic_error(
        "Duo-visibility (xx/yy) solver called for a solver that does not "
        "support duo-visibility solving");
  }

  virtual SolveResult Solve(const UniSolveData& data,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) {
    throw std::logic_error(
        "Single-visibility (xx/yy) solver called for a solver that does not "
        "support single-visibility solving");
  }

  /** Calls @ref SetSolutionWeights() for all constraints. */
  void SetDdConstraintWeights(const std::vector<std::vector<double>>& weights);

 protected:
  void Step(const std::vector<std::vector<DComplex>>& solutions,
            SolutionTensor& next_solutions) const;

  bool DetectStall(size_t iteration,
                   const std::vector<double>& step_magnitudes);

  static void MakeSolutionsFinite1Pol(
      std::vector<std::vector<DComplex>>& solutions);

  static void MakeSolutionsFinite2Pol(
      std::vector<std::vector<DComplex>>& solutions);

  static void MakeSolutionsFinite4Pol(
      std::vector<std::vector<DComplex>>& solutions);

  void PrepareConstraints();

  bool ApplyConstraints(size_t iteration, double time,
                        bool has_previously_converged, SolveResult& result,
                        SolutionTensor& next_solutions,
                        std::ostream* stat_stream) const;
  bool ApplyConstraints(size_t iteration, double time,
                        bool has_previously_converged, SolveResult& result,
                        SolutionSpan& next_solutions,
                        std::ostream* stat_stream) const;
  /**
   * Assign the solutions in nextSolutions to the solutions.
   * @returns whether the solutions have converged. Appends the current step
   * magnitude to step_magnitudes
   */
  bool AssignSolutions(std::vector<std::vector<DComplex>>& solutions,
                       SolutionTensor& new_solutions,
                       bool use_constraint_accuracy, double& avg_abs_diff,
                       std::vector<double>& step_magnitudes) const;
  bool AssignSolutions(std::vector<std::vector<DComplex>>& solutions,
                       SolutionSpan& new_solutions,
                       bool use_constraint_accuracy, double& avg_abs_diff,
                       std::vector<double>& step_magnitudes) const;

  bool ReachedStoppingCriterion(size_t iteration, bool has_converged,
                                bool constraints_satisfied,
                                const std::vector<double>& step_magnitudes) {
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
  size_t NChannelBlocks() const { return n_channel_blocks_; }
  /**
   * Total number of solutions over all directions
   * This might be different from n_directions_ when using direction-dependent
   * intervals.
   */
  size_t NSubSolutions() const { return n_sub_solutions_; }
  /**
   * Total number of visibilities over all channel blocks
   */
  size_t NVisibilities() const {
    return NChannelBlocks() * NAntennas() * NSubSolutions() *
           NSolutionPolarizations();
  }

  /**
   * Returns the number of available thread for one channel. This assumes that
   * the channel blocks is the main direction that is parallelized over (as
   * all solvers are currently doing). For example, if there are 20 threads
   * available, and the number of channel blocks is 10, then 2 is returned,
   * because each channel block can use two threads to fill up the cpu.
   * The returned value is at least 1.
   */
  size_t NSubThreads() const;

  /**
   * If useful, make a RecursiveFor object. The existance of a RecursiveFor
   * makes the StaticFor use the RecursiveFor. Using a StaticFor inside a
   * RecursiveFor allows nested parallelization, but it is slower. Therefore, if
   * no nested parallelization can be done, don't create a RecursiveFor. See
   * also @ref NSubThreads().
   */
  std::unique_ptr<aocommon::RecursiveFor> MakeOptionalRecursiveFor() const;

  /**
   * Create an LLSSolver with the given matrix dimensions.
   * Set the tolerance using 'iteration_fraction' and 'solver_precision'.
   */
  std::unique_ptr<LLSSolver> CreateLLSSolver(size_t m, size_t n,
                                             size_t nrhs) const;

 private:
  size_t n_antennas_;
  size_t n_directions_;
  size_t n_channel_blocks_;
  size_t n_sub_solutions_;

  /**
   * Calibration setup
   * @{
   */
  size_t min_iterations_;
  size_t max_iterations_;
  double accuracy_;
  double constraint_accuracy_;
  double step_size_;
  bool detect_stalling_;
  double step_diff_sigma_;

  bool phase_only_;
  std::vector<std::unique_ptr<Constraint>> constraints_;

  LLSSolverType lls_solver_type_;

  /**
   * variables for calculating mean/variance of step size vector
   * to determine stalling
   */
  size_t n_var_count_;
  double step_mean_;
  double step_var_;
  /** @} */
};

}  // namespace ddecal
}  // namespace dp3

#endif
