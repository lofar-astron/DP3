// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBase.h"

#include <algorithm>
#include <iostream>
#include <numeric>

#include <aocommon/matrix2x2.h>
#include <aocommon/staticfor.h>
#include <aocommon/xt/span.h>

#include <xsimd/xsimd.hpp>

#include <xtensor/xview.hpp>

namespace {
template <typename T>
bool IsFinite(const std::complex<T>& val) {
  return std::isfinite(val.real()) && std::isfinite(val.imag());
}
}  // namespace

namespace dp3 {
namespace ddecal {

SolverBase::SolverBase()
    : n_antennas_(0),
      n_directions_(0),
      n_channel_blocks_(0),
      min_iterations_(0),
      max_iterations_(100),
      accuracy_(1e-5),
      constraint_accuracy_(1e-4),
      step_size_(0.2),
      detect_stalling_(true),
      step_diff_sigma_(0.1),
      phase_only_(false),
      lls_solver_type_(LLSSolverType::QR) {}

void SolverBase::Initialize(
    size_t n_antennas, const std::vector<size_t>& n_solutions_per_direction,
    size_t n_channel_blocks) {
  n_directions_ = n_solutions_per_direction.size();
  n_solutions_ = std::accumulate(n_solutions_per_direction.begin(),
                                 n_solutions_per_direction.end(), 0u);
  if (!SupportsDdSolutionIntervals() && (n_solutions_ != n_directions_)) {
    throw std::runtime_error(
        "DD interval solutions not supported for the selected solver "
        "algorithm. Please "
        "select a different solver algorithm or remove "
        "ddecal.solutions_per_direction from your parset. The most common "
        "solution for this problem is to set `solveralgorithm` to "
        "`directioniterative` for the ddecal step in your parset.");
  }

  n_antennas_ = n_antennas;
  n_channel_blocks_ = n_channel_blocks;
  assert(n_solutions_ != 0);
  assert(n_directions_ != 0);
}

bool SolverBase::DetectStall(size_t iteration,
                             const std::vector<double>& step_magnitudes) {
  if (iteration == 1) {
    /* starting a new solve cycle, so reset */
    n_var_count_ = 1;
    step_mean_ = step_magnitudes[iteration - 1];
    step_var_ = 0.0;
  } else {
    n_var_count_ += 1;
    double delta_1 = step_magnitudes[iteration - 1] - step_mean_;
    step_mean_ += delta_1 / (double)n_var_count_;
    double delta_2 = step_magnitudes[iteration - 1] - step_mean_;
    step_var_ += delta_1 * delta_2;
  }

  if (iteration < 30) {
    return false;
  } else {
    /* stalled when runtime mean of step sizes
     * is less than runtime standard deviation of step sizes */
    return step_mean_ <
           step_diff_sigma_ * std::sqrt(step_var_ / (double)n_var_count_);
  }
}

void SolverBase::Step(const std::vector<std::vector<DComplex>>& solutions,
                      SolutionTensor& next_solutions) const {
  // Move the solutions towards next_solutions
  // (the moved solutions are stored in 'next_solutions')

  const size_t n_antennas = next_solutions.shape(1);
  const size_t n_solutions = next_solutions.shape(2);
  const size_t n_polarizations = next_solutions.shape(3);
  const size_t solution_size = n_antennas * n_solutions * n_polarizations;

  // Parallelizing this loop using StaticFor proved not effective:
  // When setting kPhaseOnly to true in SolverTester.h, the total loop takes
  // about 13 microseconds in the solvers/scalar test and 17 microseconds
  // in the solvers/diagonal test. With kPhaseOnly == false, the loops
  // typically take less than 5 microseconds.
  // Using a StaticFor increases the time to more than 21 microseconds per loop,
  // when using 4 threads.
  for (size_t ch_block = 0; ch_block < n_channel_blocks_; ++ch_block) {
    assert(solutions[ch_block].size() == solution_size);
    const std::complex<double>* solution = solutions[ch_block].data();
    std::complex<double>* next_solution = &next_solutions(ch_block, 0, 0, 0);

    if (!phase_only_) {
      // Using XSimd is not necessary: The compiler vectorizes this simple loop.
      for (size_t i = 0; i < solution_size; ++i) {
        next_solution[i] =
            solution[i] * (1.0 - step_size_) + next_solution[i] * step_size_;
      }
    } else {
      // In phase only mode, a step is made along the complex circle,
      // towards the shortest direction.

      using BatchComplex = xsimd::batch<std::complex<double>>;
      constexpr size_t kBatchSize = BatchComplex::size;

      const xsimd::batch<double> zero = xsimd::broadcast(0.0);
      const xsimd::batch<double> one = xsimd::broadcast(1.0);
      const xsimd::batch<double> pi = xsimd::broadcast(M_PI);
      const xsimd::batch<double> two_pi = xsimd::broadcast(2.0 * M_PI);

      size_t i = 0;
      for (; i < (solution_size / kBatchSize) * kBatchSize; i += kBatchSize) {
        const BatchComplex solution_batch =
            BatchComplex::load_unaligned(&solution[i]);
        // next_solution should be aligned, but load_aligned does not work
        BatchComplex next_solution_batch =
            BatchComplex::load_unaligned(&next_solution[i]);

        const xsimd::batch<double> phase_from = xsimd::arg(solution_batch);
        xsimd::batch<double> distance =
            xsimd::arg(next_solution_batch) - phase_from;
        distance -= xsimd::select(xsimd::gt(distance, pi), two_pi, zero);
        distance += xsimd::select(xsimd::lt(distance, -pi), two_pi, zero);

        next_solution_batch =
            xsimd::polar(one, phase_from + step_size_ * distance);
        next_solution_batch.store_unaligned(&next_solution[i]);
      }

      for (; i < solution_size; ++i) {
        const double phase_from = std::arg(solution[i]);
        double distance = std::arg(next_solution[i]) - phase_from;
        if (distance > M_PI)
          distance -= 2.0 * M_PI;
        else if (distance < -M_PI)
          distance += 2.0 * M_PI;
        next_solution[i] = std::polar(1.0, phase_from + step_size_ * distance);
      }
    }
  }
}

void SolverBase::PrepareConstraints() {
  for (std::unique_ptr<Constraint>& c : constraints_) {
    c->PrepareIteration(false, 0, false);
  }
}

bool SolverBase::ApplyConstraints(size_t iteration, double time,
                                  bool has_previously_converged,
                                  SolveResult& result,
                                  SolutionTensor& next_solutions,
                                  std::ostream* stat_stream) const {
  SolutionSpan next_solutions_span = aocommon::xt::CreateSpan(next_solutions);
  return ApplyConstraints(iteration, time, has_previously_converged, result,
                          next_solutions_span, stat_stream);
}

bool SolverBase::ApplyConstraints(size_t iteration, double time,
                                  bool has_previously_converged,
                                  SolveResult& result,
                                  SolutionSpan& next_solutions,
                                  std::ostream* stat_stream) const {
  bool constraints_satisfied = true;

  result.results.resize(constraints_.size());
  auto result_iterator = result.results.begin();

  for (const std::unique_ptr<Constraint>& c : constraints_) {
    // PrepareIteration() might change Satisfied(), and since we always want
    // to iterate at least once more when a constraint is not yet satisfied,
    // we evaluate Satisfied() before preparing.
    constraints_satisfied = c->Satisfied() && constraints_satisfied;
    c->PrepareIteration(has_previously_converged, iteration,
                        iteration + 1 >= GetMaxIterations());
    *result_iterator = c->Apply(next_solutions, time, stat_stream);
    ++result_iterator;
  }

  // If still not satisfied, at least iteration+1 constrained iterations
  // were required.
  if (!constraints_satisfied) result.constraint_iterations = iteration + 1;

  return constraints_satisfied;
}

bool SolverBase::AssignSolutions(std::vector<std::vector<DComplex>>& solutions,
                                 SolutionSpan& newSolutions,
                                 bool useConstraintAccuracy, double& avgAbsDiff,
                                 std::vector<double>& stepMagnitudes) const {
  assert(newSolutions.shape(0) == NChannelBlocks());
  assert(newSolutions.shape(1) == NAntennas());
  assert(newSolutions.shape(2) == NSolutions());
  assert(newSolutions.shape(3) == NSolutionPolarizations());
  avgAbsDiff = 0.0;
  //  Calculate the norm of the difference between the old and new solutions
  const size_t n_solution_pols = NSolutionPolarizations();
  size_t n = 0;
  const auto new_solutions_view = xt::reshape_view(
      newSolutions,
      {newSolutions.shape(0), newSolutions.shape(1) * newSolutions.shape(2),
       newSolutions.shape(3)});
  for (size_t chBlock = 0; chBlock < n_channel_blocks_; ++chBlock) {
    for (size_t i = 0; i != solutions[chBlock].size(); i += n_solution_pols) {
      // A normalized squared difference is calculated between the solutions of
      // this and the previous step:
      //   sqrt { 1/n sum over | (t1 - t0) t0^(-1) |_2 }
      // This criterion is scale independent: all solutions can be scaled
      // without affecting the number of iterations. Also, when the polarized
      // version is given scalar matrices, it will use the same number of
      // iterations as the scalar version.
      if (n_solution_pols == 1) {
        if (solutions[chBlock][i] != 0.0) {
          const double a = std::abs(
              (new_solutions_view(chBlock, i, 0) - solutions[chBlock][i]) /
              solutions[chBlock][i]);
          if (std::isfinite(a)) {
            avgAbsDiff += a;
            ++n;
          }
        }
      } else if (n_solution_pols == 2) {
        const DComplex s[2] = {solutions[chBlock][i],
                               solutions[chBlock][i + 1]};
        const DComplex sInv[2] = {1.0 / s[0], 1.0 / s[1]};
        if (IsFinite(sInv[0]) && IsFinite(sInv[1])) {
          DComplex ns[2] = {new_solutions_view(chBlock, i / 2, 0),
                            new_solutions_view(chBlock, i / 2, 1)};
          ns[0] = (ns[0] - s[0]) * sInv[0];
          ns[1] = (ns[1] - s[1]) * sInv[1];
          const double sumabs = std::abs(ns[0]) + std::abs(ns[1]);
          if (std::isfinite(sumabs)) {
            avgAbsDiff += sumabs;
            n += 2;
          }
        }
      } else {
        assert(n_solution_pols == 4);
        const aocommon::MC2x2 s(&solutions[chBlock][i]);
        aocommon::MC2x2 sInv(s);
        if (sInv.Invert()) {
          aocommon::MC2x2 ns(&new_solutions_view(chBlock, i / 4, 0));
          ns -= s;
          ns *= sInv;
          const double sumabs = std::abs(ns[0]) + std::abs(ns[1]) +
                                std::abs(ns[2]) + std::abs(ns[3]);
          if (std::isfinite(sumabs)) {
            avgAbsDiff += sumabs;
            n += 4;
          }
        }
      }
    }

    xt::adapt(solutions[chBlock]) =
        xt::flatten(xt::view(newSolutions, chBlock));
  }

  // The stepsize is taken out, so that a small stepsize won't cause
  // a premature stopping criterion.
  double stepMagnitude = (n == 0 ? 0 : avgAbsDiff / step_size_ / n);
  stepMagnitudes.emplace_back(stepMagnitude);

  if (useConstraintAccuracy)
    return stepMagnitude <= constraint_accuracy_;
  else {
    return stepMagnitude <= accuracy_;
  }
}
bool SolverBase::AssignSolutions(std::vector<std::vector<DComplex>>& solutions,
                                 SolutionTensor& newSolutions,
                                 bool useConstraintAccuracy, double& avgAbsDiff,
                                 std::vector<double>& stepMagnitudes) const {
  SolutionSpan new_solutions_span = aocommon::xt::CreateSpan(newSolutions);
  return AssignSolutions(solutions, new_solutions_span, useConstraintAccuracy,
                         avgAbsDiff, stepMagnitudes);
}
void SolverBase::MakeSolutionsFinite1Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  // Parallelizing this loop using StaticFor proved not effective:
  // In the solvers/scalar test, the total loop took 7 microseconds while
  // it took 16 microseconds, on average, when using StaticFor.
  for (std::vector<DComplex>& solution_vector : solutions) {
    // Find the average solutions for this channel
    size_t count = 0;
    double average = 0.0;
    for (const std::complex<double>& solution : solution_vector) {
      if (IsFinite(solution)) {
        average += std::abs(solution);
        ++count;
      }
    }
    // If no solution was found, replace it by the average abs value
    if (count == 0)
      average = 1.0;
    else
      average /= count;
    for (std::complex<double>& solution : solution_vector) {
      if (!IsFinite(solution)) solution = average;
    }
  }
}

void SolverBase::MakeSolutionsFinite2Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  // Parallelizing this loop using StaticFor proved not effective:
  // In the solvers/diagonal test, the total loop took 8 microseconds while
  // it took 14 microseconds, on average, when using StaticFor.
  for (std::vector<DComplex>& solution_vector : solutions) {
    // Find the average abs solution for this channel
    size_t count = 0;
    double average[2] = {0.0, 0.0};
    for (std::vector<DComplex>::iterator iter = solution_vector.begin();
         iter != solution_vector.end(); iter += 2) {
      if (IsFinite(*iter) && IsFinite(*(iter + 1))) {
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
    for (std::vector<DComplex>::iterator iter = solution_vector.begin();
         iter != solution_vector.end(); iter += 2) {
      if (!IsFinite(*iter) || !IsFinite(*(iter + 1))) {
        for (size_t p = 0; p != 2; ++p) iter[p] = average[p];
      }
    }
  }
}

void SolverBase::MakeSolutionsFinite4Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  // Parallelizing this loop using StaticFor proved not effective:
  // In the solvers/full_jones test, the total loop took 10 microseconds
  // while it took 15 microseconds, on average, when using StaticFor.
  for (std::vector<DComplex>& solution_vector : solutions) {
    // Find the average abs solution for this channel
    size_t count = 0;
    double average[4] = {0.0, 0.0, 0.0, 0.0};
    for (std::vector<DComplex>::iterator iter = solution_vector.begin();
         iter != solution_vector.end(); iter += 4) {
      if (aocommon::Matrix2x2::IsFinite(&*iter)) {
        for (size_t p = 0; p != 4; ++p) average[p] += std::abs(iter[0]);
        ++count;
      }
    }
    if (count == 0) {
      average[0] = 1.0;
      average[1] = 0.0;
      average[2] = 0.0;
      average[3] = 1.0;
    } else {
      for (size_t p = 0; p != 4; ++p) average[p] /= count;
    }
    for (std::vector<DComplex>::iterator iter = solution_vector.begin();
         iter != solution_vector.end(); iter += 4) {
      if (!aocommon::Matrix2x2::IsFinite(&*iter)) {
        for (size_t p = 0; p != 4; ++p) iter[p] = average[p];
      }
    }
  }
}

void SolverBase::GetTimings(std::ostream& os, double duration) const {
  for (const std::unique_ptr<Constraint>& constraint : constraints_) {
    constraint->GetTimings(os, duration);
  }
}

void SolverBase::SetLLSSolverType(const LLSSolverType solver_type) {
  lls_solver_type_ = solver_type;
}

std::unique_ptr<LLSSolver> SolverBase::CreateLLSSolver(
    const size_t m, const size_t n, const size_t nrhs) const {
  return LLSSolver::Make(lls_solver_type_, m, n, nrhs);
}

}  // namespace ddecal
}  // namespace dp3
