// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBase.h"

#include <aocommon/parallelfor.h>
#include <aocommon/matrix2x2.h>

#include <algorithm>
#include <iostream>
#include <numeric>

using aocommon::ParallelFor;

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
        "DD interval solutions not supported for selected Solver. Please "
        "select a different solver algorithm or remove "
        "ddecal.solutions_per_direction from your parset.");
  }

  n_antennas_ = n_antennas;
  n_channel_blocks_ = n_channel_blocks;
  assert(n_solutions_ != 0);
  assert(n_directions_ != 0);
}

bool SolverBase::DetectStall(size_t iteration,
                             const std::vector<double>& step_magnitudes) const {
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

void SolverBase::Step(const std::vector<std::vector<DComplex>>& solutions,
                      std::vector<std::vector<DComplex>>& nextSolutions) const {
  // Move the solutions towards nextSolutions
  // (the moved solutions are stored in 'nextSolutions')
  ParallelFor<size_t> loop(n_threads_);
  loop.Run(0, n_channel_blocks_, [&](size_t chBlock, size_t /*thread*/) {
    for (size_t i = 0; i != nextSolutions[chBlock].size(); ++i) {
      if (phase_only_) {
        // In phase only mode, a step is made along the complex circle,
        // towards the shortest direction.
        double phaseFrom = std::arg(solutions[chBlock][i]);
        double distance = std::arg(nextSolutions[chBlock][i]) - phaseFrom;
        if (distance > M_PI)
          distance = distance - 2.0 * M_PI;
        else if (distance < -M_PI)
          distance = distance + 2.0 * M_PI;
        nextSolutions[chBlock][i] =
            std::polar(1.0, phaseFrom + step_size_ * distance);
      } else {
        nextSolutions[chBlock][i] = solutions[chBlock][i] * (1.0 - step_size_) +
                                    nextSolutions[chBlock][i] * step_size_;
      }
    }
  });
}

void SolverBase::PrepareConstraints() {
  for (std::unique_ptr<Constraint>& c : constraints_) {
    c->PrepareIteration(false, 0, false);
  }
}

bool SolverBase::ApplyConstraints(
    size_t iteration, double time, bool has_previously_converged,
    SolveResult& result, std::vector<std::vector<DComplex>>& next_solutions,
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

bool SolverBase::AssignSolutions(
    std::vector<std::vector<DComplex>>& solutions,
    const std::vector<std::vector<DComplex>>& newSolutions,
    bool useConstraintAccuracy, double& avgAbsDiff,
    std::vector<double>& stepMagnitudes) const {
  avgAbsDiff = 0.0;
  //  Calculate the norm of the difference between the old and new solutions
  const size_t n_solution_pols = NSolutionPolarizations();
  size_t n = 0;
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
          double a =
              std::abs((newSolutions[chBlock][i] - solutions[chBlock][i]) /
                       solutions[chBlock][i]);
          if (std::isfinite(a)) {
            avgAbsDiff += a;
            ++n;
          }
        }
        solutions[chBlock][i] = newSolutions[chBlock][i];
      } else if (n_solution_pols == 2) {
        DComplex s[2] = {solutions[chBlock][i], solutions[chBlock][i + 1]};
        DComplex sInv[2] = {1.0 / s[0], 1.0 / s[1]};
        if (IsFinite(sInv[0]) && IsFinite(sInv[1])) {
          DComplex ns[2] = {newSolutions[chBlock][i],
                            newSolutions[chBlock][i + 1]};
          ns[0] = (ns[0] - s[0]) * sInv[0];
          ns[1] = (ns[1] - s[1]) * sInv[1];
          double sumabs = 0.0;
          for (size_t p = 0; p != n_solution_pols; ++p) {
            sumabs += std::abs(ns[p]);
          }
          if (std::isfinite(sumabs)) {
            avgAbsDiff += sumabs;
            n += 2;
          }
        }
        for (size_t p = 0; p != n_solution_pols; ++p) {
          solutions[chBlock][i + p] = newSolutions[chBlock][i + p];
        }
      } else {  // NPol == 4
        aocommon::MC2x2 s(&solutions[chBlock][i]), sInv(s);
        if (sInv.Invert()) {
          aocommon::MC2x2 ns(&newSolutions[chBlock][i]);
          ns -= s;
          ns *= sInv;
          double sumabs = 0.0;
          for (size_t p = 0; p != n_solution_pols; ++p) {
            sumabs += std::abs(ns[p]);
          }
          if (std::isfinite(sumabs)) {
            avgAbsDiff += sumabs;
            n += 4;
          }
        }
        for (size_t p = 0; p != n_solution_pols; ++p) {
          solutions[chBlock][i + p] = newSolutions[chBlock][i + p];
        }
      }
    }
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

void SolverBase::MakeSolutionsFinite1Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  for (std::vector<DComplex>& solVector : solutions) {
    // Find the average solutions for this channel
    size_t count = 0;
    double average = 0.0;
    for (const std::complex<double>& solution : solVector) {
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
    for (std::complex<double>& solution : solVector) {
      if (!IsFinite(solution)) solution = average;
    }
  }
}

void SolverBase::MakeSolutionsFinite2Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  for (std::vector<DComplex>& solVector : solutions) {
    // Find the average abs solution for this channel
    size_t count = 0;
    double average[2] = {0.0, 0.0};
    for (std::vector<DComplex>::iterator iter = solVector.begin();
         iter != solVector.end(); iter += 2) {
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
    for (std::vector<DComplex>::iterator iter = solVector.begin();
         iter != solVector.end(); iter += 2) {
      if (!IsFinite(*iter) || !IsFinite(*(iter + 1))) {
        for (size_t p = 0; p != 2; ++p) iter[p] = average[p];
      }
    }
  }
}

void SolverBase::MakeSolutionsFinite4Pol(
    std::vector<std::vector<DComplex>>& solutions) {
  for (std::vector<DComplex>& sol_vector : solutions) {
    // Find the average abs solution for this channel
    size_t count = 0;
    double average[4] = {0.0, 0.0, 0.0, 0.0};
    for (std::vector<DComplex>::iterator iter = sol_vector.begin();
         iter != sol_vector.end(); iter += 4) {
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
    for (std::vector<DComplex>::iterator iter = sol_vector.begin();
         iter != sol_vector.end(); iter += 4) {
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
