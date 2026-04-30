// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_TEC_CONSTRAINT_H_
#define DP3_DDECAL_TEC_CONSTRAINT_H_

#include <vector>

#include "Constraint.h"
#include "PieceWisePhaseFitter.h"

#include "base/PhaseFitter.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

class TecConstraint : public Constraint {
 public:
  enum class Mode {
    /** Solve for both a (differential) TEC and an XX/YY-common scalar */
    kTecAndCommonScalar,
    /** Solve only for a (differential) TEC value */
    kTecOnly
  };

  TecConstraint(Mode mode);

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  /** Propagate weights to the phase fitters */
  virtual void SetWeights(const std::vector<double>& weights) final override;

  /** Setter for doPhaseReference */
  void setDoPhaseReference(bool doPhaseReference) {
    do_phase_reference_ = doPhaseReference;
  }

  void Apply(SolutionSpan& solutions, double time) override;

  std::vector<ConstraintResult> GetResult() const override { return results_; };

 protected:
  virtual void initializeChild() {}

  void applyReferenceAntenna(SolutionSpan& solutions) const;

  Mode mode_;
  bool do_phase_reference_;
  std::vector<PhaseFitter> phase_fitters_;
  std::vector<double> weights_;
  std::vector<ConstraintResult> results_;
};

class ApproximateTECConstraint final : public TecConstraint {
 public:
  ApproximateTECConstraint(Mode mode)
      : TecConstraint(mode),
        finished_approximate_stage_(false),
        fitting_chunk_size_(0),
        max_approx_iters_(50) {}

  void PrepareIteration(bool has_reached_precision, size_t iteration,
                        bool final_iter) override {
    finished_approximate_stage_ =
        has_reached_precision || final_iter || iteration >= max_approx_iters_;
    for (size_t thread = 0; thread != phase_fitters_.size(); ++thread) {
      std::fill(
          phase_fitters_[thread].WeightData(),
          phase_fitters_[thread].WeightData() + phase_fitters_[thread].Size(),
          1.0);
    }
  }

  bool Satisfied() const override { return finished_approximate_stage_; }

  void Apply(SolutionSpan& solutions, double time) override;

  void SetFittingChunkSize(size_t fitting_chunk_size) {
    fitting_chunk_size_ = fitting_chunk_size;
  }

  void SetMaxApproximatingIterations(size_t max_approx_iters) {
    max_approx_iters_ = max_approx_iters;
  }

  std::vector<ConstraintResult> GetResult() const override {
    return finished_approximate_stage_ ? results_
                                       : std::vector<ConstraintResult>();
  };

 protected:
  void initializeChild() override;

 private:
  bool finished_approximate_stage_;
  std::vector<PieceWisePhaseFitter> pw_fitters_;
  std::vector<std::vector<double>> thread_data_;
  std::vector<std::vector<double>> thread_fitted_data_;
  std::vector<std::vector<double>> thread_weights_;
  size_t fitting_chunk_size_, max_approx_iters_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
