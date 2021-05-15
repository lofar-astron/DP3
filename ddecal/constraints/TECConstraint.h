// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef TEC_CONSTRAINT_H
#define TEC_CONSTRAINT_H

#include <vector>

#include "Constraint.h"
#include "PieceWisePhaseFitter.h"

#include "../../base/PhaseFitter.h"

#include <vector>
#include <ostream>

namespace dp3 {
namespace ddecal {

class TECConstraintBase : public Constraint {
 public:
  enum Mode {
    /** Solve for both a (differential) TEC and an XX/YY-common scalar */
    TECAndCommonScalarMode,
    /** Solve only for a (differential) TEC value */
    TECOnlyMode
  };

  TECConstraintBase(Mode mode);

  void Initialize(size_t nAntennas, size_t nDirections,
                  const std::vector<double>& frequencies) override;

  /** Propagate weights to the phase fitters */
  virtual void SetWeights(const std::vector<double>& weights) final override;

  /** Setter for doPhaseReference */
  void setDoPhaseReference(bool doPhaseReference) {
    _doPhaseReference = doPhaseReference;
  }

 protected:
  virtual void initializeChild() {}

  void applyReferenceAntenna(
      std::vector<std::vector<dcomplex> >& solutions) const;

  Mode _mode;
  bool _doPhaseReference;
  std::vector<PhaseFitter> _phaseFitters;
  std::vector<double> _weights;
};

class TECConstraint : public TECConstraintBase {
 public:
  TECConstraint(Mode mode) : TECConstraintBase(mode) {}

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex> >& solutions, double time,
      std::ostream* statStream) override;
};

class ApproximateTECConstraint : public TECConstraint {
 public:
  ApproximateTECConstraint(Mode mode)
      : TECConstraint(mode),
        _finishedApproximateStage(false),
        _fittingChunkSize(0),
        _maxApproxIters(50) {}

  virtual void PrepareIteration(bool hasReachedPrecision, size_t iteration,
                                bool finalIter) final override {
    _finishedApproximateStage =
        hasReachedPrecision || finalIter || iteration >= _maxApproxIters;
    for (size_t thread = 0; thread != _phaseFitters.size(); ++thread) {
      std::fill(
          _phaseFitters[thread].WeightData(),
          _phaseFitters[thread].WeightData() + _phaseFitters[thread].Size(),
          1.0);
    }
  }

  virtual bool Satisfied() const final override {
    return _finishedApproximateStage;
  }

  virtual std::vector<Result> Apply(
      std::vector<std::vector<dcomplex> >& solutions, double time,
      std::ostream* statStream) final override;

  void SetFittingChunkSize(size_t fittingChunkSize) {
    _fittingChunkSize = fittingChunkSize;
  }

  void SetMaxApproximatingIterations(size_t maxApproxIters) {
    _maxApproxIters = maxApproxIters;
  }

 protected:
  virtual void initializeChild() final override;

 private:
  bool _finishedApproximateStage;
  std::vector<PieceWisePhaseFitter> _pwFitters;
  std::vector<std::vector<double> > _threadData;
  std::vector<std::vector<double> > _threadFittedData;
  std::vector<std::vector<double> > _threadWeights;
  size_t _fittingChunkSize, _maxApproxIters;
};

}  // namespace ddecal
}  // namespace dp3

#endif
