// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#ifndef TEC_CONSTRAINT_H
#define TEC_CONSTRAINT_H

#include <vector>

#include "Constraint.h"
#include "PieceWisePhaseFitter.h"

#ifdef AOPROJECT
#include "PhaseFitter.h"
#else
#include "../DPPP/PhaseFitter.h"
#endif

#include <vector>
#include <ostream>

class TECConstraintBase : public Constraint {
 public:
  enum Mode {
    /** Solve for both a (differential) TEC and an XX/YY-common scalar */
    TECAndCommonScalarMode,
    /** Solve only for a (differential) TEC value */
    TECOnlyMode
  };

  TECConstraintBase(Mode mode);

  /** Initialize metadata with frequencies, resize some members.
   * Should be called after InitializeDimensions.
   */
  void initialize(const double* frequencies);

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

#endif
