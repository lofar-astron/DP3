#ifndef TEC_CONSTRAINT_H
#define TEC_CONSTRAINT_H

#ifdef AOPROJECT
#include "Constraint.h"
#include "PhaseFitter.h"
#include "PieceWisePhaseFitter.h"
#else
#include <DPPP_DDECal/Constraint.h>
#include <DPPP_DDECal/PieceWisePhaseFitter.h>
#include <DPPP/PhaseFitter.h>
#endif

#include <vector>

class TECConstraintBase : public Constraint
{
public:
  enum Mode {
    /** Solve for both a (differential) TEC and an XX/YY-common scalar */
    TECAndCommonScalarMode,
    /** Solve only for a (differential) TEC value */
    TECOnlyMode
  };
  
  TECConstraintBase(Mode mode);

  void initialize(size_t nAntennas, size_t nDirections, 
                  size_t nChannelBlocks, const double* frequencies);
  
protected:
  virtual void initializeChild() { }
  
  void applyReferenceAntenna(std::vector<std::vector<dcomplex> >& solutions) const;
  
  Mode _mode;
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<PhaseFitter> _phaseFitters;
};

class TECConstraint : public TECConstraintBase
{
public:
  TECConstraint(Mode mode) : TECConstraintBase(mode) { }

  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                       double time);
};

class ApproximateTECConstraint : public TECConstraint
{
public:
  ApproximateTECConstraint(Mode mode) :
    TECConstraint(mode),
    _finishedApproximateStage(false),
    _fittingChunkSize(0),
    _maxApproxIters(50)
    { }

  virtual void PrepareIteration(bool hasReachedPrecision, size_t iteration, bool finalIter) {
    _finishedApproximateStage = hasReachedPrecision || finalIter || iteration >= _maxApproxIters;
  }
  
  virtual bool Satisfied() const { return _finishedApproximateStage; }
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                       double time);
  
  void SetFittingChunkSize(size_t fittingChunkSize)
  { _fittingChunkSize = fittingChunkSize; }
  
  void SetMaxApproximatingIterations(size_t maxApproxIters)
  { _maxApproxIters = maxApproxIters; }
protected:
  virtual void initializeChild();
  
private:
  bool _finishedApproximateStage;
  std::vector<PieceWisePhaseFitter> _pwFitters;
  std::vector<std::vector<double>> _threadData;
  std::vector<std::vector<double>> _threadFittedData;
  std::vector<std::vector<double>> _threadWeights;
  size_t _fittingChunkSize, _maxApproxIters;
};

#endif

