#ifndef ROTATIONANDDIAGONAL_CONSTRAINT_H
#define ROTATIONANDDIAGONAL_CONSTRAINT_H

#ifdef AOPROJECT
#include "Constraint.h"
#else
#include <DPPP_DDECal/Constraint.h>
#endif

#include <vector>

namespace LOFAR {

class RotationAndDiagonalConstraint : public Constraint
{
public:
  RotationAndDiagonalConstraint() {};
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);

  virtual void InitializeDimensions(size_t nAntennas,
                                    size_t nDirections,
                                    size_t nChannelBlocks);

  virtual void SetWeights(std::vector<double>&);

private:
  std::vector<Constraint::Result> _resTemplate;
};

} // namespace LOFAR

#endif

