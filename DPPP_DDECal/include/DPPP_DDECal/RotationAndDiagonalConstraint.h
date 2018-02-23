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
  RotationAndDiagonalConstraint();
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);

  void initialize(size_t nAntennas, size_t nDirections, size_t nChannelBlocks);

private:
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<Constraint::Result> _resTemplate;
};

} // namespace LOFAR

#endif

