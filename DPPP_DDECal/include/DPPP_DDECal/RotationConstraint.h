#ifndef ROTATION_CONSTRAINT_H
#define ROTATION_CONSTRAINT_H

#ifdef AOPROJECT
#include "Constraint.h"
#else
#include <DPPP_DDECal/Constraint.h>
#endif

#include <vector>

namespace LOFAR {

class RotationConstraint : public Constraint
{
public:
  RotationConstraint();
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time);

  void initialize(size_t nAntennas, size_t nDirections, size_t nChannelBlocks);

  /* Compute the rotation from a 2x2 full jones solution */
  static double get_rotation(std::complex<double>* data);

private:
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<Constraint::Result> _resTemplate;
};

} // namespace LOFAR

#endif

