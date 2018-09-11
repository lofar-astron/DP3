#ifndef ROTATION_CONSTRAINT_H
#define ROTATION_CONSTRAINT_H

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace LOFAR {

class RotationConstraint : public Constraint
{
public:
  RotationConstraint() {};
  
  virtual std::vector<Result> Apply(
                    std::vector<std::vector<dcomplex> >& solutions,
                    double time, std::ostream* statStream);

  virtual void InitializeDimensions(size_t nAntennas,
                                    size_t nDirections,
                                    size_t nChannelBlocks);

  virtual void SetWeights(const std::vector<double>& weights);

  /* Compute the rotation from a 2x2 full jones solution */
  static double get_rotation(std::complex<double>* data);

private:
  std::vector<Constraint::Result> _res;
};

} // namespace LOFAR

#endif

