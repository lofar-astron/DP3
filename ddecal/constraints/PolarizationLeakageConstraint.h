#ifndef DP3_DDECAL_POLARIZATION_LEAKAGE_CONSTRAINT_H_
#define DP3_DDECAL_POLARIZATION_LEAKAGE_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>
#include <ostream>

namespace dp3::ddecal {

/**
 * Constraint that turns a solution matrix into a polarization leakage
 * matrix, i.e. of the form [ 1 a ; b 1 ].
 */
class PolarizationLeakageConstraint final : public Constraint {
 public:
  PolarizationLeakageConstraint() = default;

  std::vector<ConstraintResult> Apply(SolutionSpan& solutions, double time,
                                      std::ostream* stat_stream) override;
};

}  // namespace dp3::ddecal

#endif
