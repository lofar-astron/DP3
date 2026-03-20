#include "PolarizationLeakageConstraint.h"

namespace dp3::ddecal {

std::vector<ConstraintResult> PolarizationLeakageConstraint::Apply(
    SolutionSpan& solutions, double,
    [[maybe_unused]] std::ostream* stat_stream) {
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) == 4);  // 2x2 full jones solutions
  for (size_t ch = 0; ch < NChannelBlocks(); ++ch) {
    for (size_t ant = 0; ant < NAntennas(); ++ant) {
      for (size_t sub_solution = 0; sub_solution != NSubSolutions();
           ++sub_solution) {
        std::complex<double>* matrix = &solutions(ch, ant, sub_solution, 0);
        matrix[0] = 1.0;
        matrix[3] = 1.0;
      }
    }
  }

  return std::vector<ConstraintResult>();
}

}  // namespace dp3::ddecal
