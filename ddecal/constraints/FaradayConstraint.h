#ifndef DP3_FARADAY_ROTATION_CONSTRAINT_H_
#define DP3_FARADAY_ROTATION_CONSTRAINT_H_

#include "Constraint.h"

#include "../../base/CalType.h"
#include "../../common/PhaseLineFitter.h"

#include <vector>
#include <ostream>

namespace dp3::ddecal {

class FaradayConstraint final : public Constraint {
 public:
  /**
   * The diagonal solution type may be set to kRotational
   * to fit Faraday rotation without a diagonal.
   */
  FaradayConstraint(base::CalType diagonal_solution_type)
      : diagonal_solution_type_(diagonal_solution_type) {}

  std::vector<Result> Apply(SolutionSpan& solutions, double time,
                            std::ostream* stat_stream) final;

  void Initialize(size_t n_antennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) final;

  void SetWeights(const std::vector<double>& weights) final;

  void SetSubSolutionWeights(
      const std::vector<std::vector<double>>& sub_solution_weights) final;

 private:
  void PerformFit(SolutionSpan& solutions, size_t sub_solution, size_t antenna,
                  std::vector<common::phase_fitting::FitSample>& scratch_space);

  std::vector<Constraint::Result> results_;
  std::vector<std::vector<double>> sub_solution_weights_;
  std::vector<double> frequencies_;
  common::phase_fitting::SlopeFitRange fit_range_;
  base::CalType diagonal_solution_type_;
};

}  // namespace dp3::ddecal

#endif
