#ifndef DP3_DDECAL_TEC_DELAY_CONSTRAINT_H_
#define DP3_DDECAL_TEC_DELAY_CONSTRAINT_H_

#include <vector>

#include "Constraint.h"
#include "ConstraintResult.h"

#include "TecOffsetDelayFitting.h"

#include <vector>
#include <ostream>

namespace dp3::ddecal {

class TecOffsetDelayConstraint final : public Constraint {
 public:
  TecOffsetDelayConstraint(bool include_offset, size_t max_wraps,
                           bool do_phase_referencing);

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;

  void SetWeights(const std::vector<double>& weights) override;

  void Apply(SolutionSpan& solutions, double time) override;

  std::vector<ConstraintResult> GetResult() const override { return results_; };

 protected:
  struct ThreadData {
    std::vector<double> data;
    std::vector<double> weights;
  };
  std::vector<ThreadData> thread_data_;

  bool include_offset_ = false;
  bool do_phase_referencing_ = true;
  size_t max_wraps_ = 0;
  std::vector<double> weights_;
  double reference_frequency_ = 0.0;
  std::vector<double> relative_frequencies_;
  std::vector<ConstraintResult> results_;
};

}  // namespace dp3::ddecal

#endif
