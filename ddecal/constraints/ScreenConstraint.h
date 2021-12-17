// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SCREEN_CONSTRAINT_H_
#define DP3_DDECAL_SCREEN_CONSTRAINT_H_

#include "../../base/PhaseFitter.h"
#include "../../base/Direction.h"

#include "Constraint.h"
#include "PiercePoint.h"
#include "KLFitter.h"

#include "../../common/ParameterSet.h"

#include <cmath>
#include <complex>
#include <memory>
#include <ostream>
#include <vector>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace ddecal {

class ScreenConstraint final : public Constraint {
  static constexpr double kPhaseToTec = 1.0 / 8.4479745e9;
  static constexpr double kTecToPhase = 8.4479745e9;
  static constexpr size_t kMaxIterations =
      30;  // number of iterations to store in debug mode

 public:
  ScreenConstraint(const common::ParameterSet& parset, const string& prefix);

  void Initialize(size_t nAntennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) override;
  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions, double time,
      std::ostream* statStream) override;

  void SetCoreAntennas(const std::set<size_t>& core_antennas);
  void InitPiercePoints(const std::vector<std::array<double, 3>>& antenna_pos,
                        const std::vector<base::Direction>& source_directions);

  const std::vector<size_t>& GetCoreAntennas() const { return core_antennas_; }
  const std::vector<std::vector<PiercePoint>>& GetPiercePoints() const {
    return pierce_points_;
  }

 private:
  void SetTime(double time);
  void CalculatePiercepoints();
  void GetPpValue(const std::vector<std::vector<std::complex<double>>>&,
                  size_t solution_index, size_t direction_index,
                  double& avg_tec, double& error) const;

  std::vector<double> frequencies_;
  std::vector<double> previous_solution_;
  std::vector<double> iter_phases_;
  /// antenna positions
  /// source positions
  /// measures instance ofzo
  std::vector<std::vector<PiercePoint>>
      pierce_points_;  // temporary hold calculated piercepoints per antenna
  std::vector<KLFitter> screen_fitters_;
  std::vector<size_t> core_antennas_;
  std::vector<size_t> other_antennas_;  // has to be a vector for openmp looping
  double current_time_;
  double beta_;
  double height_;
  double order_;
  double r_diff_;
  std::string mode_;
  std::string avg_mode_;
  int debug_mode_;
  size_t iteration_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
