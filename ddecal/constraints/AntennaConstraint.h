// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ANTENNA_CONSTRAINT_H_
#define DP3_DDECAL_ANTENNA_CONSTRAINT_H_

#include "Constraint.h"

#include <set>
#include <vector>

namespace dp3 {
namespace ddecal {

/**
 * @brief This constraint averages the solutions of several groups of antennas,
 * so that antennas within the same group have equal solutions.
 *
 * The DDE solver can use this constraint e.g. to average the solutions of
 * the core antennas. Core antennas are determined by a given maximum distance
 * from a reference antenna. The reference antenna is by default the first
 * antenna. This constraint is meant to force all core stations to
 * have the same solution, thereby decreasing the noise in their solutions.
 */
class AntennaConstraint final : public Constraint {
 public:
  void SetAntennaSets(std::vector<std::set<size_t>>&& antenna_sets) {
    antenna_sets_ = std::move(antenna_sets);
  }

  const std::vector<std::set<size_t>>& GetAntennaSets() const {
    return antenna_sets_;
  }

  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<dcomplex>>& solutions,
      [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) override {
    // nGains is nPol x nSolutions (i.e., nr of gains per antenna)
    size_t n_gains = solutions.front().size() / NAntennas();
    std::vector<dcomplex> solution_per_set(n_gains);
    std::vector<size_t> solution_count_per_set(n_gains);
    for (size_t ch = 0; ch < solutions.size(); ++ch) {
      for (const std::set<size_t>& antenna_set : antenna_sets_) {
        solution_per_set.assign(n_gains, 0.0);
        solution_count_per_set.assign(n_gains, 0);
        // Calculate the sum of solutions over the set of stations
        for (size_t antenna_index : antenna_set) {
          size_t start_index = antenna_index * n_gains;
          for (size_t sol_index = 0; sol_index != n_gains; ++sol_index) {
            dcomplex value = solutions[ch][start_index + sol_index];
            if (isfinite(value)) {
              solution_per_set[sol_index] += value;
              ++solution_count_per_set[sol_index];
            }
          }
        }

        // Divide by nr of core stations to get the mean solution
        for (size_t sol_index = 0; sol_index != n_gains; ++sol_index)
          solution_per_set[sol_index] /= solution_count_per_set[sol_index];

        // Assign all core stations to the mean solution
        for (size_t antenna_index : antenna_set) {
          size_t start_index = antenna_index * n_gains;
          for (size_t sol_index = 0; sol_index != n_gains; ++sol_index)
            solutions[ch][start_index + sol_index] =
                solution_per_set[sol_index];
        }
      }
    }
    return std::vector<Constraint::Result>();
  }

 private:
  std::vector<std::set<size_t>> antenna_sets_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
