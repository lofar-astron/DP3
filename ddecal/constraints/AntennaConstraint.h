// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ANTENNA_CONSTRAINT_H_
#define DP3_DDECAL_ANTENNA_CONSTRAINT_H_

#include "Constraint.h"

#include <set>
#include <vector>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmasked_view.hpp>
#include <xtensor/xview.hpp>

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
      SolutionSpan& solutions, [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) override {
    const size_t n_channels = solutions.shape(0);
    const size_t n_directions = solutions.shape(2);
    const size_t n_polarizations = solutions.shape(3);
    for (size_t ch = 0; ch < n_channels; ++ch) {
      for (const std::set<size_t>& antenna_set : antenna_sets_) {
        xt::xtensor<dcomplex, 2> solution_per_set =
            xt::zeros<dcomplex>({n_directions, n_polarizations});
        xt::xtensor<size_t, 2> solution_count_per_set =
            xt::zeros<size_t>({n_directions, n_polarizations});

        // Declaring these variables outside the loop allows reusing memory.
        xt::xtensor<dcomplex, 2> solution_set;
        xt::xtensor<std::size_t, 2> solution_set_is_finite;

        // Calculate the sum of solutions over the set of stations
        for (size_t antenna_index : antenna_set) {
          solution_set =
              xt::view(solutions, ch, antenna_index, xt::all(), xt::all());
          solution_set_is_finite = xt::isfinite(solution_set);
          xt::masked_view(solution_set, !solution_set_is_finite)
              .fill(std::complex<double>(0.0, 0.0));

          solution_per_set += solution_set;
          solution_count_per_set += solution_set_is_finite;
        }

        // Divide by nr of core stations to get the mean solution
        solution_per_set /= xt::cast<double>(solution_count_per_set);

        // Assign all core stations to the mean solution
        for (size_t antenna_index : antenna_set) {
          xt::view(solutions, ch, antenna_index, xt::all(), xt::all()) =
              solution_per_set;
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
