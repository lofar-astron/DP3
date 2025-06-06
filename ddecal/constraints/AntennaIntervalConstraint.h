// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_ANTENNA_INTERVAL_CONSTRAINT_H_
#define DP3_DDECAL_ANTENNA_INTERVAL_CONSTRAINT_H_

#include "Constraint.h"

#include <vector>

namespace dp3::ddecal {

/**
 * This constraint allows using different intervals for different antennas.
 * It does this by averaging solutions of the same antenna but from different
 * sub-intervals together.
 *
 * Ionospheric effects vary more strongly in time and frequency for more distant
 * stations (relative to a central reference station). Thus, by using shorter
 * time intervals and less frequency smoothing for more distant groups of
 * stations, these fast variations can be captured without adding too much
 * additional free parameters for more central stations. Antenna intervals and
 * smoothing depend on the expected variation per group of antenna. For LOFAR,
 * CS can typically have at least 8 times longer intervals compared to RS.
 *
 * This constraint works together with the "direction-dependent intervals"
 * setting. If only antenna-dependent intervals are desired, the
 * solutions-per-direction setting should be set to the maximum antenna
 * averaging factor. E.g., with antenna_averaging_factors = {3, 6, 1}, the
 * solutions_per_direction is set to 6 for all directions. This will use a
 * time interval of 1/2, 1 and 1/6 solution interval for the antennas 0, 1 and
 * 2, respectively.
 *
 * Directional and antenna dependence can be combined, in which case the two
 * work cummulatively. That is, with a setting of antenna_averaging_factors =
 * {3, 6, 1} and solutions_per_direction = {6, 6, 1, 1} for direction 0 and 1,
 * the solver will divide the solution interval into 6 subsolutions. For antenna
 * 0, every 3 of those will be averaged together. For directions 2 and 3,
 * there is only 1 subsolution available, so no averaging is performed for any
 * of the antennas.
 */
class AntennaIntervalConstraint final : public Constraint {
 public:
  AntennaIntervalConstraint(
      const std::vector<size_t>& antenna_averaging_factors)
      : antenna_averaging_factors_(antenna_averaging_factors) {}

  void Initialize(size_t n_antennas,
                  const std::vector<uint32_t>& solutions_per_direction,
                  const std::vector<double>& frequencies) final;

  const std::vector<size_t>& GetIntervalsPerAntenna() const {
    return antenna_averaging_factors_;
  }

  std::vector<Constraint::Result> Apply(
      SolutionSpan& solutions, [[maybe_unused]] double time,
      [[maybe_unused]] std::ostream* stat_stream) final;

 private:
  std::vector<size_t> antenna_averaging_factors_;
};

}  // namespace dp3::ddecal

#endif
