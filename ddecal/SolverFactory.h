// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_SOLVERFACTORY_H
#define DP3_SOLVERFACTORY_H

#include "../base/CalType.h"
#include "../base/Direction.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace ddecal {

class Settings;
class SolverBase;

std::unique_ptr<SolverBase> CreateSolver(const Settings& settings,
                                         const common::ParameterSet& parset,
                                         const std::string& prefix);

/**
 * Initializes all constraints for a given solver.
 *
 * While DPInfo can distinguish between used and unused antennas, the solvers
 * work with the used antennas only. antenna_positions and antenna_names thus
 * only contain used antennas. Their ordering should match the ordering while
 * solving, when using Constraint::Apply.
 *
 * @param solver A solver that may have constraints.
 * @param settings Includes settings for the constraints.
 * @param antenna_positions For each antenna, the position.
 * @param antenna_names For each antenna, the name.
 * @param source_positions For each direction, the source position in ra/dec.
 */
void InitializeSolverConstraints(
    SolverBase& solver, const Settings& settings,
    const std::vector<std::array<double, 3>>& antenna_positions,
    const std::vector<std::string>& antenna_names,
    const std::vector<size_t>& solutions_per_direction,
    const std::vector<base::Direction>& source_positions,
    const std::vector<double>& frequencies);

}  // namespace ddecal
}  // namespace dp3

#endif  // DP3_SOLVERFACTORY_H
