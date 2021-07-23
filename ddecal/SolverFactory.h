// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_SOLVERFACTORY_H
#define DP3_SOLVERFACTORY_H

#include "../base/CalType.h"

#include <memory>
#include <string>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace ddecal {

class Settings;
class RegularSolverBase;
class BdaSolverBase;

std::unique_ptr<RegularSolverBase> CreateRegularSolver(
    const Settings& settings, const common::ParameterSet& parset,
    const std::string& prefix);

std::unique_ptr<BdaSolverBase> CreateBdaSolver(
    const Settings& settings, const common::ParameterSet& parset,
    const std::string& prefix);

}  // namespace ddecal
}  // namespace dp3

#endif  // DP3_SOLVERFACTORY_H
