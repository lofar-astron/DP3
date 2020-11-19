// Register.cc: Register steps in DPPP
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Register.h"
#include "AOFlaggerStep.h"

#include "../DPPP/DPRun.h"

// Define the function to make the AOFlaggerStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
// Also make the SlidingFlagger known.
void register_aoflagger() {
  DP3::DPPP::DPRun::registerStepCtor("aoflagger",
                                     DP3::DPPP::AOFlaggerStep::makeStep);
}
