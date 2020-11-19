// Register.cc: Register steps in DPPP
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Register.h"
#include "DDECal.h"

#include "../DPPP/DPRun.h"

// Define the function to make the DDECal 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
void register_ddecal() {
  DP3::DPPP::DPRun::registerStepCtor("ddecal", DP3::DPPP::DDECal::makeStep);
}
