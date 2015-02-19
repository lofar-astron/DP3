//# TestDyn.cc: Test of a dynamically loaded DPPP step
//# Copyright (C) 2015
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: Averager.h 30711 2015-01-15 14:36:28Z diepen $
//#
//# @author Ger van Diepen

// @file
// @brief Test of a dynamically loaded DPPP step

#include <lofar_config.h>
#include "TestDynStep.h"
#include <DPPP/DPBuffer.h>
#include <DPPP/DPRun.h>

namespace LOFAR {
  namespace DPPP {

    TestDynStep::TestDynStep (DPInput* input, const ParameterSet& pset,
                              const std::string& prefix)
      : Averager (input, pset, prefix)
    {}

    TestDynStep::~TestDynStep()
    {}

    DPStep::ShPtr TestDynStep::makeStep (DPInput* input,
                                         const ParameterSet& pset,
                                         const std::string& prefix)
      { return DPStep::ShPtr(new TestDynStep(input, pset, prefix)); }

  }
}

// Define the function to make the TestDynStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
void register_testdyndppp()
{
  LOFAR::DPPP::DPRun::registerStepCtor ("TestDynDPPP",
                                        LOFAR::DPPP::TestDynStep::makeStep);
}
