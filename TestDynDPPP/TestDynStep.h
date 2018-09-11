//# TestDynStep.h: Test of a dynamically loaded DPPP step
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

#ifndef TESTDYNDPPP_TESTDYNSTEP_H
#define TESTDYNDPPP_TESTDYNSTEP_H

// @file
// @brief Test of a dynamically loaded DPPP step

#include "../DPPP/DPStep.h"
#include "../DPPP/Averager.h"
#include "../DPPP/DPInput.h"
#include "../Common/ParameterSet.h"

namespace DP3 {
  namespace DPPP {
    // @ingroup NDPPP

    // This class is a test (and an example) of a DPStep loaded
    // dynamically from a shared library.
    // To make test life easy it uses the Averager class underneath.

    class TestDynStep: public Averager
    {
    public:
      TestDynStep (DPInput*, const ParameterSet&, const std::string&);
      virtual ~TestDynStep();
      static DPStep::ShPtr makeStep (DPInput*, const ParameterSet&,
                                     const std::string&);
    };


  }
}

// Define the function (without name mangling) to register the 'constructor'.
extern "C"
{
  void register_testdyndppp();
}

#endif
