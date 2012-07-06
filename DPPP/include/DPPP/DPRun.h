//# DPRun.h: Class to run steps like averaging and flagging on an MS
//# Copyright (C) 2010
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
//# $Id$
//#
//# @author Ger van Diepen

#ifndef DPPP_DPRUN_H
#define DPPP_DPRUN_H

// @file
// @brief Class to run steps like averaging and flagging on an MS

#include <lofar_config.h>
#include <DPPP/DPStep.h>
#include <DPPP/ParSet.h>

namespace LOFAR {
  namespace DPPP {

    // @ingroup NDPPP

    // This class contains a single static function that creates and executes
    // the steps defined in the parset file.
    // The parset file is documented on the LOFAR wiki.

    class DPRun
    {
    public:
      // Execute the stps defined in the parset file.
      static void execute (const std::string& parsetName);

    private:
      // Create the step objects.
      // It fills DPInfo object and the name of the MS being written.
      static DPStep::ShPtr makeSteps (const ParSet& parset,
                                      std::string& msName);
    };

  } //# end namespace
}

#endif
