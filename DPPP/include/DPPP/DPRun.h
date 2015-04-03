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
#include <Common/ParameterSet.h>

#include <map>

namespace LOFAR {
  namespace DPPP {

    // @ingroup NDPPP

    // This class contains a single static function that creates and executes
    // the steps defined in the parset file.
    // The parset file is documented on the LOFAR wiki.

    class DPRun
    {
    public:
      // Define the function to create a step from the given parameterset.
      typedef DPStep::ShPtr StepCtor (DPInput*, const ParameterSet&,
                                      const std::string& prefix);

      // Add a function creating a DPStep to the map.
      static void registerStepCtor (const std::string&, StepCtor*);

      // Create a step object from the given parameters.
      // It looks up the step type in theirStepMap. If not found, it will
      // try to load a shared library with that name and execute the
      // register function in it.
      static StepCtor* findStepCtor (const std::string& type);

      // Execute the steps defined in the parset file.
      // Possible parameters given at the command line are taken into account.
      static void execute (const std::string& parsetName,
                           int argc=0, char* argv[] = 0);

    private:
      // Create the step objects.
      // It fills DPInfo object and the name of the MS being written.
      static DPStep::ShPtr makeSteps (const ParameterSet& parset,
                                      std::string& msName);

      // The map to create a step object from its type name.
      static std::map<std::string, StepCtor*> theirStepMap;
    };

  } //# end namespace
}

#endif
