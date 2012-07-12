//# DPPP.cc: Program to execute steps like averaging and flagging on an MS
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

#include <lofar_config.h>
#include <DPPP/DPRun.h>
#include <DPPP/Package__Version.h>
#include <Common/LofarLogger.h>
#include <Common/SystemUtil.h>
#include <Common/Exception.h>
#include <iostream>
#include <stdexcept>

using namespace LOFAR::DPPP;
using namespace LOFAR;

// Define handler that tries to print a backtrace.
Exception::TerminateHandler t(Exception::terminate);

int main(int argc, char *argv[])
{
  try
  {
    TEST_SHOW_VERSION (argc, argv, DPPP);
    INIT_LOGGER(basename(string(argv[0])));
    // Get the name of the parset file.
    string parsetName("NDPPP.parset");
    if (argc > 1) {
      parsetName = argv[1];
    }
    // Execute the parset file.
    DPRun::execute (parsetName);
  } catch (LOFAR::Exception& err) {
    std::cerr << "LOFAR Exception detected: " << err << std::endl;
    return 1;
  }
  return 0;
}
