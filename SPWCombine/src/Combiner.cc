//# Copyright (C) 2007
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
//# @author Adriaan Renting

#include <lofar_config.h>
#include <Common/lofar_iostream.h>
#include <Common/SystemUtil.h>
#include <Common/Exceptions.h>
#include <PLC/ACCmain.h>
#include <SPWCombine/CombinerProcessControl.h>

using namespace LOFAR;

// Use a terminate handler that can produce a backtrace.
Exception::TerminateHandler t(Exception::terminate);

int main(int argc, char *argv[])
{
  try {
    INIT_LOGGER(LOFAR::basename(argv[0]));
    LOFAR::CS1::CombinerProcessControl myProcess;
    return LOFAR::ACC::PLC::ACCmain(argc, argv, &myProcess);
  } catch(Exception& ex) {
    cerr << ex << endl;
    return 1;
  }
  return 0;
}
