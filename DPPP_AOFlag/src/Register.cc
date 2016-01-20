//# Register.cc: Register steps in DPPP
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
//# $Id: AOFlaggerStep.cc 31423 2015-04-03 14:06:21Z dijkema $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP_AOFlag/Register.h>
#include <DPPP_AOFlag/AOFlaggerStep.h>
//#include <DPPP_AOFlag/SlidingFlagger.h>
#include <DPPP/DPRun.h>

// Define the function to make the AOFlaggerStep 'constructor' known.
// Its suffix must be the (lowercase) name of the package (library).
// Also make the SlidingFlagger known.
void register_aoflag()
{
  LOFAR::DPPP::DPRun::registerStepCtor ("aoflag",
                                        LOFAR::DPPP::AOFlaggerStep::makeStep);
  //  LOFAR::DPPP::DPRun::registerStepCtor ("aoflag.sliding",
  //                                        LOFAR::DPPP::SlidingFlagger::makeStep);
}
