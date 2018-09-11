//# DPStepBase.cc: Python base class for a DPStep in python
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
//# $Id: pyparameterset.cc 23074 2012-12-03 07:51:29Z diepen $

#include <lofar_config.h>
#include <PythonDPPP/DPStepBase.h>

// The C++ PythonStep must be able to call functions in Python.
// But the Python functions must be able to call the C++ functions
// (e.g., to get data, flags, etc.)
// All communication goes through DPStepBase
// Let a python step be created by a static function creating a
// PythonWorker doing the actual work. Its pointer is kept in a static
// and thereafter used to create a PythonStep.

namespace DP3 {
  namespace DPPP {

    // Define the static.
    DPStepBase* DPStepBase::theirPtr=0;

  }
}
