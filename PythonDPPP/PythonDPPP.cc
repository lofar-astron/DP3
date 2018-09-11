//# PythonDPPP.cc: Python base class for a DPStep in python
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

#include "DPStepBase.h"

#include <casacore/python/Converters/PycExcp.h>
#include <casacore/python/Converters/PycBasicData.h>
#include <casacore/python/Converters/PycValueHolder.h>
#include <casacore/python/Converters/PycRecord.h>

#include <boost/python.hpp>
#include <boost/python/args.hpp>

using namespace boost::python;

// The C++ PythonStep must be able to call functions in Python.
// But the Python functions must be able to call the C++ functions
// (e.g., to get data, flags, etc.)
// All communication goes through DPStepBase
// Let a python step be created by a static function creating a
// PythonWorker doing the actual work. Its pointer is kept in a static
// and thereafter used to create a PythonStep.

namespace DP3 {
  namespace DPPP {

    // Define the interface for the PythonStep C++ functions callable
    // from python.
    // Note that the python functions called from C++ are done using
    // the boost::python attr function.
    void dpstepbase()
    {
      class_<DPStepBase> ("_DPStepBase")
        .def (init<>())
        .def ("_getData", &DPStepBase::_getData,
              "Get the visibility data into the given array",
              (boost::python::arg("value")))
        .def ("_getFlags", &DPStepBase::_getFlags,
              "Get the flags into the given array",
              (boost::python::arg("value")))
        .def ("_getWeights", &DPStepBase::_getWeights,
              "Get the weights into the given array",
              (boost::python::arg("value")))
        .def ("_getUVW", &DPStepBase::_getUVW,
              "Get the UVW coordinates into the given array",
              (boost::python::arg("value")))
        .def ("_getModelData", &DPStepBase::_getModelData,
              "Get the model data into the given array",
              (boost::python::arg("value")))
        .def ("_processNext", &DPStepBase::_processNext,
              "Process the next step in the DPPP run",
              (boost::python::arg("values")))
        ;
    }

  }
}


// Define the python module itself.
BOOST_PYTHON_MODULE(_pythondppp)
{
  casacore::python::register_convert_excp();
  casacore::python::register_convert_basicdata();
  casacore::python::register_convert_casa_valueholder();
  casacore::python::register_convert_casa_record();

  DP3::DPPP::dpstepbase();
}
