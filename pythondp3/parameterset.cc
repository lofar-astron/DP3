// parameterset.cc: parameterset python bindings
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>

#include "../common/ParameterSet.h"

using dp3::common::ParameterSet;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

PYBIND11_MODULE(parameterset, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  py::class_<ParameterSet,
             std::shared_ptr<ParameterSet>  // holder type
             >(m, "ParameterSet")
      .def(py::init<>())
      .def("getBool", (bool (ParameterSet::*)(const std::string&) const) &
                          ParameterSet::getBool)
      .def("getBool", (bool (ParameterSet::*)(const std::string&, bool) const) &
                          ParameterSet::getBool)
      .def("getInt", (int (ParameterSet::*)(const std::string&) const) &
                         ParameterSet::getInt)
      .def("getInt", (int (ParameterSet::*)(const std::string&, int) const) &
                         ParameterSet::getInt)
      .def("getFloat", (float (ParameterSet::*)(const std::string&) const) &
                           ParameterSet::getFloat)
      .def("getFloat",
           (float (ParameterSet::*)(const std::string&, float) const) &
               ParameterSet::getFloat)
      .def("getDouble", (double (ParameterSet::*)(const std::string&) const) &
                            ParameterSet::getDouble)
      .def("getDouble",
           (double (ParameterSet::*)(const std::string&, double) const) &
               ParameterSet::getDouble)
      .def("getString",
           (std::string(ParameterSet::*)(const std::string&) const) &
               ParameterSet::getString)
      .def("getString", (std::string(ParameterSet::*)(
                            const std::string&, const std::string&) const) &
                            ParameterSet::getString);
}

}  // namespace pythondp3
}  // namespace dp3
