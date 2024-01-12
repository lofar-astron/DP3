// parameterset.cc: parameterset python bindings
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../common/ParameterSet.h"

using dp3::common::ParameterSet;
using dp3::common::ParameterValue;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

PYBIND11_MODULE(parameterset, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  py::class_<ParameterSet,
             std::shared_ptr<ParameterSet>  // holder type
             >(m, "ParameterSet")
      .def(py::init<>())
      .def(py::init<std::string>(),
           "Construct a ParameterSet from the contents of the filename given "
           "as argument.")
      .def("__str__",
           [](const ParameterSet &parameterset) {
             std::stringstream stream;
             stream << parameterset;
             return stream.str();
           })
      .def("add",
           (void(ParameterSet::*)(const std::string &, const std::string &)) &
               ParameterSet::add,
           "Add a key value pair")
      .def("get_bool", (bool(ParameterSet::*)(const std::string &) const) &
                           ParameterSet::getBool)
      .def("get_bool",
           (bool(ParameterSet::*)(const std::string &, bool) const) &
               ParameterSet::getBool)
      .def("get_int", (int(ParameterSet::*)(const std::string &) const) &
                          ParameterSet::getInt)
      .def("get_int", (int(ParameterSet::*)(const std::string &, int) const) &
                          ParameterSet::getInt)
      .def("get_float", (float(ParameterSet::*)(const std::string &) const) &
                            ParameterSet::getFloat)
      .def("get_float",
           (float(ParameterSet::*)(const std::string &, float) const) &
               ParameterSet::getFloat)
      .def("get_double", (double(ParameterSet::*)(const std::string &) const) &
                             ParameterSet::getDouble)
      .def("get_double",
           (double(ParameterSet::*)(const std::string &, double) const) &
               ParameterSet::getDouble)
      .def("get_string",
           (std::string(ParameterSet::*)(const std::string &) const) &
               ParameterSet::getString)
      .def("get_string", (std::string(ParameterSet::*)(
                             const std::string &, const std::string &) const) &
                             ParameterSet::getString)
      .def("__contains__", &ParameterSet::isDefined)
      .def("__getitem__", &ParameterSet::get);

  py::class_<ParameterValue,
             std::shared_ptr<ParameterValue>  // holder type
             >(m, "ParameterValue")
      .def("__str__", &ParameterValue::getString)
      .def("__int__", &ParameterValue::getInt)
      .def("__float__", &ParameterValue::getDouble)
      .def("__bool__", &ParameterValue::getBool)
      .def("__iter__", [](const ParameterValue &v) {
        return py::cast(v.getVector()).attr("__iter__")();
      });
}

}  // namespace pythondp3
}  // namespace dp3
