// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include <pybind11/pybind11.h>

#include <dp3/common/Fields.h>

using dp3::common::Fields;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapFields(py::module &m) {
  py::class_<Fields, std::shared_ptr<Fields>> fields(m, "Fields");
  fields.def(py::init<>())
      .def(py::init<Fields>())
      .def(py::init<Fields::Single>())
      .def_property_readonly("data", &Fields::Data)
      .def_property_readonly("flags", &Fields::Flags)
      .def_property_readonly("uvw", &Fields::Uvw)
      .def_property_readonly("weights", &Fields::Weights)
      .def(
          "update_requirements",
          [](Fields &self, const Fields &a, const Fields &b) {
            return self.UpdateRequirements(a, b);
          },
          "Updates the current object's Fields based on a step's required and "
          "provided Fields")
      .def("__str__",
           [](const Fields &a) {
             std::stringstream ss;
             ss << a;
             return ss.str();
           })
      .def("__eq__", [](const Fields &a, const Fields &b) { return a == b; })
      .def("__neq__", [](const Fields &a, const Fields &b) { return a != b; })
      .def("__or__", [](const Fields &a, const Fields &b) { return a | b; })
      .def("__ior__", [](Fields &a, const Fields &b) {
        a |= b;
        return a;
      });

  py::enum_<Fields::Single> fields_enum(fields, "Single");
  fields_enum.value("DATA", Fields::Single::kData)
      .value("FLAGS", Fields::Single::kFlags)
      .value("WEIGHTS", Fields::Single::kWeights)
      .value("UVW", Fields::Single::kUvw);

  // Custom export_values(), adaptation of the export_values() function in
  // pybind11.h For each entry in the Fields.Single enum a static property is
  // added to the Fields class. The getter function returns a new instance that
  // can be modified freely without affecting the value returned by subsequent
  // calls to the getter. This allows the usage of, for example, Fields.DATA
  // instead of the more verbose Fields(Fields.Single.DATA)
  py::dict entries = fields_enum.attr("__entries");
  for (auto kv : entries) {
    auto key = kv.first.cast<std::string>();
    auto value = kv.second[py::int_(0)];
    fields.def_property_readonly_static(
        key.c_str(), [value](py::object fields) { return fields(value); });
  }
}

}  // namespace pythondp3
}  // namespace dp3
