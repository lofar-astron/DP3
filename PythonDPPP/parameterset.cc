#include <pybind11/pybind11.h>

#include "../Common/ParameterSet.h"

namespace py = pybind11;

namespace DP3 {

PYBIND11_MODULE(parameterset, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  py::class_<ParameterSet,
             std::shared_ptr<ParameterSet>  // holder type
             >(m, "ParameterSet")
      .def(py::init<>())
      .def("getBool", (bool (ParameterSet::*)(const std::string& aKey) const) &
                          ParameterSet::getBool)
      .def("getInt", (int (ParameterSet::*)(const std::string& aKey) const) &
                         ParameterSet::getInt)
      .def("getFloat",
           (float (ParameterSet::*)(const std::string& aKey) const) &
               ParameterSet::getFloat)
      .def("getDouble",
           (double (ParameterSet::*)(const std::string& aKey) const) &
               ParameterSet::getDouble)
      .def("getString",
           (std::string(ParameterSet::*)(const std::string& aKey) const) &
               ParameterSet::getString);
}

}  // namespace DP3
