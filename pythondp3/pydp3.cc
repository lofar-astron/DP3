// pydp3.cc: python bindings to create python step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PyStep.h"
#include "PyStepImpl.h"
#include <dp3/common/Fields.h>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <casacore/measures/Measures/MDirection.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

template <typename T>
void register_cube(py::module &m, const char *name) {
  py::class_<casacore::Cube<T>>(m, name, py::buffer_protocol())
      .def_buffer([](casacore::Cube<T> &cube) -> py::buffer_info {
        return py::buffer_info(
            cube.data(),                        /* Pointer to buffer */
            sizeof(T),                          /* Size of one scalar */
            py::format_descriptor<T>::format(), /* Python struct-style format
                                                   descriptor */
            3,                                  /* Number of dimensions */
            {cube.shape()[2], cube.shape()[1],
             cube.shape()[0]}, /* Buffer dimensions */
            {sizeof(T) * cube.shape()[1] *
                 cube.shape()[0], /* Strides (in bytes) for each index */
             sizeof(T) * cube.shape()[0], sizeof(T)});
      });
}

template <typename T>
void register_matrix(py::module &m, const char *name) {
  py::class_<casacore::Matrix<T>>(m, name, py::buffer_protocol())
      .def_buffer([](casacore::Matrix<T> &matrix) -> py::buffer_info {
        return py::buffer_info(
            matrix.data(),                      /* Pointer to buffer */
            sizeof(T),                          /* Size of one scalar */
            py::format_descriptor<T>::format(), /* Python struct-style format
                                                   descriptor */
            2,                                  /* Number of dimensions */
            {matrix.shape()[1], matrix.shape()[0]}, /* Buffer dimensions */
            {sizeof(T) *
                 matrix.shape()[0], /* Strides (in bytes) for each index */
             sizeof(T)});
      });
}

template <typename T>
void register_vector(py::module &m, const char *name) {
  py::class_<casacore::Vector<T>>(m, name, py::buffer_protocol())
      .def_buffer([](casacore::Vector<T> &vector) -> py::buffer_info {
        return py::buffer_info(
            vector.data(),                      /* Pointer to buffer */
            sizeof(T),                          /* Size of one scalar */
            py::format_descriptor<T>::format(), /* Python struct-style format
                                                   descriptor */
            1,                                  /* Number of dimensions */
            {vector.shape()[0]},                /* Buffer dimensions */
            {sizeof(T)} /* Strides (in bytes) for each index */
        );
      });
}

PYBIND11_MODULE(pydp3, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  register_cube<float>(m, "Cube_float");
  register_cube<bool>(m, "Cube_bool");
  register_cube<casacore::Complex>(m, "Cube_complex_float");
  register_matrix<double>(m, "Matrix_double");
  register_vector<double>(m, "Vector_double");
  register_vector<int>(m, "Vector_int");

  py::class_<ostream_wrapper>(m, "ostream")
      .def("write", &ostream_wrapper::write);

  py::class_<StepWrapper, PyStepImpl /* <--- trampoline*/
             >(m, "Step")
      .def(py::init<>())
      .def("show", &StepWrapper::show,
           "Show step summary (stdout will be redirected to DPPP's output "
           "stream during this step)")
      .def("update_info", &StepWrapper::updateInfo, "Handle metadata")
      .def("info", &StepWrapper::info, py::return_value_policy::reference,
           "Get info object (read/write) with metadata")
      .def("finish", &StepWrapper::finish,
           "Finish processing (nextstep->finish will be called automatically")
      .def("get_next_step", &StepWrapper::get_next_step,
           py::return_value_policy::copy, "Get a reference to the next step")
      .def("process_next_step", &StepWrapper::process_next_step,
           "Process the next step")
      .def("get_count", &StepWrapper::get_count,
           "Get the number of time slots processed")
      .def("get_required_fields", &StepWrapper::getRequiredFields,
           "Get the fields required by current step")
      .def("get_provided_fields", &StepWrapper::getProvidedFields,
           "Get the fields provided by current step")
      .def_readwrite("fetch_uvw", &StepWrapper::m_fetch_uvw,
                     "Fill the UVW data in the buffer")
      .def_readwrite("fetch_weights", &StepWrapper::m_fetch_weights,
                     "Fill the weights data in the buffer");

  py::class_<DPBuffer, std::shared_ptr<DPBuffer>>(m, "DPBuffer")
      .def(py::init<>())
      .def("get_data",
           (casacore::Cube<casacore::Complex> &
            (DPBuffer::*)())(&DPBuffer::getData),
           "Get data "
           "buffer that can be read as numpy array. Shape is "
           "(nr baselines, nr channels, nr polarizations)")
      .def("get_weights",
           (casacore::Cube<float> & (DPBuffer::*)())(&DPBuffer::getWeights))
      .def("get_flags",
           (casacore::Cube<bool> & (DPBuffer::*)())(&DPBuffer::getFlags))
      .def("get_uvw",
           (casacore::Matrix<double> & (DPBuffer::*)())(&DPBuffer::getUVW))
      .def("get_exposure", &DPBuffer::getExposure,
           "Get the exposure of this buffer")
      .def("set_exposure", &DPBuffer::setExposure,
           "Set the exposure of this buffer")
      .def("get_time", &DPBuffer::getTime, "Get the time of this buffer")
      .def("set_time", &DPBuffer::setTime, "Set the time of this buffer");

  py::class_<DPInfo>(m, "DPInfo")
      .def(
          "antenna_names",
          [](DPInfo &self) -> py::array {
            // Convert casa vector of casa strings to std::vector of strings
            std::vector<std::string> names_casa = self.antennaNames();
            std::vector<std::string> names;
            for (size_t i = 0; i < names_casa.size(); ++i) {
              names.push_back(names_casa[i]);
            }
            py::array ret = py::cast(names);
            return ret;
          },
          "Get numpy array with antenna names (read only)")
      .def(
          "antenna_positions",
          [](DPInfo &self) -> py::array {
            // Convert vector of casa MPositions to std::vector of positions
            std::vector<casacore::MPosition> positions_casa = self.antennaPos();
            std::vector<std::array<double, 3>> positions;
            for (size_t i = 0; i < positions_casa.size(); ++i) {
              casacore::MVPosition position_mv = positions_casa[i].getValue();
              std::array<double, 3> position_array = {
                  position_mv(0), position_mv(1), position_mv(2)};
              positions.push_back(position_array);
            }
            py::array ret = py::cast(positions);
            return ret;
          },
          "Get a list of antenna positions in ITRF XYZ (read only)")
      .def("set_need_vis_data", &DPInfo::setNeedVisData,
           "Set whether data needs to be read before this step")
      .def("get_channel_frequencies", &DPInfo::chanFreqs,
           py::arg("baseline") = 0,
           "Get a list of channel frequencies (read only)")
      .def("get_antenna1", &DPInfo::getAnt1,
           "Get a list of first antenna numbers of correlations")
      .def("get_antenna2", &DPInfo::getAnt2,
           "Get a list of second antenna numbers of correlations")
      .def("nantenna", &DPInfo::nantenna, "Get the number of antennas")
      .def("nchan", &DPInfo::nchan, "Get the number of channels")
      .def("start_time", &DPInfo::startTime, "Get the start time")
      .def("time_interval", &DPInfo::timeInterval, "Get the time interval")
      .def("ntime", &DPInfo::ntime, "Get the total number of time slots")
      .def("ms_name", &DPInfo::msName, "Get name of measurement set");

  py::class_<Fields, std::shared_ptr<Fields>> fields(m, "Fields");
  fields.def(py::init<>())
      .def(py::init<Fields>())
      .def(py::init<Fields::Single>())
      .def_property_readonly("data", &Fields::Data)
      .def_property_readonly("flags", &Fields::Flags)
      .def_property_readonly("uvw", &Fields::Uvw)
      .def_property_readonly("weights", &Fields::Weights)
      .def_property_readonly("fullresflags", &Fields::FullResFlags)
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
      .value("FULLRESFLAGS", Fields::Single::kFullResFlags)
      .value("UVW", Fields::Single::kUvw);

  // Custom export_values()
  // For each entry in the Fields.Single enum a static property is added to the
  // Fields class. The getter function returns a new instance that can be
  // modified freely without affecting the value returned by subsequent calls to
  // the getter. This allows the usage of, for example, Fields.DATA instead of
  // the more verbose Fields(Fields.Single.DATA)
  py::dict entries = fields_enum.attr("__entries");
  for (auto kv : entries) {
    auto key = kv.first.cast<string>();
    auto value = kv.second[py::int_(0)];
    fields.def_property_readonly_static(
        key.c_str(), [value](py::object fields) { return fields(value); });
  }
}

void test() {}

}  // namespace pythondp3
}  // namespace dp3
