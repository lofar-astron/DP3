// pydp3.cc: python bindings to create python step
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include <casacore/measures/Measures/MDirection.h>

#include <dp3/base/DP3.h>
#include <dp3/steps/Step.h>
#include <dp3/common/Fields.h>

#include "../steps/InputStep.h"
#include "PyStep.h"
#include "utils.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;
using dp3::steps::Step;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapDpInfo(py::module &m);  // Defined in pydpinfo.cc

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

class PublicStep : public Step {
 public:
  using Step::info;
  using Step::infoOut;
  using Step::showTimings;
  using Step::updateInfo;
};

PYBIND11_MODULE(pydp3, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  register_cube<float>(m, "Cube_float");
  register_cube<bool>(m, "Cube_bool");
  register_cube<casacore::Complex>(m, "Cube_complex_float");
  register_matrix<double>(m, "Matrix_double");
  register_vector<double>(m, "Vector_double");
  register_vector<int>(m, "Vector_int");

  WrapDpInfo(m);

  py::class_<ostream_wrapper>(m, "ostream")
      .def("write", &ostream_wrapper::write);

  m.def("make_step", &dp3::base::MakeSingleStep);
  m.def("make_main_steps",
        [](const common::ParameterSet &parset) -> std::shared_ptr<steps::Step> {
          return dp3::base::MakeMainSteps(parset);
        });
  m.def("get_chain_required_fields", &dp3::base::GetChainRequiredFields);

  py::enum_<dp3::steps::Step::MsType>(m, "MsType")
      .value("regular", dp3::steps::Step::MsType::kRegular)
      .value("bda", dp3::steps::Step::MsType::kBda);

  py::class_<dp3::steps::Step, std::shared_ptr<dp3::steps::Step>, PyStep>(
      m, "Step")
      .def(py::init<>())
      .def(
          "show", [](const dp3::steps::Step &step) { step.show(std::cout); },
          "Show step summary. When running embedded in DP3, sys.stdout will be "
          "redirected to DP3's output stream. Overrides of show() should use "
          "print() "
          "to output a summary of the step.")
      .def(
          "show_timings",
          [](const dp3::steps::Step &step, double elapsed_time) {
            std::stringstream ss;
            step.showTimings(ss, elapsed_time);
            return ss.str();
          },
          "Show the processing time of current step. Also provides the "
          "percentage of time the current step took compared to the total "
          "elapsed time given as input (in seconds).")
      .def("__str__",
           [](const dp3::steps::Step &step) {
             std::stringstream ss;
             step.show(ss);
             return ss.str();
           })
      // Step::updateInfo is protected
      // Prepending by an underscore to indicate that this mehod is not
      // supposed to be called directly
      .def("_update_info", &PublicStep::updateInfo, "Handle metadata")
      .def("set_info", &Step::setInfo, py::return_value_policy::reference,
           "Set info object. This will call _update_info() for this step and "
           "all next steps")
      // Step::getInfoIn(), Step::getInfoOut() and Step::getInfo()
      // return const references.
      // Unfortunately there is no way to enforce constness in Python.
      // Also there is no guarantee that the Step will outlive the returned
      // info object. Although not the most efficient, the safest way
      // is returning a copy here.
      .def_property_readonly(
          "info_in", &Step::getInfoIn, py::return_value_policy::copy,
          "Get a copy of the info object containing metadata of the input")
      .def_property_readonly(
          "info_out", &Step::getInfoOut, py::return_value_policy::copy,
          "Get a copy of the info object containing metadata of the output")
      // Since a python step needs to be able to adjust its info_out in its
      // _update_info() override, _info_out returns a modifyable reference.
      .def_property_readonly("_info_out", &PublicStep::infoOut,
                             py::return_value_policy::reference,
                             "Get a modifyable reference to info object "
                             "containing metadata of the output")
      // Legacy version of info_out
      .def_property_readonly("info", &PublicStep::getInfo,
                             py::return_value_policy::reference,
                             "Get a modifyable reference to info object "
                             "containing metadata of the output")
      .def(
          "process",
          [](dp3::steps::Step &step, const base::DPBuffer &buffer) {
            return step.process(buffer);
          },
          "process buffer")
      .def("finish", &Step::finish,
           "Finish processing (nextstep->finish will be called automatically")
      .def("get_next_step", &Step::getNextStep,
           "Get a reference to the next step")
      .def("set_next_step", &Step::setNextStep,
           "Set the step that follows the current step")
      .def("get_required_fields", &Step::getRequiredFields,
           "Get the fields required by the current step")
      .def("get_provided_fields", &Step::getProvidedFields,
           "Get the fields provided by the current step");

  py::class_<DPBuffer, std::shared_ptr<DPBuffer>>(m, "DPBuffer")
      .def(py::init<>())
      .def(py::init<double, double>(),
           "Constructor with arguments: time, exposure")
      .def("get_data",
           (casacore::Cube<casacore::Complex> &
            (DPBuffer::*)())(&DPBuffer::GetCasacoreData),
           "Get data "
           "buffer that can be read as numpy array. Shape is "
           "(nr baselines, nr channels, nr polarizations)")
      .def("get_weights", (casacore::Cube<float> &
                           (DPBuffer::*)())(&DPBuffer::GetCasacoreWeights))
      .def("get_flags", (casacore::Cube<bool> &
                         (DPBuffer::*)())(&DPBuffer::GetCasacoreFlags))
      .def("get_uvw", (casacore::Matrix<double> &
                       (DPBuffer::*)())(&DPBuffer::GetCasacoreUvw))
      .def("get_exposure", &DPBuffer::getExposure,
           "Get the exposure of this buffer")
      .def(
          "set_data",
          [](DPBuffer &self,
             py::array_t<std::complex<float>, py::array::c_style> &numpy_data) {
            if (numpy_data.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            casacore::Cube<casacore::Complex> data(
                numpy_data.shape(2), numpy_data.shape(1), numpy_data.shape(0));
            std::copy_n(numpy_data.data(), data.size(), data.data());
            self.setData(data);
          },
          "set buffer data from a numpy array of complex float type.")
      .def(
          "set_weights",
          [](DPBuffer &self,
             py::array_t<float, py::array::c_style> &numpy_weights) {
            if (numpy_weights.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            casacore::Cube<float> weights(numpy_weights.shape(2),
                                          numpy_weights.shape(1),
                                          numpy_weights.shape(0));
            std::copy_n(numpy_weights.data(), weights.size(), weights.data());
            self.setWeights(weights);
          },
          "set buffer weights from a numpy array of float type.")
      .def(
          "set_flags",
          [](DPBuffer &self,
             py::array_t<bool, py::array::c_style> &numpy_flags) {
            if (numpy_flags.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            casacore::Cube<bool> flags(numpy_flags.shape(2),
                                       numpy_flags.shape(1),
                                       numpy_flags.shape(0));
            std::copy_n(numpy_flags.data(), flags.size(), flags.data());

            self.setFlags(flags);
          },
          "set buffer flags from a numpy array of boolean type.")
      .def(
          "set_uvw",
          [](DPBuffer &self,
             py::array_t<double, py::array::c_style> &numpy_uvw) {
            if (numpy_uvw.ndim() != 2) {
              throw std::runtime_error(
                  "Provided array should have two dimensions (n_baselines, "
                  "3).");
            }
            if (numpy_uvw.shape(1) != 3) {
              throw std::runtime_error(
                  "Each baseline should have 3 uvw values.");
            }
            casacore::Matrix<double> uvw(3, numpy_uvw.shape(0));

            std::copy_n(numpy_uvw.data(), uvw.size(), uvw.data());

            self.setUVW(uvw);
          },
          "set buffer uvw from a numpy array of double type.")
      .def("set_exposure", &DPBuffer::setExposure,
           "Set the exposure of this buffer")
      .def("get_time", &DPBuffer::getTime, "Get the time of this buffer")
      .def("set_time", &DPBuffer::setTime, "Set the time of this buffer");

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
    auto key = kv.first.cast<string>();
    auto value = kv.second[py::int_(0)];
    fields.def_property_readonly_static(
        key.c_str(), [value](py::object fields) { return fields(value); });
  }
}

}  // namespace pythondp3
}  // namespace dp3
