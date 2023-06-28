// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <casacore/casa/BasicSL/Complexfwd.h>

#include <dp3/base/DPBuffer.h>

using dp3::base::DPBuffer;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapDpBuffer(py::module &m) {
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
}

}  // namespace pythondp3
}  // namespace dp3
