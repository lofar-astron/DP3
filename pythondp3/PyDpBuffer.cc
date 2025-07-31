// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PyDpBuffer.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <aocommon/xt/span.h>

using dp3::base::DPBuffer;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

namespace {

/// Use Spans as datatype, and not the actual xtensor type, since pybind11
/// may copy the value. Copying a span does not hurt: It still points to the
/// same data.
template <typename T, size_t Dimensions>
void RegisterSpan(py::module& m, const char* name) {
  using SpanType = aocommon::xt::Span<T, Dimensions>;
  py::class_<SpanType>(m, name, py::buffer_protocol())
      .def_buffer([](SpanType& tensor) {
        constexpr std::size_t value_size = sizeof(T);
        auto strides = tensor.strides();
        for (auto& stride : strides) stride *= value_size;
        return py::buffer_info(tensor.data(), tensor.shape(), strides);
      });
}

/// Converts py::array_t::shape() into a shape suitable for DPBuffer.
std::array<std::size_t, 3> ConvertShape(const pybind11::ssize_t* shape) {
  return {static_cast<std::size_t>(shape[0]),
          static_cast<std::size_t>(shape[1]),
          static_cast<std::size_t>(shape[2])};
}

}  // namespace

void WrapDpBuffer(py::module& m) {
  RegisterSpan<std::complex<float>, 3>(m, "DPBuffer_Data");
  RegisterSpan<float, 3>(m, "DPBuffer_Weights");
  RegisterSpan<bool, 3>(m, "DPBuffer_Flags");
  RegisterSpan<double, 2>(m, "DPBuffer_Uvw");

  py::class_<PyDpBuffer>(m, "DPBuffer")
      .def(py::init<>())
      .def(py::init<double, double>(),
           "Constructor with arguments: time, exposure")
      .def(
          "get_time", [](const PyDpBuffer& self) { return self->GetTime(); },
          "Get the time of this buffer.")
      .def(
          "set_time",
          [](PyDpBuffer& self, double time) { self->SetTime(time); },
          "Set the time of this buffer.")
      .def(
          "get_exposure",
          [](const PyDpBuffer& self) { return self->GetExposure(); },
          "Get the exposure duration of this buffer.")
      .def(
          "set_exposure",
          [](PyDpBuffer& self, double exposure) {
            self->SetExposure(exposure);
          },
          "Set the exposure duration of this buffer.")
      .def(
          "get_data",
          [](PyDpBuffer& self, const std::string& name) {
            DPBuffer& buffer = *self;
            if (!buffer.HasData(name)) {
              throw std::runtime_error("Buffer has no data named '" + name +
                                       "'");
            }
            return aocommon::xt::CreateSpan(buffer.GetData(name));
          },
          "Get data buffer that can be used as numpy array. Shape is "
          "(nr baselines, nr channels, nr polarizations). If the 'name' "
          "argument is empty (default), get the main data buffer. "
          "Otherwise, get the extra data buffer with the given name.",
          py::arg("name") = "")
      .def(
          "get_flags",
          [](PyDpBuffer& self) {
            return aocommon::xt::CreateSpan(self->GetFlags());
          },
          "Get flags buffer that can be used as numpy array. Shape is "
          "(nr baselines, nr channels, nr polarizations).")
      .def(
          "get_weights",
          [](PyDpBuffer& self) {
            return aocommon::xt::CreateSpan(self->GetWeights());
          },
          "Get weights buffer that can be used as numpy array. Shape is "
          "(nr baselines, nr channels, nr polarizations).")
      .def(
          "get_uvw",
          [](PyDpBuffer& self) {
            return aocommon::xt::CreateSpan(self->GetUvw());
          },
          "Get UVW buffer that can be used as numpy array. Shape is "
          "(nr baselines, 3).")
      .def(
          "set_data",
          [](PyDpBuffer& self,
             py::array_t<std::complex<float>, py::array::c_style>& numpy_data) {
            if (numpy_data.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations); the input provided contains " +
                  std::to_string(numpy_data.ndim()) + " dimensions.");
            }
            const std::array<std::size_t, 3> shape =
                ConvertShape(numpy_data.shape());
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            DPBuffer& buffer = *self;
            buffer.GetData().resize(shape);
            std::copy_n(numpy_data.data(), buffer.GetData().size(),
                        buffer.GetData().data());
          },
          "Set buffer data from a 3-D numpy array of complex float type.")
      .def(
          "set_extra_data",
          [](PyDpBuffer& self, std::string name,
             py::array_t<std::complex<float>, py::array::c_style>& numpy_data) {
            if (numpy_data.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations); the input provided contains " +
                  std::to_string(numpy_data.ndim()) + " dimensions.");
            }
            const std::array<std::size_t, 3> shape =
                ConvertShape(numpy_data.shape());
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            self->GetData(name).resize(shape);
            std::copy_n(numpy_data.data(), self->GetData(name).size(),
                        self->GetData(name).data());
          },
          "Set buffer data with the given name from a 3-D numpy array of "
          "complex float type.")
      .def(
          "add_data",
          [](PyDpBuffer& self, std::string name) {
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            self->AddData(name);
          },
          "Add a new data field in the DPBuffer with the given name.")
      .def(
          "remove_data",
          [](PyDpBuffer& self, std::string name) {
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            self->RemoveData(name);
          },
          "Remove the data field in the DPBuffer with the given name.")
      .def(
          "set_weights",
          [](PyDpBuffer& self,
             py::array_t<float, py::array::c_style>& numpy_weights) {
            if (numpy_weights.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            const std::array<std::size_t, 3> shape =
                ConvertShape(numpy_weights.shape());
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            DPBuffer& buffer = *self;
            buffer.GetWeights().resize(shape);
            std::copy_n(numpy_weights.data(), buffer.GetWeights().size(),
                        buffer.GetWeights().data());
          },
          "Set buffer weights from a 3-D numpy array of float type.")
      .def(
          "set_flags",
          [](PyDpBuffer& self,
             py::array_t<bool, py::array::c_style>& numpy_flags) {
            if (numpy_flags.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            const std::array<std::size_t, 3> shape =
                ConvertShape(numpy_flags.shape());
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            DPBuffer& buffer = *self;
            buffer.GetFlags().resize(shape);
            std::copy_n(numpy_flags.data(), buffer.GetFlags().size(),
                        buffer.GetFlags().data());
          },
          "Set buffer flags from a 3-D numpy array of boolean type.")
      .def(
          "set_uvw",
          [](PyDpBuffer& self,
             py::array_t<double, py::array::c_style>& numpy_uvw) {
            if (numpy_uvw.ndim() != 2) {
              throw std::runtime_error(
                  "Provided array should have two dimensions (n_baselines, "
                  "3).");
            }
            if (numpy_uvw.shape(1) != 3) {
              throw std::runtime_error(
                  "Each baseline should have 3 uvw values.");
            }
            const std::array<std::size_t, 3> shape{
                static_cast<std::size_t>(numpy_uvw.shape(0)), 3u};
            // Potentially long running operation, so release the GIL
            py::gil_scoped_release release;
            DPBuffer& buffer = *self;
            buffer.GetUvw().resize(shape);
            std::copy_n(numpy_uvw.data(), buffer.GetUvw().size(),
                        buffer.GetUvw().data());
          },
          "Set buffer uvw from a 2-D numpy array of double type.");
}

}  // namespace pythondp3
}  // namespace dp3
