// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <dp3/base/DPBuffer.h>

using dp3::base::DPBuffer;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

namespace {

template <typename TensorType>
void RegisterTensor(py::module& m, const char* name) {
  py::class_<TensorType>(m, name, py::buffer_protocol())
      .def_buffer([](TensorType& tensor) {
        using ValueType = typename TensorType::value_type;
        constexpr std::size_t value_size = sizeof(ValueType);
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
  RegisterTensor<DPBuffer::DataType>(m, "DPBuffer_Data");
  RegisterTensor<DPBuffer::WeightsType>(m, "DPBuffer_Weights");
  RegisterTensor<DPBuffer::FlagsType>(m, "DPBuffer_Flags");
  RegisterTensor<DPBuffer::UvwType>(m, "DPBuffer_Uvw");

  py::class_<DPBuffer, std::shared_ptr<DPBuffer>>(m, "DPBuffer")
      .def(py::init<>())
      .def(py::init<double, double>(),
           "Constructor with arguments: time, exposure")
      .def("get_time", &DPBuffer::getTime, "Get the time of this buffer.")
      .def("set_time", &DPBuffer::setTime, "Set the time of this buffer.")
      .def("get_exposure", &DPBuffer::getExposure,
           "Get the exposure duration of this buffer.")
      .def("set_exposure", &DPBuffer::setExposure,
           "Set the exposure duration of this buffer.")
      .def(
          "get_data",
          [](DPBuffer& self, const std::string& name) {
            if (!self.HasData(name)) {
              throw std::runtime_error("Buffer has no data named '" + name +
                                       "'");
            }
            return self.GetData(name);
          },
          "Get data buffer that can be used as numpy array. Shape is "
          "(nr baselines, nr channels, nr polarizations). If the 'name' "
          "argument is empty (default), get the main data buffer. "
          "Otherwise, get the extra data buffer with the given name.",
          py::arg("name") = "")
      .def("get_weights",
           (DPBuffer::WeightsType & (DPBuffer::*)())(&DPBuffer::GetWeights),
           "Get weights buffer that can be used as numpy array. Shape is "
           "(nr baselines, nr channels, nr polarizations).")
      .def("get_flags",
           (DPBuffer::FlagsType & (DPBuffer::*)())(&DPBuffer::GetFlags),
           "Get flags buffer that can be used as numpy array. Shape is "
           "(nr baselines, nr channels, nr polarizations).")
      .def("get_uvw", (DPBuffer::UvwType & (DPBuffer::*)())(&DPBuffer::GetUvw),
           "Get UVW buffer that can be used as numpy array. Shape is "
           "(nr baselines, 3).")
      .def(
          "set_data",
          [](DPBuffer& self,
             py::array_t<std::complex<float>, py::array::c_style>& numpy_data) {
            if (numpy_data.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            self.ResizeData(ConvertShape(numpy_data.shape()));
            std::copy_n(numpy_data.data(), self.GetData().size(),
                        self.GetData().data());
          },
          "Set buffer data from a 3-D numpy array of complex float type.")
      .def(
          "set_weights",
          [](DPBuffer& self,
             py::array_t<float, py::array::c_style>& numpy_weights) {
            if (numpy_weights.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            self.ResizeWeights(ConvertShape(numpy_weights.shape()));
            std::copy_n(numpy_weights.data(), self.GetWeights().size(),
                        self.GetWeights().data());
          },
          "Set buffer weights from a 3-D numpy array of float type.")
      .def(
          "set_flags",
          [](DPBuffer& self,
             py::array_t<bool, py::array::c_style>& numpy_flags) {
            if (numpy_flags.ndim() != 3) {
              throw std::runtime_error(
                  "Provided array should have three dimensions (n_baselines, "
                  "n_channels, n_correlations).");
            }
            self.ResizeFlags(ConvertShape(numpy_flags.shape()));
            std::copy_n(numpy_flags.data(), self.GetFlags().size(),
                        self.GetFlags().data());
          },
          "Set buffer flags from a 3-D numpy array of boolean type.")
      .def(
          "set_uvw",
          [](DPBuffer& self,
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
            self.ResizeUvw(numpy_uvw.shape(0));
            std::copy_n(numpy_uvw.data(), self.GetUvw().size(),
                        self.GetUvw().data());
          },
          "Set buffer uvw from a 2-D numpy array of double type.");
}

}  // namespace pythondp3
}  // namespace dp3
