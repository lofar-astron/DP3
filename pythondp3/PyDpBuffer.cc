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

template <typename T, size_t dimensions>
inline std::vector<py::ssize_t> XTShapeToPyShape(
    const typename xt::xtensor<T, dimensions>::shape_type& shape) {
  std::vector<py::ssize_t> py_shape(shape.size());
  std::transform(shape.begin(), shape.end(), py_shape.begin(),
                 [](size_t dimension_size) {
                   return static_cast<py::ssize_t>(dimension_size);
                 });
  return py_shape;
}

template <typename T, size_t dimensions>
inline std::vector<py::ssize_t> XTStridesToPyStrides(
    const typename xt::xtensor<T, dimensions>::strides_type& strides) {
  std::vector<py::ssize_t> py_strides(strides.size());
  std::transform(
      strides.begin(), strides.end(), py_strides.begin(),
      [](size_t dimension_stride) {
        return static_cast<py::ssize_t>(dimension_stride * sizeof(T));
      });
  return py_strides;
}

template <typename T, size_t dimensions>
py::array_t<T> XArrayToPyArray(const xt::xtensor<T, dimensions>& source,
                               py::handle base_obj) {
  std::vector<py::ssize_t> strides =
      XTStridesToPyStrides<T, dimensions>(source.strides());
  std::vector<py::ssize_t> shape =
      XTShapeToPyShape<T, dimensions>(source.shape());

  py::array_t<T> out(shape, strides, source.data(), base_obj);

  out.attr("setflags")(false);
  return out;
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
          "Set buffer uvw from a 2-D numpy array of double type.")
      .def(
          "get_solution",
          [](PyDpBuffer& self) {
            py::dict solution_per_directions;
            for (auto& [direction_name, direction_solutions] :
                 self->GetSolution()) {
              solution_per_directions[py::str(direction_name)] = std::move(
                  XArrayToPyArray(direction_solutions, py::cast(&self)));
              ;
            }

            return solution_per_directions;
          },
          "Get gain solutions from the buffer.\n"
          "\n"
          "Returns\n"
          "-------\n"
          "solutions : dict[str, numpy.ndarray]\n"
          "    Mapping from direction name to complex-valued gain solutions.\n"
          "    Each value is a numpy.ndarray with dtype complex128 and shape\n"
          "    (n_channels, n_antennas, n_polarizations).\n",
          py::return_value_policy::reference_internal)
      .def(
          "set_solution",
          [](PyDpBuffer& self, const py::dict& solution_per_direction) {
            DPBuffer::SolutionType solutions;

            for (auto& [py_direction_name, py_solution_per_directions] :
                 solution_per_direction) {
              const std::string direction_name =
                  py::cast<std::string>(py_direction_name);

              py::array_t<std::complex<double>,
                          py::array::c_style | py::array::forcecast>
                  solution_per_direction =
                      py::cast<py::array>(py_solution_per_directions);

              constexpr ssize_t kExpectedDimensions = 3;
              if (solution_per_direction.ndim() != kExpectedDimensions) {
                throw std::runtime_error(
                    "Each solution array must have shape "
                    "(n_channels, n_antennas, n_polarizations).");
              }

              const size_t n_channels =
                  static_cast<size_t>(solution_per_direction.shape(0));
              const size_t n_antennas =
                  static_cast<size_t>(solution_per_direction.shape(1));
              const size_t n_polarizations =
                  static_cast<size_t>(solution_per_direction.shape(2));

              xt::xtensor<std::complex<double>, 3> solution(
                  {n_channels, n_antennas, n_polarizations});
              std::copy_n(solution_per_direction.data(), solution.size(),
                          solution.data());

              solutions.emplace(direction_name, std::move(solution));
            }

            self->SetSolution(solutions);
          },
          "Set gain solutions in the buffer.\n"
          "\n"
          "Parameters\n"
          "----------\n"
          "solutions : dict[str, numpy.ndarray]\n"
          "    Mapping from direction name to gain solutions.\n"
          "    Each value must be a complex-valued numpy.ndarray with dtype\n"
          "    complex128 and shape (n_channels, n_antennas, "
          "n_polarizations).\n");
}

}  // namespace pythondp3
}  // namespace dp3
