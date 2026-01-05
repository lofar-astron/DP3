// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "ddecal/constraints/SmoothnessConstraint.h"
#include "ddecal/constraints/TECConstraint.h"

using dp3::ddecal::Constraint;
using dp3::ddecal::ConstraintResult;
using dp3::ddecal::SmoothnessConstraint;
using dp3::ddecal::TECConstraint;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

namespace {
// pybind11 cannot convert string arguments to enums.
TECConstraint::Mode TecModeFromString(const std::string& mode) {
  if (mode == "tec_and_common_scalar") {
    return TECConstraint::Mode::kTecAndCommonScalar;
  } else if (mode == "tec_only") {
    return TECConstraint::Mode::kTecOnly;
  } else {
    throw std::invalid_argument("Invalid TEC mode: " + mode);
  }
}
}  // namespace

PYBIND11_MODULE(fitters, m) {
  m.doc() = "DP3 spectral fitters";

  // Expose the class of the return variable of the fitter functions to Python.
  py::class_<ConstraintResult>(m, "Result",
                               "Constraint result for a single type. "
                               "The values are in result.values.")
      .def_readwrite("values", &ConstraintResult::vals,
                     "Result values as a flat list. .shape defines the shape.")
      .def_readwrite("weights", &ConstraintResult::weights,
                     "Weights for the values as a flat list. "
                     ".shape defines the shape.")
      .def_readwrite("axes", &ConstraintResult::axes,
                     "Comma-separated dimension names.")
      .def_readwrite("shape", &ConstraintResult::dims,
                     "Shape of .values and .weights, as a list with the size "
                     "of each dimension.")
      .def_readwrite("name", &ConstraintResult::name,
                     "Identifier for the result type.")
      .def("__str__",
           [](const ConstraintResult& result) {
             std::stringstream stream;
             stream << "[Result name:" << result.name + " axes:" << result.axes
                    << " shape:";
             if (result.dims.empty()) {
               stream << "(none)";
             } else {
               stream << result.dims.front();
               for (std::size_t i = 1; i < result.dims.size(); ++i) {
                 stream << "," << result.dims[i];
               }
             }
             stream << "]";
             return stream.str();
           })
      .def("__repr__", [](const ConstraintResult& result) {
        std::stringstream stream;
        stream << "<dp3.fitters.Result name=" << result.name + " axes="
               << result.axes << " values=[";
        for (const double& value : result.vals) {
          stream << value << ",";
        }
        stream << "] weights=[";
        for (const double& weight : result.weights) {
          stream << weight << ",";
        }
        stream << "]>";
        return stream.str();
      });

  // Create bindings for the internal fitter functions, which allows using these
  // functions directly from Python.
  m.def(
      "fit_tec",
      [](py::array_t<std::complex<double>, py::array::c_style>& gains,
         const std::vector<double>& frequencies, double time,
         const std::string& mode) {
        if ((gains.ndim() != 1 && gains.ndim() != 2) ||
            static_cast<size_t>(gains.shape(0)) != frequencies.size()) {
          throw std::invalid_argument("Incorrect gains shape");
        }

        const std::size_t n_antennas = (gains.ndim() == 1) ? 1 : gains.shape(1);
        const std::size_t kNSolutions = 1;
        const std::vector<uint32_t> kNSolutionsPerDirection = {1};

        TECConstraint constraint(TecModeFromString(mode));
        constraint.Initialize(n_antennas, kNSolutionsPerDirection, frequencies);
        // With a single antenna, the result is always zero if phase referencing
        // is enabled: That single antenna is then used as reference.
        if (n_antennas == 1) constraint.setDoPhaseReference(false);

        const std::array<std::size_t, 4> shape{frequencies.size(), n_antennas,
                                               kNSolutions, 1};
        aocommon::xt::Span<std::complex<double>, 4> gains_span =
            aocommon::xt::CreateSpan(gains.mutable_data(), shape);

        std::vector<ConstraintResult> results =
            constraint.Apply(gains_span, time, nullptr);

        // Remove unused axes from the result.
        for (ConstraintResult& result : results) {
          assert(result.axes == "ant,dir,freq");
          assert((result.dims == std::vector<size_t>{n_antennas, 1, 1}));
          result.axes = "ant";
          result.dims = {n_antennas};
        }

        return results;
      },
      py::arg("gains"), py::arg("frequencies"), py::arg("time") = 0.0,
      py::arg("mode") = "tec_and_common_scalar",
      R"(Apply a TEC fitter to complex gains.

The fitter fits a function through the phases of the complex numbers of each
antenna, and overwrites the gains by their constrained values. The function to
which the phases are fitted is -8.44797245e9 * alpha / freq if 'mode' is set to
'tec_only', or -8.44797245e9 * alpha / freq + beta if 'mode' is set to
'tec_and_common_scalar'. The values for alpha (tec_only) or alpha and beta
(tec_and_common_scalar) are stored in the returned data structure.

When using multiple antennas, the fitter uses the phase of the first antenna
as reference. The alpha, beta and error value for the first antenna are then
zero. The phases for the other antennas are relative to the first antenna.
With a single antenna, this phase referencing is disabled.

Parameters:
  gains: 1-D or 2-D numpy array with complex values:
        (frequencies) or (frequencies, antennas)
        If there is one dimension, the number of antennas is always one.
        The number of frequencies (first dimension) must match the length
        of the 'frequencies' argument.
        The TEC fitter adjusts the values in this parameter.
  frequencies: Channel frequencies in Hz.
  time: Optional. Time in seconds. The TEC fitter does not use this parameter.
  mode: Optional. Fitting mode: "tec_and_common_scalar" (default) or "tec_only".

Returns: A list of two or three dp3.fitters.Result items. When the mode is
        "tec_and_common_scalar", it contains: [TEC, common scalar phase, error],
        and otherwise [TEC, error].)");

  m.def(
      "fit_smooth",
      [](py::array_t<std::complex<double>, py::array::c_style>& gains,
         const std::vector<double>& frequencies, double time,
         const double bandwidth_hz, const double bandwidth_ref_frequency_hz,
         double spectral_exponent, bool kernel_truncation) {
        if ((gains.ndim() != 1 && gains.ndim() != 2) ||
            static_cast<size_t>(gains.shape(0)) != frequencies.size()) {
          throw std::invalid_argument("Incorrect gains shape");
        }
        if (bandwidth_hz == 0.0) {
          throw std::invalid_argument("bandwidth_hz may not be zero");
        }

        const std::size_t n_channel_blocks = frequencies.size();
        const std::size_t n_antennas = (gains.ndim() == 1) ? 1 : gains.shape(1);
        const std::size_t kNDirections = 1;
        const std::size_t kNPolarizations = 1;
        const std::vector<uint32_t> kNSolutionsPerDirection(kNDirections, 1);

        const std::vector<double> weights(n_channel_blocks, 1.0);
        std::vector<double> antenna_distance_factors{1.0};

        SmoothnessConstraint constraint(bandwidth_hz,
                                        bandwidth_ref_frequency_hz,
                                        spectral_exponent, kernel_truncation);
        constraint.Initialize(n_antennas, kNSolutionsPerDirection, frequencies);
        constraint.SetAntennaFactors(std::move(antenna_distance_factors));
        constraint.SetWeights(weights);

        const std::array<std::size_t, 4> shape{n_channel_blocks, n_antennas,
                                               kNDirections, kNPolarizations};
        aocommon::xt::Span<std::complex<double>, 4> gains_span =
            aocommon::xt::CreateSpan(gains.mutable_data(), shape);

        constraint.Apply(gains_span, time, nullptr);
      },
      py::arg("gains"), py::arg("frequencies"), py::arg("time") = 0.0,
      py::arg("bandwidth_hz"), py::arg("bandwidth_ref_frequency_hz") = 0.0,
      py::arg("spectral_exponent") = -1.0, py::arg("kernel_truncation") = true,
      R"(Apply a Smoothness fitter to complex gains.

The fitter smooths a series of possibly irregularly gridded values by a given
Gaussian kernel. The Gaussian kernel is trimmed off at 3 sigma, and further
defined by two bandwidth parameters.

Parameters:
  gains: 2-D (or 1-D) numpy array with complex values of shape
          (frequencies, antennas) or (frequencies).
        If the number of antennas is one, it can be provided as a 1-D array.
        The number of frequencies (first dimension) must match the length
        of the 'frequencies' argument.
        The Smoothness fitter adjusts the values in this parameter.
  frequencies: Array defining the frequency of each channel, in Hz.
  time: Optional. Time in seconds. The fitter does not use this parameter.
  bandwidth_hz: Size, in frequency units (Hz), of the Gaussian kernel
                (smoothing strength) that is used for smoothing. May not be zero.
  bandwidthRefFrequencyHz: Optional. Kernel size over frequency. May be zero
                          (default) to have a constant kernel size over frequency.

Returns: None. The actual return is written in-place to the gains argument.)");
}

}  // namespace pythondp3
}  // namespace dp3
