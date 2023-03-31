// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "../ddecal/constraints/TECConstraint.h"

using dp3::ddecal::Constraint;
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

  py::class_<Constraint::Result>(m, "Result",
                                 "Constraint result for a single type. "
                                 "The values are in result.values.")
      .def_readwrite("values", &Constraint::Result::vals,
                     "Result values as a flat list. .shape defines the shape.")
      .def_readwrite("weights", &Constraint::Result::weights,
                     "Weights for the values as a flat list. "
                     ".shape defines the shape.")
      .def_readwrite("axes", &Constraint::Result::axes,
                     "Comma-separated dimension names.")
      .def_readwrite("shape", &Constraint::Result::dims,
                     "Shape of .values and .weights, as a list with the size "
                     "of each dimension.")
      .def_readwrite("name", &Constraint::Result::name,
                     "Identifier for the result type.")
      .def("__str__",
           [](const Constraint::Result& result) {
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
      .def("__repr__", [](const Constraint::Result& result) {
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

        std::vector<Constraint::Result> results =
            constraint.Apply(gains_span, time, nullptr);

        // Remove unused axes from the result.
        for (Constraint::Result& result : results) {
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
  spectrum: 1-D or 2-D numpy array with complex values:
            (frequencies) or (frequencies, antennas)
            If there is one dimension, the number of antennas is always one.
            The number of frequencies (first dimension) must match the length
            of the 'frequencies' argument.
            The TEC fitter adjusts the values in this parameter.
  frequencies: Channel frequencies in Hz.
  time: Optional. Time in seconds. The TEC fitter does not use this parameter.
  mode: Optional. Fitting mode: "tec_and_common_scalar" (default) or "tec_only".

Returns: Result with the values for tec or tec and common scalar phase, along
         with an error.)");
}

}  // namespace pythondp3
}  // namespace dp3
