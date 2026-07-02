// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/beammode.h>
#include <EveryBeam/elementresponse.h>
#include <EveryBeam/load.h>
#include <EveryBeam/options.h>
#include <EveryBeam/stationnode.h>
#include <EveryBeam/telescope/telescope.h>

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapEveryBeam(py::module& m) {
  py::module everybeam_module =
      m.def_submodule("everybeam", "EveryBeam bindings bundled with DP3.");

  py::enum_<everybeam::BeamMode>(everybeam_module, "BeamMode",
                                 py::module_local())
      .value("none", everybeam::BeamMode::kNone)
      .value("full", everybeam::BeamMode::kFull)
      .value("array_factor", everybeam::BeamMode::kArrayFactor)
      .value("element", everybeam::BeamMode::kElement);

  py::enum_<everybeam::BeamNormalisationMode>(
      everybeam_module, "BeamNormalisationMode", py::module_local())
      .value("none", everybeam::BeamNormalisationMode::kNone)
      .value("pre_applied", everybeam::BeamNormalisationMode::kPreApplied)
      .value("pre_applied_or_full",
             everybeam::BeamNormalisationMode::kPreAppliedOrFull)
      .value("amplitude", everybeam::BeamNormalisationMode::kAmplitude)
      .value("full", everybeam::BeamNormalisationMode::kFull);

  py::enum_<everybeam::ElementResponseModel>(
      everybeam_module, "ElementResponseModel", py::module_local())
      .value("default_", everybeam::ElementResponseModel::kDefault)
      .value("hamaker", everybeam::ElementResponseModel::kHamaker)
      .value("hamaker_lba", everybeam::ElementResponseModel::kHamakerLba)
      .value("lobes", everybeam::ElementResponseModel::kLOBES)
      .value("oskar_dipole", everybeam::ElementResponseModel::kOSKARDipole)
      .value("skala40_wave",
             everybeam::ElementResponseModel::kOSKARSphericalWave)
      .value("skala40_spherical",
             everybeam::ElementResponseModel::kOSKARSphericalWave)
      .value("ska_mid_analytical",
             everybeam::ElementResponseModel::kSkaMidAnalytical)
      .value("aartfaac_inner", everybeam::ElementResponseModel::kAartfaacInner)
      .value("aartfaac_outer", everybeam::ElementResponseModel::kAartfaacOuter)
      .value("lwa", everybeam::ElementResponseModel::kLwa)
      .value("oskar_dipole_cos",
             everybeam::ElementResponseModel::kOSKARDipoleCos)
      .value("ska_low_feko", everybeam::ElementResponseModel::kSkaLowFeko);

  py::enum_<everybeam::TelescopeType>(everybeam_module, "TelescopeType",
                                      py::module_local())
      .value("AARTFAAC", everybeam::TelescopeType::kAARTFAAC)
      .value("LOFAR", everybeam::TelescopeType::kLofarTelescope)
      .value("OSKAR", everybeam::TelescopeType::kOSKARTelescope)
      .value("SKA_MID", everybeam::TelescopeType::kSkaMidTelescope);

  py::class_<everybeam::Options>(everybeam_module, "Options",
                                 py::module_local())
      .def(py::init<>())
      .def_readwrite("coeff_path", &everybeam::Options::coeff_path)
      .def_readwrite("beam_normalisation_mode",
                     &everybeam::Options::beam_normalisation_mode)
      .def_readwrite("use_channel_frequency",
                     &everybeam::Options::use_channel_frequency)
      .def_readwrite("data_column_name", &everybeam::Options::data_column_name)
      .def_readwrite("element_response_model",
                     &everybeam::Options::element_response_model)
      .def_readwrite("beam_mode", &everybeam::Options::beam_mode)
      .def_readwrite("frequency_interpolation",
                     &everybeam::Options::frequency_interpolation);

  using Axes = everybeam::StationCoordinateSystem::Axes;
  py::class_<Axes>(everybeam_module, "StationCoordinateSystemAxes",
                   py::module_local())
      .def(py::init([](const std::array<double, 3>& p,
                       const std::array<double, 3>& q,
                       const std::array<double, 3>& r) {
             return Axes{p, q, r};
           }),
           py::arg("p") = std::array<double, 3>{1.0, 0.0, 0.0},
           py::arg("q") = std::array<double, 3>{0.0, 1.0, 0.0},
           py::arg("r") = std::array<double, 3>{0.0, 0.0, 1.0})
      .def_readwrite("p", &Axes::p)
      .def_readwrite("q", &Axes::q)
      .def_readwrite("r", &Axes::r);

  py::class_<everybeam::StationCoordinateSystem>(
      everybeam_module, "StationCoordinateSystem", py::module_local())
      .def(py::init([](const std::array<double, 3>& origin, const Axes& axes) {
             everybeam::StationCoordinateSystem coordinate_system;
             coordinate_system.origin = origin;
             coordinate_system.axes = axes;
             return coordinate_system;
           }),
           py::arg("origin") = std::array<double, 3>{0.0, 0.0, 0.0},
           py::arg("axes") =
               Axes{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}})
      .def_readwrite("origin", &everybeam::StationCoordinateSystem::origin)
      .def_readwrite("axes", &everybeam::StationCoordinateSystem::axes);

  py::class_<everybeam::StationNode>(everybeam_module, "StationNode",
                                     py::module_local())
      .def(py::init<const everybeam::StationCoordinateSystem&,
                    const std::string&>(),
           py::arg("coordinate_system") = everybeam::StationCoordinateSystem(),
           py::arg("name") = "")
      .def(
          "add_child_node",
          [](everybeam::StationNode& self, everybeam::StationNode child,
             const std::array<double, 3>& position) {
            self.AddChild(std::move(child), position);
          },
          py::arg("child"), py::arg("position"))
      .def("add_child_element",
           py::overload_cast<std::array<double, 3>, bool, bool>(
               &everybeam::StationNode::AddChild),
           py::arg("position"), py::arg("is_x_enabled") = true,
           py::arg("is_y_enabled") = true);

  py::class_<everybeam::telescope::Telescope,
             std::shared_ptr<everybeam::telescope::Telescope>>(
      everybeam_module, "Telescope", py::module_local());

  everybeam_module.def(
      "create_telescope",
      [](everybeam::TelescopeType telescope_type,
         const everybeam::Options& options,
         const everybeam::StationNode& station_tree,
         const std::vector<std::array<double, 2>>& delay_directions,
         const std::array<double, 2>& tile_beam_direction,
         const std::array<double, 2>& preapplied_beam_direction,
         everybeam::BeamMode preapplied_beam_mode,
         const std::vector<double>& dish_diameters, double reference_frequency,
         const std::vector<int>& mwa_delay_factors)
          -> std::shared_ptr<everybeam::telescope::Telescope> {
        return std::shared_ptr<everybeam::telescope::Telescope>(
            everybeam::CreateTelescope(telescope_type, options, station_tree,
                                       delay_directions, tile_beam_direction,
                                       preapplied_beam_direction,
                                       preapplied_beam_mode, dish_diameters,
                                       reference_frequency, mwa_delay_factors));
      },
      "Create an EveryBeam Telescope object from specified metadata.",
      py::arg("telescope_type"), py::arg("options"), py::arg("station_tree"),
      py::arg("delay_directions") = std::vector<std::array<double, 2>>(),
      py::arg("tile_beam_direction") = std::array<double, 2>{0.0, 0.0},
      py::arg("preapplied_beam_direction") = std::array<double, 2>{0.0, 0.0},
      py::arg("preapplied_beam_mode") = everybeam::BeamMode::kNone,
      py::arg("dish_diameters") = std::vector<double>(),
      py::arg("reference_frequency") = 0.0,
      py::arg("mwa_delay_factors") = std::vector<int>());
}

}  // namespace pythondp3
}  // namespace dp3
