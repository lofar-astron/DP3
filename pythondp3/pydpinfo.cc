// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <casacore/casa/Quanta/Quantum.h>

#include <dp3/base/DPInfo.h>

using dp3::base::DPInfo;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapDpInfo(py::module& m) {
  py::class_<DPInfo>(m, "DPInfo")
      .def(py::init<>())

      /* Antenna properties */
      .def(
          "set_antennas",
          [](DPInfo& self, const std::vector<std::string>& names,
             const std::vector<double>& diameters,
             const std::vector<std::array<double, 3>>& positions,
             const std::vector<int>& antenna1_indices,
             const std::vector<int>& antenna2_indices) {
            std::vector<casacore::MPosition> positions_casa;
            positions_casa.reserve(positions.size());
            for (const std::array<double, 3>& position : positions) {
              casacore::Vector<double> vector(position.size());
              std::copy_n(position.data(), position.size(), vector.data());
              casacore::Quantum<casacore::Vector<double>> quantum(vector, "m");
              positions_casa.emplace_back(quantum, casacore::MPosition::ITRF);
            }
            self.setAntennas(names, diameters, positions_casa, antenna1_indices,
                             antenna2_indices);
          },
          R"(Set antenna properties.
Parameters:
  names: List with the names for each antenna.
  diameters: List with the diameters for each antenna, in meters.
  positions: List with the antenna positions in ITRF XYZ format.
             An antenna position is a list with x, y, z values in meters.
  first_antenna_indices: List with the index of the first antenna for each
                         baseline.
  second_antenna_indices: List with the index of the second antenna for
                          each baseline.)",
          py::arg("names"), py::arg("diameters"), py::arg("positions"),
          py::arg("first_antenna_indices"), py::arg("second_antenna_indices"))
      .def_property_readonly("antenna_names", &DPInfo::antennaNames,
                             "A list with antenna names (read only)")
      .def_property_readonly(
          "antenna_positions",
          [](DPInfo& self) {
            // Convert vector of casa MPositions to std::vector of positions
            std::vector<casacore::MPosition> positions_casa = self.antennaPos();
            std::vector<std::array<double, 3>> positions;
            for (size_t i = 0; i < positions_casa.size(); ++i) {
              casacore::MVPosition position_mv = positions_casa[i].getValue();
              std::array<double, 3> position_array = {
                  position_mv(0), position_mv(1), position_mv(2)};
              positions.push_back(position_array);
            }
            return positions;
          },
          "A list of antenna positions in ITRF XYZ format (read only)")
      .def_property_readonly("first_antenna_indices", &DPInfo::getAnt1,
                             "A list with index of the first antenna for each "
                             "baseline (read-only)")
      .def_property_readonly("second_antenna_indices", &DPInfo::getAnt2,
                             "A list with index of the second antenna for each "
                             "baseline (read-only)")
      .def_property_readonly("n_antenna", &DPInfo::nantenna,
                             "The number of antennas (read-only)")

      /* Other properties */
      // TODO(AST-1097): Add more bindings and use proper property bindings. */
      .def("get_channel_frequencies", &DPInfo::chanFreqs,
           py::arg("baseline") = 0,
           "Get a list of channel frequencies (read only)")
      .def("nchan", &DPInfo::nchan, "Get the number of channels")
      .def("start_time", &DPInfo::startTime, "Get the start time")
      .def("time_interval", &DPInfo::timeInterval, "Get the time interval")
      .def("ntime", &DPInfo::ntime, "Get the total number of time slots")
      .def("ms_name", &DPInfo::msName, "Get name of measurement set");
}

}  // namespace pythondp3
}  // namespace dp3
