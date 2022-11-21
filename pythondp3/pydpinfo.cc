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
      .def(py::init<unsigned int, unsigned int, unsigned int, std::string>(),
           py::arg("n_correlations") = 0, py::arg("n_original_channels") = 0,
           py::arg("start_channel") = 0, py::arg("antenna_set") = "")
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

      /* Channel properties. Not all properties are bound in python. */
      .def(
          "set_channels",
          [](DPInfo& self, std::vector<double> frequencies,
             std::vector<double> widths) {
            self.setChannels(std::move(frequencies), std::move(widths));
          },
          R"(Set (basic) channel properties".
Parameters:
  frequencies: List with the frequencies for each channel.
  widths: List with the widths of each channel.
          Typically, all channels have the same width.)",
          py::arg("frequencies"), py::arg("widths"))
      .def_property_readonly("n_channels", &DPInfo::nchan, "Number of channels")
      .def_property_readonly("original_n_channels", &DPInfo::origNChan,
                             "Original number of channels")
      .def_property_readonly(
          "start_channel", &DPInfo::startchan,
          "Number of the first channel. Channel numbers start at 0.")
      .def_property_readonly("n_correlations", &DPInfo::ncorr,
                             "Number of correlations")
      .def_property_readonly("antenna_set", &DPInfo::antennaSet,
                             "Antenna set (LOFAR specific)")
      .def_property_readonly(
          "channel_frequencies", [](DPInfo& self) { return self.chanFreqs(0); },
          R"(A list of channel frequencies (read only).
When using BDA, this property has the values for the first baseline.

See Also
--------
bda_channel_frequencies
)")
      .def_property_readonly(
          "bda_channel_frequencies", &DPInfo::BdaChanFreqs,
          R"(A list with, for each baseline, a list of channel frequencies (read only).
When not using BDA, there is always a single frequency list.)")
      .def_property_readonly(
          "channel_widths", [](DPInfo& self) { return self.chanWidths(0); },
          R"(A list of channel widths (read only).
When using BDA, this property has the values for the first baseline.

See Also
--------
bda_channel_widths)")
      .def_property_readonly(
          "bda_channel_widths", &DPInfo::BdaChanWidths,
          R"(A list with, for each baseline, a list of channel widths (read only).
When not using BDA, there is always a single list of channel widths.)")

      /* Time properties */
      .def("set_times", &DPInfo::setTimes,
           R"(Set basic time properties.
Parameters:
  first_time: Centroid time of the first time slot (mjd seconds).
  last_time: Centroid time of the last time slot (mjd seconds).
  time_interval: Time interval between two time slots (seconds)",
           py::arg("first_time"), py::arg("last_time"),
           py::arg("time_interval"))
      .def_property_readonly(
          "start_time", &DPInfo::startTime,
          R"(The time point at which the measurements started (mjd seconds).
The start time equals the first time minus half a time interval)")
      .def_property_readonly(
          "first_time", &DPInfo::firstTime,
          "The centroid time of the first time slot (mjd seconds)")
      .def_property_readonly(
          "last_time", &DPInfo::lastTime,
          "The centroid time of the last time slot (mjd seconds)")
      .def_property_readonly(
          "time_interval", &DPInfo::timeInterval,
          "The time interval between two time slots (seconds)")
      .def_property_readonly("n_times", &DPInfo::ntime,
                             "The number of time slots")

      /* Array information */
      .def_property(
          "phase_center",
          [](DPInfo& self) {
            return std::array<double, 2>{
                self.phaseCenter().getValue().get()[0],
                self.phaseCenter().getValue().get()[1]};
          },
          [](DPInfo& self, const std::array<double, 2>& ra_dec) {
            const casacore::Quantity ra(ra_dec[0], "rad");
            const casacore::Quantity dec(ra_dec[1], "rad");
            const casacore::MDirection direction(ra, dec,
                                                 casacore::MDirection::J2000);
            self.setPhaseCenter(direction);
          },
          "The phase center for the measurements ([ra, dec] in radians)")

      /* Other properties */
      .def_property("ms_name", &DPInfo::msName, &DPInfo::setMsName,
                    "The name of the measurement set");
}

}  // namespace pythondp3
}  // namespace dp3
