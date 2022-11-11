// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <dp3/base/DPInfo.h>

using dp3::base::DPInfo;

namespace py = pybind11;

namespace dp3 {
namespace pythondp3 {

void WrapDpInfo(py::module &m) {
  py::class_<DPInfo>(m, "DPInfo")
      .def(
          "antenna_names",
          [](DPInfo &self) -> py::array {
            // Convert casa vector of casa strings to std::vector of strings
            std::vector<std::string> names_casa = self.antennaNames();
            std::vector<std::string> names;
            for (size_t i = 0; i < names_casa.size(); ++i) {
              names.push_back(names_casa[i]);
            }
            py::array ret = py::cast(names);
            return ret;
          },
          "Get numpy array with antenna names (read only)")
      .def(
          "antenna_positions",
          [](DPInfo &self) -> py::array {
            // Convert vector of casa MPositions to std::vector of positions
            std::vector<casacore::MPosition> positions_casa = self.antennaPos();
            std::vector<std::array<double, 3>> positions;
            for (size_t i = 0; i < positions_casa.size(); ++i) {
              casacore::MVPosition position_mv = positions_casa[i].getValue();
              std::array<double, 3> position_array = {
                  position_mv(0), position_mv(1), position_mv(2)};
              positions.push_back(position_array);
            }
            py::array ret = py::cast(positions);
            return ret;
          },
          "Get a list of antenna positions in ITRF XYZ (read only)")
      .def("get_channel_frequencies", &DPInfo::chanFreqs,
           py::arg("baseline") = 0,
           "Get a list of channel frequencies (read only)")
      .def("get_antenna1", &DPInfo::getAnt1,
           "Get a list of first antenna numbers of correlations")
      .def("get_antenna2", &DPInfo::getAnt2,
           "Get a list of second antenna numbers of correlations")
      .def("nantenna", &DPInfo::nantenna, "Get the number of antennas")
      .def("nchan", &DPInfo::nchan, "Get the number of channels")
      .def("start_time", &DPInfo::startTime, "Get the start time")
      .def("time_interval", &DPInfo::timeInterval, "Get the time interval")
      .def("ntime", &DPInfo::ntime, "Get the total number of time slots")
      .def("ms_name", &DPInfo::msName, "Get name of measurement set");
}

}  // namespace pythondp3
}  // namespace dp3
