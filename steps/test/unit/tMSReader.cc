// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../MSReader.h"
#include "../../../common/ParameterSet.h"

#include <EveryBeam/load.h>
#include <EveryBeam/telescope/phasedarray.h>

#include <casacore/casa/Arrays/Vector.h>
#include <vector>

using dp3::steps::MSReader;

BOOST_AUTO_TEST_SUITE(msreader)

// Test reading a LOFAR measurement set
BOOST_AUTO_TEST_CASE(read_lofar) {
  const casacore::MeasurementSet ms("tNDPPP-generic.MS");
  dp3::common::ParameterSet parset;
  MSReader msreader(ms, parset, "");

  std::vector<std::string> ant_vec = {"CS001HBA0", "CS002HBA0"};
  casacore::Vector<casacore::String> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      msreader.GetTelescope(everybeam::ElementResponseModel::kHamaker, false);
  const std::vector<size_t> station_indices =
      MSReader::SelectStationIndices(telescope.get(), ant_names);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(phasedarray.GetStation(station_indices[i])->GetName(),
                      ant_vec[i]);
  }
}

// Reading an OSKAR measurement set
BOOST_AUTO_TEST_CASE(read_oskar) {
  const casacore::MeasurementSet ms("tOSKAR.in_MS");
  dp3::common::ParameterSet parset;
  MSReader msreader(ms, parset, "");

  std::vector<std::string> ant_vec = {"s0012", "s0013", "s0015"};
  casacore::Vector<casacore::String> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      msreader.GetTelescope(
          everybeam::ElementResponseModel::kOSKARSphericalWave, false);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);
  std::vector<size_t> station_indices =
      MSReader::SelectStationIndices(telescope.get(), ant_names);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(phasedarray.GetStation(station_indices[i])->GetName(),
                      ant_vec[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
