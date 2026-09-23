// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "base/Telescope.h"

#include <boost/test/unit_test.hpp>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <EveryBeam/everybeam.h>
#include <EveryBeam/telescope/phasedarray.h>

using dp3::base::GetStationIndices;
using dp3::base::GetTelescope;

BOOST_AUTO_TEST_SUITE(telescope)

BOOST_AUTO_TEST_CASE(read_lofar) {
  const std::string kMsName = "tNDPPP-generic.MS";
  const casacore::MeasurementSet kMs(kMsName);
  const std::vector<std::string> kAntennaNames = {"CS001HBA0", "CS002HBA0"};

  everybeam::Options everybeam_options;
  everybeam_options.element_response_model =
      everybeam::ElementResponseModel::kHamaker;
  everybeam_options.use_channel_frequency = false;

  std::unique_ptr<everybeam::Telescope> telescope =
      everybeam::LoadTelescope(kMs, everybeam_options);
  const std::vector<size_t> station_indices =
      GetStationIndices(*telescope, kAntennaNames);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(telescope->GetStationName(station_indices[i]),
                      kAntennaNames[i]);
  }
}

BOOST_AUTO_TEST_CASE(read_oskar) {
  const std::string kMsName = "tOSKAR.in_MS";
  const casacore::MeasurementSet kMs(kMsName);
  const std::vector<std::string> kAntennaNames = {"s0012", "s0013", "s0015"};

  everybeam::Options everybeam_options;
  everybeam_options.element_response_model =
      everybeam::ElementResponseModel::kOSKARSphericalWave;
  everybeam_options.use_channel_frequency = true;

  std::unique_ptr<everybeam::Telescope> telescope =
      everybeam::LoadTelescope(kMs, everybeam_options);

  const std::vector<size_t> station_indices =
      GetStationIndices(*telescope, kAntennaNames);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(telescope->GetStationName(station_indices[i]),
                      kAntennaNames[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
