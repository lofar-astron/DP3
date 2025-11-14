// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Telescope.h"

#include <boost/test/unit_test.hpp>

#include <EveryBeam/telescope/phasedarray.h>

using dp3::base::GetTelescope;
using dp3::base::SelectStationIndices;

// Support both old and new return value types of PhasedArray::GetStation().
// TODO(AST-1001): Remove support for the old type in 2023.
namespace {
[[maybe_unused]] const std::string& GetStationName(
    const std::shared_ptr<const everybeam::Station>& station) {
  return station->GetName();
}

[[maybe_unused]] const std::string& GetStationName(
    const everybeam::Station& station) {
  return station.GetName();
}
}  // namespace

BOOST_AUTO_TEST_SUITE(telescope)

BOOST_AUTO_TEST_CASE(read_lofar) {
  const std::string kMsName = "tNDPPP-generic.MS";
  const std::vector<std::string> kAntennaNames = {"CS001HBA0", "CS002HBA0"};

  std::unique_ptr<everybeam::telescope::Telescope> telescope = GetTelescope(
      kMsName, everybeam::ElementResponseModel::kHamaker, false, "");
  const std::vector<size_t> station_indices =
      SelectStationIndices(*telescope, kAntennaNames);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(
        GetStationName(phasedarray.GetStation(station_indices[i])),
        kAntennaNames[i]);
  }
}

BOOST_AUTO_TEST_CASE(read_oskar) {
  const std::string kMsName = "tOSKAR.in_MS";
  const std::vector<std::string> kAntennaNames = {"s0012", "s0013", "s0015"};

  std::unique_ptr<everybeam::telescope::Telescope> telescope = GetTelescope(
      kMsName, everybeam::ElementResponseModel::kOSKARSphericalWave, true, "");
  const std::vector<size_t> station_indices =
      SelectStationIndices(*telescope, kAntennaNames);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(
        GetStationName(phasedarray.GetStation(station_indices[i])),
        kAntennaNames[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
