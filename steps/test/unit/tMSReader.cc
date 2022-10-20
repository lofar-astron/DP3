// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../MSReader.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/Telescope.h"

#include <EveryBeam/load.h>
#include <EveryBeam/telescope/phasedarray.h>

#include <casacore/casa/Arrays/Vector.h>
#include <vector>

using dp3::steps::MSReader;

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

BOOST_AUTO_TEST_SUITE(msreader)

BOOST_AUTO_TEST_CASE(provided_fields) {
  const dp3::common::Fields kFieldsToRead(
      dp3::common::Fields::Single::kFullResFlags);
  const dp3::common::Fields kWeightsField(
      dp3::common::Fields::Single::kWeights);
  // TODO(AST-1061): Use an MS that can be autoweighted, and re-add
  // 'autoweight' to kAutoweightSettings.
  const casacore::MeasurementSet ms("tNDPPP_tmp.MS");
  const dp3::common::ParameterSet parset;
  const std::vector<std::string> kAutoweightSettings = {/*"autoweight",*/
                                                        "forceautoweight"};

  MSReader msreader(ms, parset, "");
  BOOST_TEST(msreader.getProvidedFields() == dp3::common::Fields());
  msreader.setFieldsToRead(kFieldsToRead);
  BOOST_TEST(msreader.getProvidedFields() == kFieldsToRead);

  for (const std::string& setting : kAutoweightSettings) {
    dp3::common::ParameterSet parset_autoweight;
    parset_autoweight.add(setting, "true");
    MSReader autoweight(ms, parset_autoweight, "");
    BOOST_TEST(autoweight.getProvidedFields() == kWeightsField);
    autoweight.setFieldsToRead(kFieldsToRead);
    BOOST_TEST(autoweight.getProvidedFields() ==
               (kFieldsToRead | kWeightsField));
  }
}

// Test reading a LOFAR measurement set
BOOST_AUTO_TEST_CASE(read_lofar) {
  const std::string kMsName = "tNDPPP-generic.MS";
  const casacore::MeasurementSet ms(kMsName);
  dp3::common::ParameterSet parset;
  MSReader msreader(ms, parset, "");

  std::vector<std::string> ant_vec = {"CS001HBA0", "CS002HBA0"};
  std::vector<std::string> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      dp3::common::GetTelescope(
          kMsName, everybeam::ElementResponseModel::kHamaker, false);
  const std::vector<size_t> station_indices =
      MSReader::SelectStationIndices(telescope.get(), ant_names);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(
        GetStationName(phasedarray.GetStation(station_indices[i])), ant_vec[i]);
  }
}

// Reading an OSKAR measurement set
BOOST_AUTO_TEST_CASE(read_oskar) {
  const std::string kMsName = "tOSKAR.in_MS";
  const casacore::MeasurementSet ms(kMsName);
  dp3::common::ParameterSet parset;
  MSReader msreader(ms, parset, "");

  std::vector<std::string> ant_vec = {"s0012", "s0013", "s0015"};
  std::vector<std::string> ant_names(ant_vec.size());

  for (size_t i = 0; i < ant_vec.size(); ++i) {
    ant_names[i] = ant_vec[i];
  }

  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      dp3::common::GetTelescope(
          kMsName, everybeam::ElementResponseModel::kOSKARSphericalWave, false);
  const everybeam::telescope::PhasedArray& phasedarray =
      static_cast<const everybeam::telescope::PhasedArray&>(*telescope);
  std::vector<size_t> station_indices =
      MSReader::SelectStationIndices(telescope.get(), ant_names);

  for (size_t i = 0; i < station_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(
        GetStationName(phasedarray.GetStation(station_indices[i])), ant_vec[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
