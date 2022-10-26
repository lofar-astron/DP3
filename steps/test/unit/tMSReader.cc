// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSReader.h"

#include <vector>

#include <boost/test/unit_test.hpp>

#include "../../../common/ParameterSet.h"

using dp3::steps::MSReader;

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

BOOST_AUTO_TEST_SUITE_END()
