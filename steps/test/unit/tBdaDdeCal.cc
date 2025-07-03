// tBdaDdeCal.cc: Test program for class BdaDdeCal
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../BdaDdeCal.h"

#include <fstream>

#include <boost/test/unit_test.hpp>

#include <dp3/common/Fields.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"
#include "../../MSBDAReader.h"
#include "../../MultiResultStep.h"

#include "tDdeCalCommon.h"

using dp3::base::BdaBuffer;
using dp3::common::Fields;
using dp3::common::test::FixtureDirectory;
using dp3::steps::BdaDdeCal;
using dp3::steps::BDAResultStep;
using dp3::steps::MSBDAReader;
using dp3::steps::test::kTrueFalseRange;

namespace {

dp3::common::ParameterSet CreateMinimalParameterSet() {
  dp3::common::ParameterSet parset;
  parset.add("directions", "[[" + dp3::steps::test::kPredictDirection + "]]");
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("h5parm", "test.h5");
  return parset;
}

/// Fixture that creates a minimal working BdaDdeCal step.
struct BdaDdeCalFixture {
  BdaDdeCalFixture() : bdaddecal() {
    dp3::common::ParameterSet parset = CreateMinimalParameterSet();
    bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(parset, "");
  }

  std::shared_ptr<BdaDdeCal> bdaddecal;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(bdaddecal)

void CheckChannelBlockIndex(std::shared_ptr<dp3::steps::BdaDdeCal>& bdaddecal,
                            size_t n_channels, size_t n_channel_blocks,
                            const std::vector<size_t>& expected_result) {
  for (size_t i = 0; i < n_channels; ++i) {
    if (n_channels >= n_channel_blocks) {
      BOOST_TEST(bdaddecal->GetChanBlockIndex(
                     i, n_channels, n_channel_blocks) == expected_result[i]);
    } else {
      BOOST_REQUIRE_THROW(
          bdaddecal->GetChanBlockIndex(i, n_channels, n_channel_blocks),
          std::runtime_error);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(channel_block_mapping_16_channels, BdaDdeCalFixture) {
  size_t n_channel_blocks = 4;
  size_t n_channels = 16;
  std::vector<size_t> expected_result = {0, 0, 0, 0, 1, 1, 1, 1,
                                         2, 2, 2, 2, 3, 3, 3, 3};
  CheckChannelBlockIndex(bdaddecal, n_channels, n_channel_blocks,
                         expected_result);
}

BOOST_FIXTURE_TEST_CASE(channel_block_mapping_10_channels, BdaDdeCalFixture) {
  size_t n_channel_blocks = 4;
  size_t n_channels = 10;
  std::vector<size_t> expected_result = {0, 0, 0, 1, 1, 2, 2, 2, 3, 3};

  CheckChannelBlockIndex(bdaddecal, n_channels, n_channel_blocks,
                         expected_result);
}

BOOST_FIXTURE_TEST_CASE(channel_block_mapping_4_channels, BdaDdeCalFixture) {
  size_t n_channel_blocks = 4;
  size_t n_channels = 4;
  std::vector<size_t> expected_result = {0, 1, 2, 3};
  CheckChannelBlockIndex(bdaddecal, n_channels, n_channel_blocks,
                         expected_result);
}

BOOST_FIXTURE_TEST_CASE(channel_block_mapping_3_channels, BdaDdeCalFixture) {
  size_t n_channel_blocks = 4;
  size_t n_channels = 3;
  std::vector<size_t> expected_result = {};
  CheckChannelBlockIndex(bdaddecal, n_channels, n_channel_blocks,
                         expected_result);
}

BOOST_FIXTURE_TEST_CASE(get_required_fields, BdaDdeCalFixture) {
  using dp3::steps::Step;
  const dp3::common::Fields kExpectedFields =
      Step::kFlagsField | Step::kWeightsField | Step::kUvwField;
  BOOST_TEST(bdaddecal->getRequiredFields() == kExpectedFields);
}

BOOST_AUTO_TEST_CASE(get_required_fields_only_predict) {
  dp3::common::ParameterSet parset = CreateMinimalParameterSet();
  parset.add("onlypredict", "true");
  const BdaDdeCal bdaddecal(parset, "");

  const dp3::common::Fields kExpectedFields = dp3::steps::Step::kUvwField;
  BOOST_TEST(bdaddecal.getRequiredFields() == kExpectedFields);
}

BOOST_AUTO_TEST_CASE(get_required_fields_correct_time_smearing) {
  using dp3::steps::Step;
  dp3::common::ParameterSet parset = CreateMinimalParameterSet();
  parset.add("correcttimesmearing", "10");
  auto bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(parset, "");

  const dp3::common::Fields kExpectedFields =
      Step::kFlagsField | Step::kWeightsField | Step::kUvwField;
  BOOST_TEST(bdaddecal->getRequiredFields() == kExpectedFields);
}

BOOST_FIXTURE_TEST_CASE(get_provided_fields, BdaDdeCalFixture) {
  const dp3::common::Fields kExpectedFields = {};
  BOOST_TEST(bdaddecal->getProvidedFields() == kExpectedFields);
}

BOOST_AUTO_TEST_CASE(get_provided_fields_subtract) {
  dp3::common::ParameterSet parset = CreateMinimalParameterSet();
  parset.add("subtract", "true");
  auto bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(parset, "");

  const dp3::common::Fields kExpectedFields = dp3::steps::Step::kDataField;
  BOOST_TEST(bdaddecal->getProvidedFields() == kExpectedFields);
}

struct BdaMsFixture : public FixtureDirectory {
  BdaMsFixture() : FixtureDirectory() {
    ExtractResource("tNDPPP-bda.MS.tgz");
    casacore::MeasurementSet ms("tNDPPP-bda.MS");
    const dp3::common::ParameterSet kEmptyParset;
    reader = std::make_shared<MSBDAReader>(ms, kEmptyParset, "");
  }

  /// Reader for the extracted MS. Tests can make it generate a valid DPInfo
  /// object for BdaDdeCal or use it in a step chain, for example.
  std::shared_ptr<MSBDAReader> reader;
};

BOOST_DATA_TEST_CASE_F(BdaMsFixture, info_directions_no_reuse, kTrueFalseRange,
                       keep_model_data) {
  using dp3::steps::test::ddecal::TestInfoDirectionsWithoutReuse;
  TestInfoDirectionsWithoutReuse<BdaDdeCal>(*reader, keep_model_data);
}

BOOST_DATA_TEST_CASE_F(BdaMsFixture, info_directions_with_reuse,
                       kTrueFalseRange, keep_model_data) {
  using dp3::steps::test::ddecal::TestInfoDirectionsWithReuse;
  TestInfoDirectionsWithReuse<BdaDdeCal>(*reader, keep_model_data);
}

BOOST_DATA_TEST_CASE_F(BdaMsFixture, keep_or_discard_model_data,
                       kTrueFalseRange* kTrueFalseRange, only_predict,
                       keep_model_data) {
  // Generate a skymodel for testing, similar to tBdaDdeCal.py.
  const std::string kSkymodelFileName = "test.skymodel";
  std::ofstream skymodel_file(kSkymodelFileName);
  skymodel_file
      << "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, "
         "PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\n"
      << "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \n"
      << "ra_off, POINT, 16:58:28.205000, +63.44.34.314000, 1, , , , , \n"
      << "radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 1, , , , , \n";
  skymodel_file.close();

  dp3::common::ParameterSet parset;
  if (only_predict) {
    parset.add("ddecal.onlypredict", "true");
  } else {
    parset.add("ddecal.h5parm", "test.h5parm");
  }
  if (keep_model_data) parset.add("ddecal.keepmodel", "true");

  parset.add("ddecal.sourcedb", kSkymodelFileName);
  parset.add("ddecal.solint", "2");
  parset.add("ddecal.nchan", "42");

  auto ddecal = std::make_shared<BdaDdeCal>(parset, "ddecal.");
  BOOST_TEST(parset.unusedKeys().empty());

  reader->setFieldsToRead(ddecal->getRequiredFields());
  auto output = std::make_shared<BDAResultStep>();
  dp3::steps::test::Execute({reader, ddecal, output});

  std::vector<std::string> expected_data_names;
  if (only_predict) {
    expected_data_names.push_back("");
  }
  if (keep_model_data) {
    expected_data_names.push_back("ddecal.center");
    expected_data_names.push_back("ddecal.ra_off");
    expected_data_names.push_back("ddecal.radec_off");
  }

  std::vector<std::unique_ptr<BdaBuffer>> result = output->Extract();
  BOOST_TEST(!result.empty());
  for (std::unique_ptr<BdaBuffer>& buffer : result) {
    const std::vector<std::string> data_names = buffer->GetDataNames();
    BOOST_CHECK_EQUAL_COLLECTIONS(data_names.begin(), data_names.end(),
                                  expected_data_names.begin(),
                                  expected_data_names.end());
    if (only_predict) {
      // When only_predict is enabled, flags and weights are not in
      // BdaDdeCal::getRequiredFields(), so the reader does not read them, and
      // BdaDdeCal also does not output these fields.
      BOOST_TEST(!buffer->GetFlags());
      BOOST_TEST(!buffer->GetWeights());
    } else {
      // When only predict is disabled, BdaDdeCal should output the flags and
      // weights it received as input.
      BOOST_TEST(buffer->GetFlags());
      BOOST_TEST(buffer->GetWeights());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
