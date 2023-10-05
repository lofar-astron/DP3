// tBdaDdeCal.cc: Test program for class BdaDdeCal
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <aocommon/recursivefor.h>

#include "tPredict.h"
#include "tStepCommon.h"

#include "../../BdaDdeCal.h"
#include <dp3/common/Fields.h>
#include "../../../common/ParameterSet.h"

using dp3::steps::BdaDdeCal;

namespace {

dp3::common::ParameterSet CreateMinimalParameterSet() {
  dp3::common::ParameterSet parset;
  parset.add("msin", "");
  parset.add("directions", "[[" + dp3::steps::test::kPredictDirection + "]]");
  parset.add("sourcedb", dp3::steps::test::kPredictSourceDB);
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
  const dp3::common::Fields kExpectedFields = dp3::steps::Step::kUvwField;
  BOOST_TEST(bdaddecal->getRequiredFields() == kExpectedFields);
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

BOOST_AUTO_TEST_SUITE_END()
