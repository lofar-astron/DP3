// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <aocommon/threadpool.h>

#include "tPredict.h"
#include "tStepCommon.h"

#include "../../BdaDdeCal.h"
#include <dp3/common/Fields.h>
#include "../../../common/ParameterSet.h"

using dp3::steps::BdaDdeCal;

namespace {

dp3::common::ParameterSet CreateMinimalPararameterSet() {
  dp3::common::ParameterSet parset;
  parset.add("msin", "");
  parset.add("directions", "[[" + dp3::steps::test::kPredictDirection + "]]");
  parset.add("sourcedb", dp3::steps::test::kPredictSourceDB);
  parset.add("h5parm", "test.h5");
  return parset;
}

/// Fixture that creates a minimal working BdaDdeCal step.
struct BdaDdeCalFixture {
  BdaDdeCalFixture() : input(), bdaddecal() {
    dp3::common::ParameterSet parset = CreateMinimalPararameterSet();
    bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(&input, parset, "");
  }

  dp3::steps::MockInput input;
  std::shared_ptr<BdaDdeCal> bdaddecal;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(bdaddecal)

/// Helper forwarder to improve test failure output.
static void TestShow(
    const std::string& expected,
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  dp3::steps::MockInput input;
  const dp3::common::ParameterSet parset =
      dp3::steps::test::CreateParameterSet(parameters);
  const BdaDdeCal bdaddecal(&input, parset, "prefix.");
  BOOST_CHECK_EQUAL(expected, dp3::steps::test::Show(bdaddecal));

  BOOST_CHECK_MESSAGE(parset.unusedKeys().empty(),
                      "Not all keys are used, is there a typo?");
}

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

BOOST_AUTO_TEST_CASE(show_default) {
  TestShow(
      R"(BdaDdeCal prefix.
  mode (constraints):  diagonal
  directions:          [[center]]
  solver algorithm:    directionsolve
  H5Parm:              tDDECal.MS/instrument.h5
  subtract model:      false
  solution interval:   0 s
  #channels/block:     1
  #channel blocks:     18446744073709551615
  tolerance:           0.0001
  max iter:            50
  flag unconverged:    false
     diverged only:    false
  propagate solutions: false
       converged only: false
  detect stalling:     true
  step size:           0.2
Model steps for direction [center]
Predict
BDAExpander prefix.
OnePredict prefix.
  sourcedb:                tDDECal.MS/sky.txt
   number of patches:      1
   patches clustered:      false
   number of components:   1
   absolute orientation:   false
   all unpolarized:        true
   correct freq smearing:  false
  apply beam:              false
  operation:               replace
  threads:                 )" +
          std::to_string(aocommon::system::ProcessorCount()) + R"(
BDAAverager prefix.
  timebase:        0s
  max interval:    0s
  frequencybase:   0
  min channels:    1
  max time factor: 1
  max freq factor: 1

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky.txt"}});
}

BOOST_AUTO_TEST_CASE(show_modified) {
  TestShow(
      R"(BdaDdeCal prefix.
  mode (constraints):  diagonal
  directions:          [[center]]
  solver algorithm:    hybrid
  H5Parm:              tDDECal.MS/instrument.h5
  subtract model:      true
  solution interval:   0 s
  #channels/block:     44
  #channel blocks:     18446744073709551615
  tolerance:           1e-05
  max iter:            49
  flag unconverged:    true
     diverged only:    true
  propagate solutions: true
       converged only: true
  detect stalling:     true
  step size:           0.2
  coreconstraint:      45.123
  smoothnessconstraint:46.123
  smoothnessreffrequency:47.123
  smoothnessrefdistance:48.123
  tecscreen.coreconstraint:49.123
Model steps for direction [center]
Predict
BDAExpander prefix.
OnePredict prefix.
  sourcedb:                tDDECal.MS/sky.txt
   number of patches:      1
   patches clustered:      false
   number of components:   1
   absolute orientation:   false
   all unpolarized:        true
   correct freq smearing:  false
  apply beam:              false
  operation:               replace
  threads:                 )" +
          std::to_string(aocommon::system::ProcessorCount()) + R"(
BDAAverager prefix.
  timebase:        0s
  max interval:    0s
  frequencybase:   0
  min channels:    1
  max time factor: 1
  max freq factor: 1

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.propagatesolutions", "true"},
       {"prefix.propagateconvergedonly", "true"},
       {"prefix.flagunconverged", "true"},
       {"prefix.flagdivergedonly", "true"},
       {"prefix.subtract", "true"},
       {"prefix.solveralgorithm", "hybrid"},
       {"prefix.sourcedb", "tDDECal.MS/sky.txt"},
       {"prefix.nchan", "44"},
       {"prefix.coreconstraint", "45.123"},
       {"prefix.smoothnessconstraint", "46.123"},
       {"prefix.smoothnessreffrequency", "47.123"},
       {"prefix.smoothnessrefdistance", "48.123"},
       {"prefix.tecscreen.coreconstraint", "49.123"},
       {"prefix.maxiter", "49"}});
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
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset = CreateMinimalPararameterSet();
  parset.add("correcttimesmearing", "10");
  auto bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(&input, parset, "");

  const dp3::common::Fields kExpectedFields =
      Step::kFlagsField | Step::kWeightsField | Step::kFullResFlagsField |
      Step::kUvwField;
  BOOST_TEST(bdaddecal->getRequiredFields() == kExpectedFields);
}

BOOST_FIXTURE_TEST_CASE(get_provided_fields, BdaDdeCalFixture) {
  const dp3::common::Fields kExpectedFields = {};
  BOOST_TEST(bdaddecal->getProvidedFields() == kExpectedFields);
}

BOOST_AUTO_TEST_CASE(get_provided_fields_subtract) {
  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset = CreateMinimalPararameterSet();
  parset.add("subtract", "true");
  auto bdaddecal = std::make_shared<dp3::steps::BdaDdeCal>(&input, parset, "");

  const dp3::common::Fields kExpectedFields = dp3::steps::Step::kDataField;
  BOOST_TEST(bdaddecal->getProvidedFields() == kExpectedFields);
}

BOOST_AUTO_TEST_SUITE_END()
