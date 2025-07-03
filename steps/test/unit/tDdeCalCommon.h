// Common test code for tDdeCal.cc and tBdaDdeCal.cc.
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "tPredict.h"
#include "tStepCommon.h"

namespace dp3::steps::test::ddecal {

namespace {

/// Creates a basic parameter set for testing DdeCal / BdaDdeCal, when
/// the step gets input from tNDPPP-generic.MS / tNDPPP-bda.MS.
dp3::common::ParameterSet DdeParameterSet() {
  return CreateParameterSet(
      {{"directions", "[[" + dp3::steps::test::kPredictDirection + "]]"},
       {"sourcedb", dp3::steps::test::kPredictSkymodel},
       {"h5parm", "test.h5"},
       {"solint", "2"},
       {"nchan", "42"}});
}

}  // namespace

template <class DdeStep>
void TestInfoDirectionsWithoutReuse(dp3::steps::Step& reader,
                                    bool keep_model_data) {
  reader.setInfo(dp3::base::DPInfo());

  dp3::common::ParameterSet parset = DdeParameterSet();
  if (keep_model_data) parset.add("keepmodel", "true");

  auto bdaddecal = std::make_shared<DdeStep>(parset, "");
  BOOST_TEST(parset.unusedKeys().empty());
  const dp3::base::DPInfo& info_out = bdaddecal->setInfo(reader.getInfoOut());

  const std::map<std::string, dp3::base::Direction>& directions =
      info_out.GetDirections();
  if (keep_model_data) {
    BOOST_TEST_REQUIRE(directions.size() == 1);
    BOOST_TEST(directions.begin()->first ==
               dp3::steps::test::kPredictDirection);
  } else {
    BOOST_TEST_REQUIRE(directions.empty());
  }
}

template <class DdeStep>
void TestInfoDirectionsWithReuse(dp3::steps::Step& reader,
                                 bool keep_model_data) {
  const std::vector<std::string> kReusePatterns = {
      "[dir3,dir1,dir2]", "[*]", "[d?r?]", "[dir1]", "[*2]", "[???3]"};

  // When keep_model_data is true, the output DPInfo should contain everything.
  const std::vector<std::string> kAllDirections = {
      dp3::steps::test::kPredictDirection, "dir1", "dir2", "dir3"};

  // When keep_model_data is false, the output DPInfo should only contain the
  // directions that do not match the reuse pattern.
  const std::vector<std::vector<std::string>> kRemainingDirections = {
      {}, {}, {}, {"dir2", "dir3"}, {"dir1", "dir3"}, {"dir1", "dir2"}};

  reader.setInfo(dp3::base::DPInfo());

  for (std::size_t i = 0; i < kReusePatterns.size(); ++i) {
    dp3::common::ParameterSet parset = DdeParameterSet();
    parset.add("keepmodel", keep_model_data ? "true" : "false");
    parset.add("reusemodel", kReusePatterns[i]);

    dp3::base::DPInfo info_in = reader.getInfoOut();
    std::map<std::string, dp3::base::Direction>& input_directions =
        info_in.GetDirections();
    input_directions.emplace("dir1", dp3::base::Direction());
    input_directions.emplace("dir2", dp3::base::Direction());
    input_directions.emplace("dir3", dp3::base::Direction());

    auto bdaddecal = std::make_shared<DdeStep>(parset, "");
    BOOST_TEST(parset.unusedKeys().empty());
    const dp3::base::DPInfo& info_out = bdaddecal->setInfo(info_in);

    const std::map<std::string, dp3::base::Direction>& output_directions =
        info_out.GetDirections();
    const std::vector<std::string>& expected_directions =
        keep_model_data ? kAllDirections : kRemainingDirections[i];

    BOOST_TEST_REQUIRE(output_directions.size() == expected_directions.size());
    auto output_direction = output_directions.begin();
    for (const std::string& expected_direction : expected_directions) {
      BOOST_TEST(output_direction->first == expected_direction);
      ++output_direction;
    }
  }
}

}  // namespace dp3::steps::test::ddecal