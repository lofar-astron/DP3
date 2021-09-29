// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "tPredict.h"
#include "../../BdaGroupPredict.h"
#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::BdaGroupPredict;

namespace {
class BdaPredictFixture {
 public:
  BdaPredictFixture() : input_(), parset_(), predict_() {
    parset_.add("sourcedb", dp3::steps::test::kPredictSourceDB);
    predict_ = std::make_shared<BdaGroupPredict>(input_, parset_, "");
  }

  void SetInfo() {
    dp3::base::DPInfo info;
    info.init(1, 0, 1, 10, 0.0, 1.0, "", "");

    const std::vector<int> kAnt1{0, 0, 1};
    const std::vector<int> kAnt2{1, 2, 2};
    const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
    const std::vector<double> kAntDiam(3, 1.0);
    const std::vector<casacore::MPosition> kAntPos(3);
    info.set(kAntNames, kAntDiam, kAntPos, kAnt1, kAnt2);

    std::vector<std::vector<double>> chan_freqs{
        {10.0e6}, {9.0e6, 10.0e6, 11.0e6}, {10.0e6}};
    std::vector<std::vector<double>> chan_widths{
        {3.0e6}, {1.0e6, 1.0e6, 1.0e6}, {3.0e6}};

    info.set(std::move(chan_freqs), std::move(chan_widths));
    predict_->setInfo(info);
  }

 protected:
  dp3::steps::MockInput input_;
  dp3::common::ParameterSet parset_;
  std::shared_ptr<dp3::steps::BdaGroupPredict> predict_;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(bdapredict)

BOOST_FIXTURE_TEST_CASE(constructor, BdaPredictFixture) {
  // Nothing to do: The fixture calls the constructor.
}

BOOST_FIXTURE_TEST_CASE(getfirstdirection_noinfo, BdaPredictFixture) {
  BOOST_CHECK_THROW(predict_->GetFirstDirection(), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(getfirstdirection, BdaPredictFixture) {
  SetInfo();

  const dp3::base::Direction first_direction = predict_->GetFirstDirection();

  BOOST_CHECK_CLOSE(first_direction.ra,
                    dp3::steps::test::kExpectedFirstDirection.ra, 1.0e-3);
  BOOST_CHECK_CLOSE(first_direction.dec,
                    dp3::steps::test::kExpectedFirstDirection.dec, 1.0e-3);
}

BOOST_AUTO_TEST_SUITE_END()
