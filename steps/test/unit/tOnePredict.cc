// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "tPredict.h"
#include "../../OnePredict.h"
#include "../../../common/ParameterSet.h"

#include "mock/MockInput.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::OnePredict;

namespace {
class OnePredictFixture {
 public:
  OnePredictFixture() : input_(), predict_() {
    dp3::common::ParameterSet parset;
    parset.add("sourcedb", dp3::steps::test::kPredictSourceDB);
    predict_ = std::make_shared<OnePredict>(&input_, parset, "");
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

    std::vector<double> chan_freqs(1, 10.0e6);
    std::vector<double> chan_widths(1, 3.0e6);

    info.set(std::move(chan_freqs), std::move(chan_widths));
    predict_->setInfo(info);
  }

 protected:
  dp3::steps::MockInput input_;
  std::shared_ptr<dp3::steps::OnePredict> predict_;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(onepredict)

BOOST_FIXTURE_TEST_CASE(constructor, OnePredictFixture) {
  // Nothing to do: The fixture calls the constructor.
}

BOOST_FIXTURE_TEST_CASE(getfirstdirection, OnePredictFixture) {
  const std::pair<double, double> first_direction =
      predict_->GetFirstDirection();

  BOOST_CHECK_CLOSE(first_direction.first,
                    dp3::steps::test::kExpectedFirstDirection.first, 1.0e-3);
  BOOST_CHECK_CLOSE(first_direction.second,
                    dp3::steps::test::kExpectedFirstDirection.second, 1.0e-3);
}

BOOST_AUTO_TEST_SUITE_END()
