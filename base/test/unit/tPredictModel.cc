// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../PredictModel.h"

#include <boost/test/unit_test.hpp>

constexpr size_t kNThreads = 3;
constexpr size_t kNCorrelations = 4;
constexpr size_t kNChannels = 128;
constexpr size_t kNBaselines = 10 * 9 / 2;
constexpr bool kIncludeBeam = true;

BOOST_AUTO_TEST_SUITE(predict_model)

BOOST_AUTO_TEST_CASE(construct) {
  dp3::base::PredictModel buffer(kNThreads, kNCorrelations, kNChannels,
                                 kNBaselines, kIncludeBeam);
  for (size_t i = 0; i != kNThreads; ++i) {
    BOOST_CHECK(buffer.GetModel(i).shape()[0] == kNBaselines);
    BOOST_CHECK(buffer.GetModel(i).shape()[1] == kNChannels);
    BOOST_CHECK(buffer.GetModel(i).shape()[2] == kNCorrelations);
    BOOST_CHECK(buffer.GetPatchModel(i).shape()[0] == kNBaselines);
    BOOST_CHECK(buffer.GetPatchModel(i).shape()[1] == kNChannels);
    BOOST_CHECK(buffer.GetPatchModel(i).shape()[2] == kNCorrelations);
  }
}

BOOST_AUTO_TEST_SUITE_END()
