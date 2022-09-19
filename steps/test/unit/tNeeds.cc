// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Needs.h"

#include <boost/test/unit_test.hpp>

using dp3::steps::Needs;

BOOST_AUTO_TEST_SUITE(needs)

namespace {
void CheckNeeds(const Needs& needs, bool data, bool flags, bool weights,
                bool full_res_flags, bool uvw) {
  BOOST_CHECK(needs.Data() == data);
  BOOST_CHECK(needs.Flags() == flags);
  BOOST_CHECK(needs.Weights() == weights);
  BOOST_CHECK(needs.FullResFlags() == full_res_flags);
  BOOST_CHECK(needs.Uvw() == uvw);
}
}  // namespace

BOOST_AUTO_TEST_CASE(constructor_default) {
  const Needs needs;
  CheckNeeds(needs, false, false, false, false, false);
}

BOOST_AUTO_TEST_CASE(constructor_single_need) {
  const Needs data(Needs::Single::kData);
  const Needs flags(Needs::Single::kFlags);
  const Needs weights(Needs::Single::kWeights);
  const Needs full_res_flags(Needs::Single::kFullResFlags);
  const Needs uvw(Needs::Single::kUvw);
  CheckNeeds(data, true, false, false, false, false);
  CheckNeeds(flags, false, true, false, false, false);
  CheckNeeds(weights, false, false, true, false, false);
  CheckNeeds(full_res_flags, false, false, false, true, false);
  CheckNeeds(uvw, false, false, false, false, true);
}

BOOST_AUTO_TEST_CASE(combine) {
  Needs needs;
  needs |= Needs(Needs::Single::kFlags);
  CheckNeeds(needs, false, true, false, false, false);
  needs |= Needs(Needs::Single::kFullResFlags);
  CheckNeeds(needs, false, true, false, true, false);

  const Needs data(Needs::Single::kData);
  const Needs uvw(Needs::Single::kUvw);
  CheckNeeds(data | uvw, true, false, false, false, true);
}

BOOST_AUTO_TEST_SUITE_END()
