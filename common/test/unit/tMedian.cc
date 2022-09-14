// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Median.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(median)

BOOST_AUTO_TEST_CASE(empty) {
  std::vector<int> data = {};
  BOOST_CHECK_EQUAL(dp3::common::Median(data), 0);
}

BOOST_AUTO_TEST_CASE(unique) {
  for (int n = 1; n <= 10; ++n) {
    std::vector<int> data(n);
    for (int i = 0; i < n; ++i) {
      data[i] = i + 1;
    }

    do {
      std::vector<int> copy = data;
      int median = dp3::common::Median(copy);
      if (n % 2 == 0) {
        BOOST_CHECK_EQUAL(median, n / 2);
      } else {
        BOOST_CHECK_EQUAL(median, n / 2 + 1);
      }
    } while (std::next_permutation(data.begin(), data.end()));
  }
}

BOOST_AUTO_TEST_SUITE_END()
