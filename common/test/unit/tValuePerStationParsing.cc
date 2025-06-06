// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../ValuePerStationParsing.h"
#include "../../../base/test/LoggerFixture.h"

using dp3::base::test::LoggerFixture;
using dp3::common::ParseValuePerStation;

BOOST_AUTO_TEST_SUITE(value_per_station_parsing)

BOOST_AUTO_TEST_CASE(empty) {
  ParseValuePerStation<double>({}, {}, {});

  std::vector<int> values{1337, 1338};
  const std::vector<std::string> station_names{"station_a", "station_b"};
  ParseValuePerStation<int>(values, {}, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 1338);

  const std::vector<std::string> pattern{"[]:3"};
  ParseValuePerStation<int>(values, pattern, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 1338);
}

BOOST_AUTO_TEST_CASE(assign_single_station) {
  std::vector<int> values{1337, 1338, 1339};
  const std::vector<std::string> station_names{"station_a", "station_b",
                                               "station_c"};
  const std::vector<std::string> pattern_all{"station_a:3", "station_b:7",
                                             "station_c:1"};

  const std::vector pattern_a{pattern_all[0]};
  ParseValuePerStation<int>(values, pattern_a, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 1338);
  BOOST_CHECK_EQUAL(values[2], 1339);

  const std::vector pattern_b{pattern_all[1]};
  ParseValuePerStation<int>(values, pattern_b, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 1339);

  const std::vector pattern_c{pattern_all[2]};
  ParseValuePerStation<int>(values, pattern_c, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 1);

  values = {1337, 1338, 1339};
  ParseValuePerStation<int>(values, pattern_all, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 1);
}

BOOST_AUTO_TEST_CASE(assign_multiple_stations) {
  std::vector<int> values{1337, 1338, 1339};
  const std::vector<std::string> station_names{"station_a", "station_b",
                                               "station_c"};
  const std::vector<std::string> pattern_all{
      "[station_b,station_c]:3", "[station_a,station_c,station_b]:7",
      "[station_c,station_a]:1"};

  const std::vector pattern_1{pattern_all[0]};
  ParseValuePerStation<int>(values, pattern_1, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 3);
  BOOST_CHECK_EQUAL(values[2], 3);

  const std::vector pattern_2{pattern_all[1]};
  ParseValuePerStation<int>(values, pattern_2, station_names);
  BOOST_CHECK_EQUAL(values[0], 7);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 7);

  const std::vector pattern_3{pattern_all[2]};
  ParseValuePerStation<int>(values, pattern_3, station_names);
  BOOST_CHECK_EQUAL(values[0], 1);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 1);

  values = {1337, 1338, 1339};
  ParseValuePerStation<int>(values, pattern_all, station_names);
  BOOST_CHECK_EQUAL(values[0], 1);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 1);
}

BOOST_AUTO_TEST_CASE(assign_with_wildcards) {
  std::vector<int> values{1337, 1338, 1339};
  const std::vector<std::string> station_names{"station_a", "station_b",
                                               "station_c"};
  const std::vector<std::string> pattern_all{
      "[station_*]:3", "[station_*,^station_a]:7", "station_[a-b]*:1"};

  const std::vector pattern_1{pattern_all[0]};
  ParseValuePerStation<int>(values, pattern_1, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 3);
  BOOST_CHECK_EQUAL(values[2], 3);

  const std::vector pattern_2{pattern_all[1]};
  ParseValuePerStation<int>(values, pattern_2, station_names);
  BOOST_CHECK_EQUAL(values[0], 3);
  BOOST_CHECK_EQUAL(values[1], 7);
  BOOST_CHECK_EQUAL(values[2], 7);

  const std::vector pattern_3{pattern_all[2]};
  ParseValuePerStation<int>(values, pattern_3, station_names);
  BOOST_CHECK_EQUAL(values[0], 1);
  BOOST_CHECK_EQUAL(values[1], 1);
  BOOST_CHECK_EQUAL(values[2], 7);

  values = {1337, 1338, 1339};
  ParseValuePerStation<int>(values, pattern_all, station_names);
  BOOST_CHECK_EQUAL(values[0], 1);
  BOOST_CHECK_EQUAL(values[1], 1);
  BOOST_CHECK_EQUAL(values[2], 7);
}

BOOST_AUTO_TEST_CASE(non_existing_stations,
                     *boost::unit_test::fixture<LoggerFixture>()) {
  // The other tests use int: let's use doubles here to increase tested
  // scenarios.
  std::vector<double> values{1337, 1338, 1339};
  const std::vector<std::string> station_names{"a", "b", "c"};
  const std::vector<std::string> pattern_all{"z:3", "[b,z,y]:5", "[z*, c*]:7"};

  const std::vector pattern_1{pattern_all[0]};
  ParseValuePerStation<double>(values, pattern_1, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 1338);
  BOOST_CHECK_EQUAL(values[2], 1339);

  const std::vector pattern_2{pattern_all[1]};
  ParseValuePerStation<double>(values, pattern_2, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 5);
  BOOST_CHECK_EQUAL(values[2], 1339);

  const std::vector pattern_3{pattern_all[2]};
  ParseValuePerStation<double>(values, pattern_3, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 5);
  BOOST_CHECK_EQUAL(values[2], 7);

  values = {1337, 1338, 1339};
  ParseValuePerStation<double>(values, pattern_all, station_names);
  BOOST_CHECK_EQUAL(values[0], 1337);
  BOOST_CHECK_EQUAL(values[1], 5);
  BOOST_CHECK_EQUAL(values[2], 7);
}

BOOST_AUTO_TEST_SUITE_END()
