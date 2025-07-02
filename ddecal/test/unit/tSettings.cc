// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Settings.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace dp3::ddecal {

BOOST_AUTO_TEST_SUITE(ddecal_settings)

BOOST_AUTO_TEST_CASE(solution_direction_map) {
  std::vector<size_t> result = GetSolutionToDirectionVector({});
  BOOST_TEST(result.empty());

  result = GetSolutionToDirectionVector({1, 1, 1});
  BOOST_TEST((result == std::vector<size_t>{0, 1, 2}));

  result = GetSolutionToDirectionVector({3, 1, 1, 2});
  BOOST_TEST((result == std::vector<size_t>{0, 0, 0, 1, 2, 3, 3}));
}

struct InputDirectionsFixture {
  InputDirectionsFixture() {
    // Test prefix removal using names with and without dots.
    input_directions.emplace("prefix.dots.foo", base::Direction());
    input_directions.emplace("without_dot", base::Direction());
    parset.add("h5parm", "foo.h5");  // Mandatory setting.
  }

  std::map<std::string, base::Direction> input_directions;
  common::ParameterSet parset;
};

BOOST_FIXTURE_TEST_CASE(get_reused_directions_no_input_directions,
                        InputDirectionsFixture) {
  input_directions.clear();

  const Settings no_reuse(parset, "");
  BOOST_TEST(no_reuse.GetReusedDirections(input_directions).empty());

  parset.add("reusemodel", "[foo]");
  const Settings reuse_foo(parset, "");
  BOOST_TEST(parset.unusedKeys().empty());
  BOOST_CHECK_THROW(reuse_foo.GetReusedDirections(input_directions),
                    std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(get_reused_directions_match_one,
                        InputDirectionsFixture) {
  parset.add("reusemodel", "[prefix.dots.foo]");
  const Settings settings(parset, "");

  const std::vector<std::pair<std::string, std::string>> reused_directions =
      settings.GetReusedDirections(input_directions);
  BOOST_TEST_REQUIRE(reused_directions.size() == 1);
  BOOST_TEST(reused_directions[0].first == "prefix.dots.foo");
  BOOST_TEST(reused_directions[0].second == "dots.foo");
}

BOOST_DATA_TEST_CASE_F(InputDirectionsFixture, get_reused_directions_match_two,
                       boost::unit_test::data::make(
                           {"[prefix.dots.foo,without_dot]",  // Same order.
                            "[without_dot,prefix.dots.foo]",  // Opposite order.
                            "[*]"}),  // Use match expression.
                       patterns) {
  parset.add("reusemodel", patterns);
  const Settings settings(parset, "");

  const std::vector<std::pair<std::string, std::string>> reused_directions =
      settings.GetReusedDirections(input_directions);
  BOOST_TEST_REQUIRE(reused_directions.size() == 2);
  // GetReusedDirections processes directions in the order of the input
  // directions, which is sorted, since it's a map. The return value is
  // therefore also sorted.
  BOOST_TEST(reused_directions[0].first == "prefix.dots.foo");
  BOOST_TEST(reused_directions[0].second == "dots.foo");
  BOOST_TEST(reused_directions[1].first == "without_dot");
  BOOST_TEST(reused_directions[1].second == "without_dot");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace dp3::ddecal
