#include <boost/test/unit_test.hpp>

#include "../../Settings.h"

namespace dp3::ddecal {

BOOST_AUTO_TEST_SUITE(ddecal_settings)

BOOST_AUTO_TEST_CASE(solution_direction_map) {
  std::vector<size_t> result = GetSolutionToDirectionVector({});
  BOOST_CHECK(result.empty());

  result = GetSolutionToDirectionVector({1, 1, 1});
  BOOST_CHECK((result == std::vector<size_t>{0, 1, 2}));

  result = GetSolutionToDirectionVector({3, 1, 1, 2});
  BOOST_CHECK((result == std::vector<size_t>{0, 0, 0, 1, 2, 3, 3}));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace dp3::ddecal
