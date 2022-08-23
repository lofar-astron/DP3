#include "../../baseline_indices.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(baselineutils)

BOOST_AUTO_TEST_CASE(nbaselines) {
  for (size_t n_antennas = 2; n_antennas < 10; ++n_antennas) {
    const size_t n_baselines = n_antennas * (n_antennas - 1) / 2 + n_antennas;
    BOOST_TEST(n_baselines == dp3::common::ComputeNBaselines(n_antennas));
  }
}

/*
 * Helper function to test the row-major baseline utilities
 */
void testComputeBaselineRowMajor(size_t n_antennas) {
  const dp3::common::BaselineOrder order =
      dp3::common::BaselineOrder::kRowMajor;
  // Iterate all baselines
  size_t baseline_count = 0;
  for (size_t antenna1 = 0; antenna1 < n_antennas; ++antenna1) {
    for (size_t antenna2 = antenna1; antenna2 < n_antennas; ++antenna2) {
      // Check reference and computed baseline
      const std::pair<size_t, size_t> baseline_reference = {antenna1, antenna2};
      const std::pair<size_t, size_t> baseline =
          dp3::common::ComputeBaseline(baseline_count, n_antennas, order);
      BOOST_TEST(baseline_reference.first == baseline.first);
      BOOST_TEST(baseline_reference.second == baseline.second);

      // Check baseline indices
      const size_t baseline_index = dp3::common::ComputeBaselineIndex(
          antenna1, antenna2, n_antennas, order);
      BOOST_TEST(baseline_count == baseline_index);
      ++baseline_count;
    }

    // Check baseline indices for antenna1, it should be n_antennas
    const std::vector<size_t> baselineIndices =
        dp3::common::ComputeBaselineList(antenna1, n_antennas, order);
    BOOST_TEST(baselineIndices.size() == n_antennas);

    // Use the baseline indices to lookup the baseline and check whether
    // antenna1 occurs exactly n_antennas + 1 times (correlations + 1
    // autocorrelation).
    size_t antenna_count = 0;
    for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
      const size_t baseline_index = dp3::common::ComputeBaselineIndex(
          antenna1, antenna2, n_antennas, order);
      BOOST_TEST(baselineIndices[antenna2] == baseline_index);
      const std::pair<size_t, size_t> baseline =
          dp3::common::ComputeBaseline(baseline_index, n_antennas, order);
      antenna_count += baseline.first == antenna1;
      antenna_count += baseline.second == antenna1;
    }
    BOOST_TEST(antenna_count == n_antennas + 1);
  }
}

BOOST_AUTO_TEST_CASE(rowmajor) {
  testComputeBaselineRowMajor(5);
  testComputeBaselineRowMajor(10);
}

/*
 * Helper function to test the column-major baseline utilities
 */
void testComputeBaselineColumnMajor(size_t n_antennas) {
  const dp3::common::BaselineOrder order =
      dp3::common::BaselineOrder::kColumnMajor;
  // Iterate all baselines
  size_t baseline_count = 0;
  for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
    for (size_t antenna1 = 0; antenna1 <= antenna2; ++antenna1) {
      // Check reference and computed baseline
      std::pair<size_t, size_t> baseline_reference = {antenna1, antenna2};
      std::pair<size_t, size_t> baseline =
          dp3::common::ComputeBaseline(baseline_count, n_antennas, order);
      BOOST_TEST(baseline_reference.first == baseline.first);
      BOOST_TEST(baseline_reference.second == baseline.second);

      // Check baseline indices
      size_t baseline_index = dp3::common::ComputeBaselineIndex(
          antenna1, antenna2, n_antennas, order);
      BOOST_TEST(baseline_count == baseline_index);
      ++baseline_count;
    }

    // Check baseline indices for antenna2, it should be n_antennas
    const std::vector<size_t> baseline_indices =
        dp3::common::ComputeBaselineList(antenna2, n_antennas, order);
    BOOST_TEST(baseline_indices.size() == n_antennas);

    // Use the baseline indices to lookup the baseline and check whether
    // antenna1 occurs exactly n_antennas + 1 times (correlations + 1
    // autocorrelation).
    size_t antenna_count = 0;
    for (size_t antenna1 = 0; antenna1 < n_antennas; ++antenna1) {
      const size_t baseline_index = dp3::common::ComputeBaselineIndex(
          antenna1, antenna2, n_antennas, order);
      BOOST_TEST(baseline_indices[antenna1] == baseline_index);
      const std::pair<size_t, size_t> baseline =
          dp3::common::ComputeBaseline(baseline_index, n_antennas, order);
      antenna_count += baseline.first == antenna2;
      antenna_count += baseline.second == antenna2;
    }
    BOOST_TEST(antenna_count == n_antennas + 1);
  }
}

BOOST_AUTO_TEST_CASE(colmajor) {
  testComputeBaselineColumnMajor(5);
  testComputeBaselineColumnMajor(10);
}

BOOST_AUTO_TEST_SUITE_END()