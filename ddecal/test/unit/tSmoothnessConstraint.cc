// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <vector>
#include <complex>

#include "../../constraints/SmoothnessConstraint.h"

#include <aocommon/xt/span.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <xtensor/xview.hpp>

using dp3::ddecal::SmoothnessConstraint;

namespace {
// We use slightly more than 2 MHz as bandwidth, to make sure the kernel
// is at least 2 MHz wide (and thus covers 5 channels).
const double bandwidth_hz = 2.01e6;
const double bandwidth_ref_frequency_hz = 0.0;
const double spectral_exponent = -1.0;
const double solution_time = 0.0;
const size_t kNChannelBlocks = 5;

SmoothnessConstraint makeConstraint(size_t n_antennas, size_t n_directions) {
  const std::vector<double> frequencies{1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6};
  const std::vector<double> weights(kNChannelBlocks, 1.0);
  std::vector<double> antenna_distance_factors{1.0};

  SmoothnessConstraint c(bandwidth_hz, bandwidth_ref_frequency_hz,
                         spectral_exponent);
  c.Initialize(n_antennas, std::vector<uint32_t>(n_directions, 1), frequencies);
  c.SetAntennaFactors(std::move(antenna_distance_factors));
  c.SetWeights(weights);
  return c;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(smoothness_constraint)

BOOST_AUTO_TEST_CASE(simple_case) {
  const size_t kNAntennas = 1;
  const size_t kNDirections = 3;
  const size_t kNPolarizations = 2;

  // The columns form the polarizations and directions.
  // Each polarization (2x) and each direction (3x) is individually smoothed
  // The first three columns test smoothing a delta function.
  // The fourth and fifth columns are constant, and when the input is constant,
  // the smoother should not modify it.
  // The sixth column tests smoothing a fixed slope.
  // The imaginary part of the solution is always zero in this test.

  dp3::ddecal::SolutionTensor solutions{
      {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 5.0}}},  // channel block 0
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 4.0}}},  // channel block 1
      {{{1.0, 0.0}, {0.0, 1.0}, {0.0, 3.0}}},  // channel block 2
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}}},  // channel block 3
      {{{0.0, 0.0}, {10.0, 1.0}, {0.0, 1.0}}}  // channel block 4
  };
  BOOST_REQUIRE(solutions.shape(0) == kNChannelBlocks);
  BOOST_REQUIRE(solutions.shape(1) == kNAntennas);
  BOOST_REQUIRE(solutions.shape(2) == kNDirections);
  BOOST_REQUIRE(solutions.shape(3) == kNPolarizations);

  SmoothnessConstraint constraint = makeConstraint(kNAntennas, kNDirections);
  dp3::ddecal::SolutionSpan solutions_span =
      aocommon::xt::CreateSpan(solutions);
  constraint.Apply(solutions_span, solution_time, nullptr);

  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{0.000121798, 0.902597}, {0.0, 1.0}, {0.0, 4.90247}}},
      {{{0.0886568, 0.0886568}, {0.0, 1.0}, {0.0, 3.99978}}},
      {{{0.822484, 0.000110987}, {0.00110987, 1.0}, {0.0, 3}}},
      {{{0.0886568, 0.0}, {0.886568, 1.0}, {0.0, 2.00022}}},
      {{{0.000121798, 0.0}, {9.02597, 1.0}, {0.0, 1.09753}}},
  };

  BOOST_CHECK(xt::allclose(solutions, reference_solutions));
}

BOOST_AUTO_TEST_CASE(dd_factors) {
  const size_t kNAntennas = 1;
  const size_t kNDirections = 3;

  dp3::ddecal::SolutionTensor solutions{
      {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 5.0}}},  // channel block 0
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 4.0}}},  // channel block 1
      {{{1.0, 0.0}, {0.0, 1.0}, {0.0, 3.0}}},  // channel block 2
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}}},  // channel block 3
      {{{0.0, 0.0}, {10.0, 1.0}, {0.0, 1.0}}}  // channel block 4
  };

  SmoothnessConstraint constraint = makeConstraint(kNAntennas, kNDirections);
  std::vector<double> factors{4.0, 2.0, 1.0};
  BOOST_REQUIRE_EQUAL(factors.size(), kNDirections);
  constraint.SetDdSmoothingFactors(std::move(factors));
  dp3::ddecal::SolutionSpan solutions_span =
      aocommon::xt::CreateSpan(solutions);
  constraint.Apply(solutions_span, solution_time, nullptr);

  // In comparison with the 'simple_case' test above, the values
  // in the first and second column are smoothed by a higher value
  // (4x and 2x, respectively).
  // TODO it seems that the kernel smoother does not enlarge its
  // window when a value > 1 is specified (hence the zero values
  // in the second column).
  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{0.209986, 0.366484}, {0.0, 1.0}, {0.0, 4.90247}}},
      {{{0.262608, 0.262608}, {0.0, 1.0}, {0.0, 3.99978}}},
      {{{0.257334, 0.147445}, {0.456402, 1.0}, {0.0, 3}}},
      {{{0.262608, 0.0}, {2.542337, 1.0}, {0.0, 2.00022}}},
      {{{0.234536, 0.0}, {5.949716, 1.0}, {0.0, 1.09753}}},
  };

  auto result_iter = solutions_span.begin();
  for (const std::complex<double>& ref : reference_solutions) {
    BOOST_CHECK_CLOSE_FRACTION(ref.real(), result_iter->real(), 1e-3);
    BOOST_CHECK_LT(std::abs(ref.imag()), 1e-3);
    ++result_iter;
  }

  BOOST_CHECK(xt::allclose(solutions, reference_solutions));
}

BOOST_AUTO_TEST_CASE(flagged_channels) {
  const size_t kNAntennas = 1;
  const size_t kNDirections = 6;
  // Using a single polarization makes XTensor ignore the last index value.
  // Indexing using (ch, 0, 0, solution_index), where solution_index > 1,
  // does not work anymore in that case.
  const size_t kNPolarizations = 1;

  SmoothnessConstraint constraint = makeConstraint(kNAntennas, kNDirections);
  // Flagged channels should be ignored by the constraint, but still
  // provide a value in the output.
  // The first two channels are flagged:
  constraint.SetWeights({0.0, 0.0, 1.0, 1.0, 1.0});

  dp3::ddecal::SolutionTensor solutions{
      {{{0.0}, {1.0}, {0.0}, {1.0}, {0.0}, {5.0}}},  // channel block 0
      {{{0.0}, {0.0}, {0.0}, {1.0}, {0.0}, {4.0}}},  // channel block 1
      {{{1.0}, {0.0}, {0.0}, {1.0}, {0.0}, {3.0}}},  // channel block 2
      {{{0.0}, {0.0}, {0.0}, {1.0}, {0.0}, {2.0}}},  // channel block 3
      {{{0.0}, {0.0}, {10.0}, {1.0}, {0.0}, {1.0}}}  // channel block 4
  };
  BOOST_REQUIRE(solutions.shape(0) == kNChannelBlocks);
  BOOST_REQUIRE(solutions.shape(1) == kNAntennas);
  BOOST_REQUIRE(solutions.shape(2) == kNDirections);
  BOOST_REQUIRE(solutions.shape(3) == kNPolarizations);

  dp3::ddecal::SolutionSpan solutions_span =
      aocommon::xt::CreateSpan(solutions);
  constraint.Apply(solutions_span, solution_time, nullptr);

  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{1.0}, {0.0}, {0.0}, {1.0}, {0.0}, {3.0}}},
      {{{0.99875}, {0.0}, {0.0, 1.0}, {0.0}, {2.99875}}},
      {{{0.902597}, {0.0}, {0.00121798}, {1.0}, {0.0}, {2.90247}}},
      {{{0.0886666}, {0.0}, {0.886666}, {1.0}, {0.0}, {2.0}}},
      {{{0.000121798}, {0.0}, {9.02597}, {1.0}, {0.0}, {1.09753}}},
  };

  BOOST_CHECK(xt::allclose(solutions, reference_solutions));
}

BOOST_AUTO_TEST_SUITE_END()
