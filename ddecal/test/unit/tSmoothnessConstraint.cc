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

SmoothnessConstraint MakeConstraint(size_t n_antennas, size_t n_directions,
                                    size_t n_channel_blocks, bool truncate) {
  const std::vector<double> frequencies{1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6};
  const std::vector<double> weights(n_channel_blocks, 1.0);
  std::vector<double> antenna_distance_factors{1.0};

  SmoothnessConstraint c(bandwidth_hz, bandwidth_ref_frequency_hz,
                         spectral_exponent, truncate);
  c.Initialize(n_antennas, std::vector<uint32_t>(n_directions, 1), frequencies);
  c.SetAntennaFactors(std::move(antenna_distance_factors));
  c.SetWeights(weights);
  return c;
}

void TestConstraint(dp3::ddecal::SolutionTensor& input_solutions,
                    const dp3::ddecal::SolutionTensor& reference_solutions,
                    const std::vector<double>& weights,
                    const std::vector<double>& dd_factors,
                    bool truncate_kernel) {
  const size_t n_channel_blocks = input_solutions.shape(0);
  const size_t n_antennas = input_solutions.shape(1);
  const size_t n_directions = input_solutions.shape(2);
  // size_t n_polarizations = input_solutions.shape(3);

  SmoothnessConstraint constraint = MakeConstraint(
      n_antennas, n_directions, n_channel_blocks, truncate_kernel);
  // Flagged channels should be ignored by the constraint, but still
  // provide a value in the output.
  // The first two channels are flagged:
  if (!weights.empty()) constraint.SetWeights(weights);
  if (!dd_factors.empty()) {
    BOOST_REQUIRE_EQUAL(dd_factors.size(), n_directions);
    constraint.SetDdSmoothingFactors(dd_factors);
  }

  dp3::ddecal::SolutionSpan solutions_span =
      aocommon::xt::CreateSpan(input_solutions);
  constraint.Apply(solutions_span, solution_time, nullptr);

  auto input_iterator = input_solutions.begin();
  for (const std::complex<double>& reference : reference_solutions) {
    BOOST_CHECK_CLOSE_FRACTION(reference.real(), input_iterator->real(), 1e-4);
    BOOST_CHECK_CLOSE_FRACTION(reference.imag(), input_iterator->imag(), 1e-4);
    ++input_iterator;
  }

  BOOST_CHECK(xt::allclose(input_solutions, reference_solutions));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(smoothness_constraint)

BOOST_AUTO_TEST_CASE(simple_case) {
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
  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{0.000121798, 0.902597}, {0.0, 1.0}, {0.0, 4.90247}}},
      {{{0.0886568, 0.0886568}, {0.0, 1.0}, {0.0, 3.99978}}},
      {{{0.822484, 0.000110987}, {0.00110987, 1.0}, {0.0, 3}}},
      {{{0.0886568, 0.0}, {0.886568, 1.0}, {0.0, 2.00022}}},
      {{{0.000121798, 0.0}, {9.02597, 1.0}, {0.0, 1.09753}}},
  };

  TestConstraint(solutions, reference_solutions, {}, {}, true);
}

BOOST_AUTO_TEST_CASE(dd_factors) {
  dp3::ddecal::SolutionTensor solutions{
      {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 5.0}}},  // channel block 0
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 4.0}}},  // channel block 1
      {{{1.0, 0.0}, {0.0, 1.0}, {0.0, 3.0}}},  // channel block 2
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}}},  // channel block 3
      {{{0.0, 0.0}, {10.0, 1.0}, {0.0, 1.0}}}  // channel block 4
  };

  // In comparison with the 'simple_case' test above, the values
  // in the first and second column are smoothed by a higher value
  // (4x and 2x, respectively).
  // The zero values in the second column are because the kernel
  // is truncated (see the no_truncation test). This is quite
  // extreme in this example, nice for testing/visualization, but
  // note that ddecal actually doesn't allow dd-weights > 1 to
  // avoid this issue.
  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{0.209986, 0.366484}, {0.0, 1.0}, {0.0, 4.90247}}},
      {{{0.262608, 0.262608}, {0.0, 1.0}, {0.0, 3.99978}}},
      {{{0.257334, 0.147445}, {0.456402, 1.0}, {0.0, 3}}},
      {{{0.262608, 0.0}, {2.542337, 1.0}, {0.0, 2.00022}}},
      {{{0.234536, 0.0}, {5.949716, 1.0}, {0.0, 1.09753}}},
  };

  const std::vector<double> factors{4.0, 2.0, 1.0};
  TestConstraint(solutions, reference_solutions, {}, factors, true);
}

BOOST_AUTO_TEST_CASE(no_truncation) {
  dp3::ddecal::SolutionTensor solutions{
      {{{0.0, 1.0}, {0.0, 1.0}, {0.0, 5.0}}},  // channel block 0
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 4.0}}},  // channel block 1
      {{{1.0, 0.0}, {0.0, 1.0}, {0.0, 3.0}}},  // channel block 2
      {{{0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}}},  // channel block 3
      {{{0.0, 0.0}, {10.0, 1.0}, {0.0, 1.0}}}  // channel block 4
  };

  // As in the dd_factors test, we use the dd weights, as these
  // make it easier to evalute the truncation. The reference can
  // be compared to the dd_factors test to see the effect of truncation.
  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{0.202006, 0.352558}, {0.000799632, 1.0}, {0.0, 4.90247}}},
      {{{0.241765, 0.241765}, {0.0294459, 1.0}, {0.0, 3.99978}}},
      {{{0.257334, 0.147445}, {0.456402, 1.0}, {0.0, 3}}},
      {{{0.241765, 0.0793708}, {2.534851, 1.0}, {0.0, 2.00022}}},
      {{{0.202006, 0.0379986}, {5.925774, 1.0}, {0.0, 1.09753}}},
  };

  const std::vector<double> factors{4.0, 2.0, 1.0};
  TestConstraint(solutions, reference_solutions, {}, factors, false);
}

BOOST_AUTO_TEST_CASE(flagged_channels) {
  // Using a single polarization makes XTensor ignore the last index value.
  // Indexing using (ch, 0, 0, solution_index), where solution_index > 1,
  // does not work anymore in that case.
  dp3::ddecal::SolutionTensor solutions{
      {{{0.0}, {1.0}, {0.0}, {1.0}, {0.0}, {5.0}}},  // channel block 0
      {{{0.0}, {0.0}, {0.0}, {1.0}, {0.0}, {4.0}}},  // channel block 1
      {{{1.0}, {0.0}, {0.0}, {1.0}, {0.0}, {3.0}}},  // channel block 2
      {{{0.0}, {0.0}, {0.0}, {1.0}, {0.0}, {2.0}}},  // channel block 3
      {{{0.0}, {0.0}, {10.0}, {1.0}, {0.0}, {1.0}}}  // channel block 4
  };

  const dp3::ddecal::SolutionTensor reference_solutions{
      {{{1.0}, {0.0}, {0.0}, {1.0}, {0.0}, {3.0}}},
      {{{0.99875}, {0.0}, {0.0, 1.0}, {0.0}, {2.99875}}},
      {{{0.902597}, {0.0}, {0.00121798}, {1.0}, {0.0}, {2.90247}}},
      {{{0.0886666}, {0.0}, {0.886666}, {1.0}, {0.0}, {2.0}}},
      {{{0.000121798}, {0.0}, {9.02597}, {1.0}, {0.0}, {1.09753}}},
  };

  // Flagged channels should be ignored by the constraint, but still
  // provide a value in the output.
  // The first two channels are flagged:
  const std::vector<double> weights{0.0, 0.0, 1.0, 1.0, 1.0};
  TestConstraint(solutions, reference_solutions, weights, {}, true);
}

BOOST_AUTO_TEST_SUITE_END()
