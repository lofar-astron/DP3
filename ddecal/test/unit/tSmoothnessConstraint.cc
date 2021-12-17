// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <vector>
#include <complex>

#include "../../constraints/SmoothnessConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::SmoothnessConstraint;

namespace {
// We use slightly more than 2 MHz as bandwidth, to make sure the kernel
// is at least 2 MHz wide (and thus covers 5 channels).
const double bandwidth_hz = 2.01e6;
const double bandwidth_ref_frequency_hz = 0.0;
const double solution_time = 0.0;
const size_t n_antennas = 1;
const size_t n_directions = 3;
const size_t n_solution_polarizations = 2;
const size_t n_channelblocks = 5;

SmoothnessConstraint makeConstraint() {
  const std::vector<double> frequencies{1e6, 2e6, 3e6, 4e6, 5e6};
  std::vector<double> weights{1, 1, 1, 1, 1};
  std::vector<double> antenna_distance_factors{1};

  SmoothnessConstraint c(bandwidth_hz, bandwidth_ref_frequency_hz);
  c.SetNThreads(1);
  c.Initialize(n_antennas, std::vector<uint32_t>(n_directions, 1), frequencies);
  c.SetDistanceFactors(std::move(antenna_distance_factors));
  c.SetWeights(weights);
  return c;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(smoothness_constraint)

BOOST_AUTO_TEST_CASE(simple_cases) {
  // The columns form the polarizations and directions.
  // Each polarization (2x) and each direction (3x) is individually smoothed
  // The first three columns test smoothing a delta function.
  // The fourth and fifth columns are constant, and when the input is constant,
  // the smoother should not modify it.
  std::vector<std::vector<std::complex<double>>> solutions{
      {0.0, 1.0, 0.0, 1.0, 0.0, 5.0},  // channel block 0
      {0.0, 0.0, 0.0, 1.0, 0.0, 4.0},  // channel block 1
      {1.0, 0.0, 0.0, 1.0, 0.0, 3.0},  // channel block 2
      {0.0, 0.0, 0.0, 1.0, 0.0, 2.0},  // channel block 3
      {0.0, 0.0, 10.0, 1.0, 0.0, 1.0}  // channel block 4
  };

  SmoothnessConstraint constraint = makeConstraint();
  constraint.Apply(solutions, solution_time, nullptr);

  const std::vector<std::vector<std::complex<double>>> reference_solutions{
      {0.000121798, 0.902597, 0.0, 1.0, 0.0, 4.90247},
      {0.0886568, 0.0886568, 0.0, 1.0, 0.0, 3.99978},
      {0.822484, 0.000110987, 0.00110987, 1.0, 0.0, 3},
      {0.0886568, 0.0, 0.886568, 1.0, 0.0, 2.00022},
      {0.000121798, 0.0, 9.02597, 1.0, 0.0, 1.09753},
  };

  for (size_t pol_dir = 0; pol_dir != n_solution_polarizations * n_directions;
       ++pol_dir) {
    for (size_t cb = 0; cb != n_channelblocks; ++cb) {
      BOOST_CHECK_CLOSE(solutions[cb][pol_dir].real(),
                        reference_solutions[cb][pol_dir].real(), 1e-3);
      BOOST_CHECK_CLOSE(solutions[cb][pol_dir].imag(),
                        reference_solutions[cb][pol_dir].imag(), 1e-3);
    }
  }
}

BOOST_AUTO_TEST_CASE(flagged_channels) {
  SmoothnessConstraint constraint = makeConstraint();
  // Flagged channels should be ignored by the constraint, but still
  // provide a value in the output.
  // The first two channels are flagged:
  constraint.SetWeights({0.0, 0.0, 1.0, 1.0, 1.0});

  std::vector<std::vector<std::complex<double>>> solutions{
      {0.0, 1.0, 0.0, 1.0, 0.0, 5.0},  // channel block 0
      {0.0, 0.0, 0.0, 1.0, 0.0, 4.0},  // channel block 1
      {1.0, 0.0, 0.0, 1.0, 0.0, 3.0},  // channel block 2
      {0.0, 0.0, 0.0, 1.0, 0.0, 2.0},  // channel block 3
      {0.0, 0.0, 10.0, 1.0, 0.0, 1.0}  // channel block 4
  };

  constraint.Apply(solutions, solution_time, nullptr);

  const std::vector<std::vector<std::complex<double>>> reference_solutions{
      {1.0, 0.0, 0.0, 1.0, 0.0, 3.0},
      {0.99875, 0.0, 0.0, 1.0, 0.0, 2.99875},
      {0.902597, 0.0, 0.00121798, 1.0, 0.0, 2.90247},
      {0.0886666, 0.0, 0.886666, 1.0, 0.0, 2},
      {0.000121798, 0.0, 9.02597, 1.0, 0.0, 1.09753},
  };

  for (size_t pol_dir = 0; pol_dir != n_solution_polarizations * n_directions;
       ++pol_dir) {
    for (size_t cb = 0; cb != n_channelblocks; ++cb) {
      BOOST_CHECK_CLOSE(solutions[cb][pol_dir].real(),
                        reference_solutions[cb][pol_dir].real(), 1e-3);
      BOOST_CHECK_CLOSE(solutions[cb][pol_dir].imag(),
                        reference_solutions[cb][pol_dir].imag(), 1e-3);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
