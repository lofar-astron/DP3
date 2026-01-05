// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ddecal/constraints/AntennaIntervalConstraint.h"

#include <complex>
#include <numeric>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <xtensor/xview.hpp>

using dp3::ddecal::AntennaIntervalConstraint;
using dp3::ddecal::Constraint;
using dp3::ddecal::ConstraintResult;

BOOST_AUTO_TEST_SUITE(antenna_interval_constraint)

namespace {
void CheckValue(dp3::ddecal::SolutionTensor& solutions, size_t channel,
                size_t antenna, size_t sub_solution, double expected) {
  BOOST_CHECK_CLOSE_FRACTION(
      solutions(channel, antenna, sub_solution, 0).real(), expected, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(
      solutions(channel, antenna, sub_solution, 1).real(), expected + 1.0,
      1e-6);
}

}  // namespace

BOOST_AUTO_TEST_CASE(single_direction) {
  constexpr size_t kNChannels = 2;
  constexpr size_t kNAntennas = 3;
  constexpr size_t kNSubSolutions = 5;
  constexpr size_t kNPolarizations = 2;

  std::vector<size_t> intervals_per_antenna{1, 3, 2};
  AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));

  constraint.Initialize(kNAntennas, {5u}, {100.0e6, 120.0e6});

  dp3::ddecal::SolutionTensor solutions(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});

  std::iota(solutions.begin(), solutions.end(), 0.0);

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(solutions);
  std::vector<ConstraintResult> constraint_result =
      constraint.Apply(onesolution_span, 0.0, nullptr);
  BOOST_CHECK(constraint_result.empty());

  for (size_t channel = 0; channel != kNChannels; ++channel) {
    const size_t channel_start =
        channel * kNAntennas * kNSubSolutions * kNPolarizations;
    size_t value = channel_start;

    // Antenna 0 is not averaged, so check if it is equal to the input.
    for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
         ++sub_solution) {
      for (size_t polarization = 0; polarization != kNPolarizations;
           ++polarization) {
        BOOST_CHECK_CLOSE_FRACTION(
            solutions(channel, 0, sub_solution, polarization).real(), value,
            1e-6);
        ++value;
      }
    }

    // Antenna 1 subsolutions 0-2 have been averaged over 3 subsolutions, with
    // values 10, 12 and 14.
    CheckValue(solutions, channel, 1, 0, channel_start + 12.0);
    CheckValue(solutions, channel, 1, 1, channel_start + 12.0);
    CheckValue(solutions, channel, 1, 2, channel_start + 12.0);
    // Antenna 1 subsolutions 3-4 have been averaged over 2 subsolutions, with
    // values 16, 18.
    CheckValue(solutions, channel, 1, 3, channel_start + 17.0);
    CheckValue(solutions, channel, 1, 4, channel_start + 17.0);
    // Antenna 2 subsolutions 0-1 have been averaged over 2 subsolutions, with
    // values 20, 22.
    CheckValue(solutions, channel, 2, 0, channel_start + 21.0);
    CheckValue(solutions, channel, 2, 1, channel_start + 21.0);
    // Antenna 2 subsolutions 2-3 over values 24, 26.
    CheckValue(solutions, channel, 2, 2, channel_start + 25.0);
    CheckValue(solutions, channel, 2, 3, channel_start + 25.0);
    // Antenna 2 subsolution 4 is not averaged and has value 28.
    CheckValue(solutions, channel, 2, 4, channel_start + 28.0);
  }
}

BOOST_AUTO_TEST_CASE(multiple_directions) {
  constexpr size_t kNChannels = 2;
  constexpr size_t kNAntennas = 3;
  constexpr size_t kNSubSolutions = 8;
  constexpr size_t kNPolarizations = 2;

  std::vector<size_t> intervals_per_antenna{1, 3, 2};
  AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));

  std::vector<uint32_t> solutions_per_direction{3, 5};
  constraint.Initialize(kNAntennas, solutions_per_direction,
                        {100.0e6, 120.0e6});

  dp3::ddecal::SolutionTensor solutions(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});

  std::iota(solutions.begin(), solutions.end(), 0.0);

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(solutions);
  std::vector<ConstraintResult> constraint_result =
      constraint.Apply(onesolution_span, 0.0, nullptr);
  BOOST_CHECK(constraint_result.empty());

  for (size_t channel = 0; channel != kNChannels; ++channel) {
    const size_t channel_start =
        channel * kNAntennas * kNSubSolutions * kNPolarizations;
    size_t value = channel_start;

    // Direction 0 has 3 subsolutions
    // antenna 0 is not extra averaged, so should have 3 independent values
    // Direction 1 has 5 subsolutions (subsolutions 3-7)
    // antenna 0 is not extra averaged, so should have 5 independent values
    for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
         ++sub_solution) {
      for (size_t polarization = 0; polarization != kNPolarizations;
           ++polarization) {
        BOOST_CHECK_CLOSE_FRACTION(
            solutions(channel, 0, sub_solution, polarization).real(), value,
            1e-6);
        ++value;
      }
    }

    // Direction 0 has 3 subsolutions
    // Antenna 1 subsolutions 0-2 have been averaged over 3 subsolutions, with
    // values 16, 18 and 20.
    CheckValue(solutions, channel, 1, 0, channel_start + 18.0);
    CheckValue(solutions, channel, 1, 1, channel_start + 18.0);
    CheckValue(solutions, channel, 1, 2, channel_start + 18.0);
    // Direction 1 has 5 subsolutions. Averaging this over 3 subsolutions will
    // average over (22, 24, 26) and (28, 30).
    CheckValue(solutions, channel, 1, 3, channel_start + 24.0);
    CheckValue(solutions, channel, 1, 4, channel_start + 24.0);
    CheckValue(solutions, channel, 1, 5, channel_start + 24.0);
    CheckValue(solutions, channel, 1, 6, channel_start + 29.0);
    CheckValue(solutions, channel, 1, 7, channel_start + 29.0);

    // Direction 0 has 3 subsolutions
    // Antenna 2 subsolutions 0-1 have been averaged over 2 subsolutions
    // (32, 34), and 3 is 36.
    CheckValue(solutions, channel, 2, 0, channel_start + 33.0);
    CheckValue(solutions, channel, 2, 1, channel_start + 33.0);
    CheckValue(solutions, channel, 2, 2, channel_start + 36.0);
    // Direction 1 has 5 subsolutions. Averaging this over 2 subsolutions:
    // (38, 40) ; (42, 44) and 46.
    CheckValue(solutions, channel, 2, 3, channel_start + 39.0);
    CheckValue(solutions, channel, 2, 4, channel_start + 39.0);
    CheckValue(solutions, channel, 2, 5, channel_start + 43.0);
    CheckValue(solutions, channel, 2, 6, channel_start + 43.0);
    CheckValue(solutions, channel, 2, 7, channel_start + 46.0);
  }
}

BOOST_AUTO_TEST_CASE(no_averaging) {
  constexpr size_t kNChannels = 7;
  constexpr size_t kNAntennas = 3;
  constexpr size_t kNSubSolutions = 9;
  constexpr size_t kNPolarizations = 4;

  std::vector<size_t> intervals_per_antenna(kNAntennas, 1);
  AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));

  std::vector<uint32_t> solutions_per_direction{3, 1, 4, 1};
  constraint.Initialize(
      kNAntennas, solutions_per_direction,
      {100.0e6, 105.0e6, 110.0e6, 115.0e6, 120.0e6, 125.0e6, 130.0e6});

  dp3::ddecal::SolutionTensor solutions(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});

  std::iota(solutions.begin(), solutions.end(), 0.0);

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(solutions);
  std::vector<ConstraintResult> constraint_result =
      constraint.Apply(onesolution_span, 0.0, nullptr);
  BOOST_CHECK(constraint_result.empty());

  double value = 0.0;
  for (size_t channel = 0; channel != kNChannels; ++channel) {
    for (size_t antenna = 0; antenna != kNAntennas; ++antenna) {
      for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
           ++sub_solution) {
        for (size_t polarization = 0; polarization != kNPolarizations;
             ++polarization) {
          BOOST_CHECK_CLOSE_FRACTION(
              solutions(channel, antenna, sub_solution, polarization).real(),
              value, 1e-6);
          ++value;
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(incorrect_zero_factor) {
  auto function = []() {
    constexpr size_t kNAntennas = 3;
    std::vector<size_t> intervals_per_antenna{1, 0, 2};
    AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));
    std::vector<uint32_t> solutions_per_direction{3, 5};
    constraint.Initialize(kNAntennas, solutions_per_direction,
                          {100.0e6, 120.0e6});
  };
  BOOST_CHECK_THROW(function(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(incorrect_high_factor) {
  auto function = []() {
    constexpr size_t kNAntennas = 3;
    std::vector<size_t> intervals_per_antenna{1, 6, 2};
    AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));
    std::vector<uint32_t> solutions_per_direction{3, 5};
    constraint.Initialize(kNAntennas, solutions_per_direction,
                          {100.0e6, 120.0e6});
  };
  BOOST_CHECK_THROW(function(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(incorrect_number_of_factors) {
  auto function = []() {
    constexpr size_t kNAntennas = 3;
    // intentionally wrong number of elements:
    std::vector<size_t> intervals_per_antenna{1, 6};
    AntennaIntervalConstraint constraint(std::move(intervals_per_antenna));
    std::vector<uint32_t> solutions_per_direction{3, 5};
    constraint.Initialize(kNAntennas, solutions_per_direction,
                          {100.0e6, 120.0e6});
  };
  BOOST_CHECK_THROW(function(), std::runtime_error);
}
BOOST_AUTO_TEST_SUITE_END()
