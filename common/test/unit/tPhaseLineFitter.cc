// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "common/PhaseLineFitter.h"

#include <random>

#include <boost/test/unit_test.hpp>

#include <aocommon/constants.h>

using dp3::common::phase_fitting::FitSample;
using dp3::common::phase_fitting::UnwrapAndFit;

namespace {

std::vector<FitSample> MakeSimpleData(double slope) {
  std::vector<FitSample> data;
  for (size_t i = 0; i != 100; ++i) {
    // add one very bad sample with zero weight to checking weighting
    // All calculations are performed in module 2pi, hence a value of 1e10 mod
    // 2pi is not doing much in the tests with noise. However, in the tests
    // without noise a value of 1e10 mod 2pi should still cause a bias, so this
    // still tests whether weighting works in the exact fits.
    const double x = i + 30;
    if (i == 10)
      data.emplace_back(FitSample{x, 1.0e10, 0.0});
    else
      data.emplace_back(
          FitSample{x, std::fmod(double(x * slope), 2.0 * M_PI), 1.0});
  }
  return data;
}

std::vector<FitSample> MakeFaradayData(double faraday_depth) {
  constexpr double kStartFrequency = 120.0e6;
  constexpr double kBandwidth = 60.0e6;
  constexpr double kNSamples = 250;
  std::vector<FitSample> data;
  for (size_t i = 0; i != kNSamples; ++i) {
    const double frequency =
        kStartFrequency + kBandwidth * (kNSamples - i - 1) / kNSamples;
    const double wavelength = aocommon::kSpeedOfLight / frequency;
    const double x = wavelength * wavelength;
    data.push_back(FitSample{x, std::fmod(x * faraday_depth, 2.0 * M_PI), 1.0});
  }
  return data;
}

std::vector<FitSample> MakeTecData(double slope) {
  constexpr double kStartFrequency = 115.0e6;
  constexpr double kBandwidth = 75.0e6;
  constexpr double kNSamples = 375;
  std::vector<FitSample> data;
  for (size_t i = 0; i != kNSamples; ++i) {
    const double frequency =
        kStartFrequency + kBandwidth * (kNSamples - i - 1) / kNSamples;
    const double x = 1.0 / frequency;
    data.push_back(FitSample{x, std::fmod(x * slope, 2.0 * M_PI), 1.0});
  }
  return data;
}

std::vector<double> MakeTrialSlopes() {
  std::vector<double> slopes;
  const size_t n_trials = 100;
  const double max_slope = 3.0;
  for (size_t i = 0; i != n_trials; ++i) {
    const double slope = i * 2 * max_slope / n_trials - max_slope;
    slopes.emplace_back(slope);
  }
  return slopes;
}

void AddNoise(std::vector<FitSample>& data, double stddev) {
  std::mt19937 rnd;
  std::normal_distribution<float> distribution(0.0, stddev);
  for (FitSample& sample : data) {
    sample.y += distribution(rnd);
  }
}

void AddNoiseWithWeights(std::vector<FitSample>& data, double stddev) {
  std::mt19937 rnd;
  std::normal_distribution<float> distribution(0.0, stddev);
  for (FitSample& sample : data) {
    double noise = distribution(rnd);
    sample.y += noise;
    // Give less weight to noisy values. See weighted_noisy_fit test.
    sample.weight = 1.0 / (noise * noise);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(phase_line_fitter)

BOOST_AUTO_TEST_CASE(exact_unwrap_and_fit) {
  for (size_t i = 0; i != 20; ++i) {
    const double slope = i * 0.1 - 1.0;

    const std::vector<FitSample> data = MakeSimpleData(slope);
    const double result = UnwrapAndFit(data, slope);
    constexpr double kSlopeTolerance = 1.0e-4;
    BOOST_CHECK_CLOSE_FRACTION(result, slope, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_CASE(close_unwrap_and_fit) {
  const std::vector<double> trial_slopes = MakeTrialSlopes();
  for (double slope : trial_slopes) {
    const std::vector<FitSample> data = MakeSimpleData(slope);
    // The max x value is 100, so adding or subtracting 0.01 from the slope
    // should not cause extra wraping
    const double result_a = UnwrapAndFit(data, slope + 0.01);
    const double result_b = UnwrapAndFit(data, slope - 0.01);
    constexpr double kSlopeTolerance = 1.0e-4;
    if (std::fabs(slope) < kSlopeTolerance) {
      BOOST_CHECK_LT(std::abs(result_a), kSlopeTolerance);
      BOOST_CHECK_LT(std::abs(result_b), kSlopeTolerance);
    } else {
      BOOST_CHECK_CLOSE_FRACTION(result_a, slope, kSlopeTolerance);
      BOOST_CHECK_CLOSE_FRACTION(result_b, slope, kSlopeTolerance);
    }
  }
}

BOOST_AUTO_TEST_CASE(noiseless_fit) {
  const std::vector<double> trial_slopes = MakeTrialSlopes();
  for (double slope : trial_slopes) {
    const std::vector<FitSample> data = MakeSimpleData(slope);
    const double result = FitSlope(data);
    constexpr double kSlopeTolerance = 1.0e-4;
    if (std::fabs(slope) < kSlopeTolerance)
      BOOST_CHECK_LT(std::abs(result), kSlopeTolerance);
    else
      BOOST_CHECK_CLOSE_FRACTION(result, slope, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_CASE(noisy_fit) {
  const std::vector<double> trial_slopes = MakeTrialSlopes();
  for (double slope : trial_slopes) {
    std::vector<FitSample> data = MakeSimpleData(slope);
    AddNoise(data, 0.1);
    const double result = FitSlope(data);
    constexpr double kSlopeTolerance = 5.0e-2;
    if (std::fabs(slope) < 1e-4)
      BOOST_CHECK_LT(std::abs(result), kSlopeTolerance);
    else
      BOOST_CHECK_CLOSE_FRACTION(result, slope, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_CASE(weighted_noisy_fit) {
  std::vector<double> trial_slopes = MakeTrialSlopes();
  for (double slope : trial_slopes) {
    std::vector<FitSample> data = MakeSimpleData(slope);
    AddNoiseWithWeights(data, 0.1);
    const double result = FitSlope(data);
    // The weighting downweights values with high noise values. Therefore,
    // these results should be more accurate than the results from the
    // 'noisy_fit' test.
    constexpr double kSlopeTolerance = 1.0e-4;
    if (std::fabs(slope) < kSlopeTolerance)
      BOOST_CHECK_LT(std::abs(result), kSlopeTolerance);
    else
      BOOST_CHECK_CLOSE_FRACTION(result, slope, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_CASE(faraday_function) {
  // The difference between this test and the noisy_fit test, is that a
  // different x axis is used. This test uses a grid that corresponds with a
  // faraday depth fit, i.e. the x-axis is the wavelength^2 value.
  for (double faraday_depth = -100.0; faraday_depth < 100.5;
       faraday_depth += 4) {
    std::vector<FitSample> data = MakeFaradayData(faraday_depth);

    const dp3::common::phase_fitting::SlopeFitRange range_full =
        dp3::common::phase_fitting::GetRange(data);
    BOOST_CHECK_GT(range_full.max_wraps, 100);
    BOOST_CHECK_LT(range_full.wrap_step, 1);
    // The slope must be inside the searched area
    BOOST_CHECK_GT(range_full.wrap_step * range_full.max_wraps,
                   std::abs(faraday_depth));

    constexpr int max_wrap_limit = 200;
    const dp3::common::phase_fitting::SlopeFitRange range_limited =
        dp3::common::phase_fitting::GetRange(
            data, aocommon::OptionalNumber<int>(max_wrap_limit));
    BOOST_CHECK_EQUAL(range_limited.max_wraps, max_wrap_limit);
    BOOST_CHECK_LT(range_limited.wrap_step, 1);
    BOOST_CHECK_GT(range_limited.wrap_step * range_limited.max_wraps,
                   std::abs(faraday_depth));

    AddNoise(data, 0.1);
    const double result = FitSlope(data, range_limited);
    constexpr double kSlopeTolerance = 0.05;
    if (std::fabs(faraday_depth) < 1e-4)
      BOOST_CHECK_LT(std::abs(result), kSlopeTolerance);
    else
      BOOST_CHECK_CLOSE_FRACTION(result, faraday_depth, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_CASE(tec_function) {
  // This value comes from how TEC is defined. I don't know its origin,
  // I copied it from AlphaToTEC() in base/PhaseFitter.h.
  constexpr double kTecToSlopeFactor = -8.44797245e9;
  for (double tec_value = -0.3; tec_value < 0.305; tec_value += 0.03) {
    const double slope = tec_value * kTecToSlopeFactor;
    std::vector<FitSample> data = MakeTecData(slope);

    const dp3::common::phase_fitting::SlopeFitRange range =
        dp3::common::phase_fitting::GetRange(data,
                                             aocommon::OptionalNumber<int>(25));
    BOOST_CHECK_EQUAL(range.max_wraps, 25);
    // The slope must be inside the searched area
    BOOST_CHECK_GT(range.wrap_step * range.max_wraps, std::abs(slope));

    const double result = FitSlope(data, range);
    constexpr double kSlopeTolerance = 0.05;
    if (std::fabs(slope) < kSlopeTolerance)
      BOOST_CHECK_LT(std::abs(result), kSlopeTolerance);
    else
      BOOST_CHECK_CLOSE_FRACTION(result, slope, kSlopeTolerance);
  }
}

BOOST_AUTO_TEST_SUITE_END()
