// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_PHASE_LINE_FITTER_H_
#define DP3_COMMON_PHASE_LINE_FITTER_H_

#include <cmath>
#include <limits>
#include <vector>

#include <aocommon/optionalnumber.h>

namespace dp3::common::phase_fitting {

struct FitSample {
  /**
   * The grid x value, e.g. the wavelength^2 value when fitting the
   * Faraday value.
   */
  double x;
  /**
   * The phase value of this sample in radians.
   */
  double y;
  double weight;
};

struct SlopeFitRange {
  int max_wraps;
  double wrap_step;
};

inline double SlopeModelCost(const std::vector<FitSample>& data, double slope) {
  double cost = 0.0;
  double weight_sum = 0.0;
  for (const FitSample& sample : data) {
    const double y_model = slope * sample.x;
    double sample_cost = std::fmod(std::fabs(y_model - sample.y), 2.0 * M_PI);
    if (sample_cost > M_PI) sample_cost = 2.0 * M_PI - sample_cost;
    sample_cost *= sample.weight;
    cost += sample_cost;
    weight_sum += sample.weight;
  }
  if (weight_sum == 0.0)
    return 0.0;
  else
    return cost / weight_sum;
}

inline double UnwrapAndFit(const std::vector<FitSample>& data,
                           double unwrapping_slope) {
  // Minimize: sum w_i (slope * x_i - y_i)^2
  // Therefore:
  // sum w_i 2 x_i (slope * x_i - y_i) = 0
  // slope = (sum w_i x_i y_i) / (sum w_i x_i^2)
  double numerator = 0.0;
  double denominator = 0.0;
  for (const FitSample& sample : data) {
    const double y_model = unwrapping_slope * sample.x;
    double difference = std::fmod(sample.y - y_model, 2.0 * M_PI);
    if (difference >= M_PI)
      difference -= 2.0 * M_PI;
    else if (difference < -M_PI)
      difference += 2.0 * M_PI;
    const double unwrapped_y = y_model + difference;
    const double w_times_x = sample.x * sample.weight;
    numerator += w_times_x * unwrapped_y;
    denominator += w_times_x * sample.x;
  }
  if (denominator != 0.0)
    return numerator / denominator;
  else
    return 0.0;
}

/**
 * Computes the range of slopes that needs to be searched through. This only
 * uses the @c x values of the @p data parameter, and the result can thus be
 * reused for multiple fits if the @c x values stay the same. This function
 * requires the data to be sorted by x-value.
 *
 * @param wrap_count_limit if set, limits the number of wraps that are searched
 * through. The number of wraps is counted from the origin up to the largest @c
 * x value. If unset, the number of wraps to be searched through is calculated
 * from the data, such that any wrapping that can be resolved, is resolved. This
 * number can be very high for dense x-values far from zero (e.g. when fitting a
 * slope between 150-151 MHz with 1 kHz resolution), hence for these cases it
 * might be useful to limit the number of wraps searched through.
 */
SlopeFitRange GetRange(const std::vector<FitSample>& data,
                       aocommon::OptionalNumber<int> wrap_count_limit = {}) {
  if (data.empty()) return SlopeFitRange{0, 0.0};

  double min_delta = std::numeric_limits<double>::max();
  for (size_t i = 0; i != data.size() - 1; ++i) {
    const double delta = data[i + 1].x - data[i].x;
    assert(delta >= 0.0);  // 'data' should be ordered on 'x'.
    if (delta > 0.0 && delta < min_delta) min_delta = delta;
  }
  const double min_x = data.front().x;
  const double max_x = data.back().x;
  const double max_abs_x = max_x < 0.0 ? -min_x : max_x;

  // The search must go from 0 wraps to the maximum number of wraps that can
  // still lead to finding the right wrapping.
  // In theory, it is still possible to find the right slope value when the
  // slope is so steep that there is a wrap in the small grid distance. This can
  // lead to a large nr of wraps; e.g. in an extreme case of fitting a slope
  // through two values around 150 MHz separated by 1 Hz it would cause 150e6
  // values to search through.
  int max_wraps = std::ceil(max_abs_x / min_delta);

  // There is another upper bound to the max nr of wraps: with a low nr of
  // values, there is a fixed, low nr of ways that the values can have been
  // wrapped. This leads to an exponential bound (~ 2^n_samples). However, in
  // this case, the search space is not uniform, so we add a bit more space. It
  // is only relevant for a low nr of samples.
  if (data.size() < 20) {
    constexpr int extra_space = 10;
    max_wraps = std::min(max_wraps, extra_space * (2 << data.size()));
  }

  if (wrap_count_limit) {
    max_wraps = std::min(max_wraps, *wrap_count_limit);
  }

  // To not miss any possible correct unwrapping, we need to step through every
  // slope for which the largest x-value increases by 2pi. Given y = slope * x,
  // we can dy = dslope * x. Therefore, dslope = dy / max_abs_x = 2 pi /
  // max_abs_x. However, with 2pi steps it seems some slopes are missed (unit
  // tests fail), but with 'pi' steps it works.
  const double wrap_step = M_PI / max_abs_x;

  return SlopeFitRange{max_wraps, wrap_step};
}

/**
 * Same as the other @ref FitSlope() overload, but with a parameter to provide
 * the range. If performing the same fit multiple times, the range can be
 * calculated once for all fits.
 */
inline double FitSlope(const std::vector<FitSample>& data,
                       const SlopeFitRange& range) {
  double best_cost = std::numeric_limits<double>::max();
  double best_estimate = 0.0;
  for (int i = -range.max_wraps; i <= range.max_wraps; ++i) {
    double slope = range.wrap_step * i;
    const double estimate = UnwrapAndFit(data, slope);
    const double cost = SlopeModelCost(data, estimate);
    if (cost < best_cost) {
      best_cost = cost;
      best_estimate = estimate;
    }
  }
  // Perform one more fit iteration on the best estimate. Because a better
  // estimate is known, the wrapping of noisy values might slightly change, so
  // refit.
  return UnwrapAndFit(data, best_estimate);
}

/**
 * Fits phase data to a line, accounting for phase wraps. It
 * optimizes slope s to approximately minimize the least squares
 * value between the data and the model function:
 * y = ( s x ) modulo 2pi
 * by minimizing the cost function:
 * w_i (y_i - s * x_i)^2
 * This function requires the data to be sorted by x-value.
 * The data should have finite values (even for zero-weighted
 * data).
 *
 * While searching for the right unwrapping, the function minimizes the
 * sum of absolute values instead of least squares. This is slightly
 * more stable. It does a brute for search of all possible wraps.
 * Once the right unwrapping is found, the values are unwrapped and
 * a least-squares fit of a line to the data is performed.
 *
 * Because the grid of x-value does not need to be regular, this function
 * is suitable for finding:
 * - delay values by using frequency as x-value;
 * - TEC values by using wavelength or 1/nu as x-value;
 * - or Faraday-depth values by using wavelength^2 as x-value.
 */
inline double FitSlope(const std::vector<FitSample>& data) {
  return FitSlope(data, GetRange(data));
}

}  // namespace dp3::common::phase_fitting

#endif
