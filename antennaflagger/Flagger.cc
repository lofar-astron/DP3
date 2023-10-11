// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <set>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xmasked_view.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xview.hpp>

// See Flagger.h for algorithm details.
#include "Flagger.h"

namespace {
xt::xtensor<int, 1> FindOutliers(float sigma, size_t max_iterations,
                                 xt::xtensor<float, 1>&& data) {
  // All non finite numbers are flagged as outlier
  xt::xtensor<int, 1> flags = !xt::isfinite(data);

  // If all the data is infinite, return the flags
  // without detecting outliers
  if (xt::sum(flags)() == static_cast<int>(flags.size())) {
    return flags;
  }

  // Take the median of the finite numbers only, otherwise the median is NaN.
  const float median = xt::median(xt::filter(data, !flags));
  // Overwrite NaN's with median to workaround limitation of xt::stddev.
  xt::masked_view(data, flags) = median;

  for (size_t i = 0; i < max_iterations; ++i) {
    const float stddev = xt::stddev(data)();

    size_t outlier_count = 0;
    const float lower_bound = median - (sigma * stddev);
    const float upper_bound = median + (sigma * stddev);

    for (size_t j = 0; j < data.size(); ++j) {
      const float value = data[j];
      if (!flags(j) && (value < lower_bound || value > upper_bound)) {
        flags(j) = true;
        data(j) = median;
        ++outlier_count;
      }
    }

    if (outlier_count == 0) {
      break;
    }
  }

  return flags;
}
}  // namespace

namespace dp3 {
namespace antennaflagger {

xt::xtensor<std::complex<float>, 2> Flagger::ComputeStatsStdDev(
    const dp3::base::DPBuffer::DataType& data) {
  const size_t n_baselines = data.shape(0);
  const size_t n_correlations = data.shape(2);
  xt::xtensor<std::complex<float>, 2> stats({n_baselines, n_correlations});
  xt::real(stats) = xt::stddev(xt::real(data), {1});
  xt::imag(stats) = xt::stddev(xt::imag(data), {1});
  return stats;
}

xt::xtensor<std::complex<float>, 2> Flagger::ComputeStatsSumSquare(
    const dp3::base::DPBuffer::DataType& data) {
  const size_t n_baselines = data.shape(0);
  const size_t n_correlations = data.shape(2);
  xt::xtensor<std::complex<float>, 2> stats({n_baselines, n_correlations});
  xt::real(stats) = xt::sum(xt::square(xt::real(data)), {1});
  xt::imag(stats) = xt::sum(xt::square(xt::imag(data)), {1});
  return stats;
}

xt::xtensor<std::complex<float>, 3> Flagger::GroupStats(
    size_t n_stations, size_t n_antennas_per_station,
    const xt::xtensor<std::complex<float>, 2>& stats_baseline,
    common::BaselineOrder baseline_order) {
  const size_t n_correlations = stats_baseline.shape(1);
  const size_t n_antennas = n_stations * n_antennas_per_station;

  xt::xtensor<std::complex<float>, 3> stats_antenna(
      {n_stations, n_antennas_per_station, n_correlations});

  for (size_t station = 0; station < n_stations; ++station) {
    for (size_t antenna1 = 0; antenna1 < n_antennas_per_station; ++antenna1) {
      const size_t antenna2 = station * n_antennas_per_station + antenna1;
      const std::vector<size_t> baseline_indices =
          common::ComputeBaselineList(antenna2, n_antennas, baseline_order);
      xt::view(stats_antenna, station, antenna1, xt::all()) += xt::nansum(
          xt::view(stats_baseline, xt::keep(baseline_indices), xt::all()), 0);
    }
  }

  stats_antenna /= float(n_antennas_per_station);

  return stats_antenna;
}

xt::xtensor<int, 2> Flagger::ComputeAntennaFlags(
    float sigma, int max_iterations,
    const xt::xtensor<std::complex<float>, 3>& stats) {
  const size_t n_stations = stats.shape(0);
  const size_t n_antennas_per_station = stats.shape(1);
  const size_t n_correlations = stats.shape(2);

  xt::xtensor<int, 2> flags({n_stations, n_antennas_per_station}, false);

  for (size_t station = 0; station < n_stations; ++station) {
    for (size_t cor = 0; cor < n_correlations; ++cor) {
      const xt::xtensor<std::complex<float>, 1> check_sum =
          xt::log(1.0f + xt::view(stats, station, xt::all(), cor));

      xt::view(flags, station, xt::all()) |=
          FindOutliers(sigma, max_iterations, xt::real(check_sum)) |
          FindOutliers(sigma, max_iterations, xt::imag(check_sum));
    }
  }

  return flags;
}

xt::xtensor<int, 2> Flagger::ComputeStationFlags(
    float sigma, int max_iterations,
    const xt::xtensor<std::complex<float>, 3>& stats) {
  const size_t n_stations = stats.shape(0);
  const size_t n_correlations = stats.shape(2);

  // In case of n_correlations == 4, use only XX and YY (index 0 and 3),
  // otherwise use only XX (index 0).
  const size_t n_correlations_out = n_correlations == 4 ? 2 : 1;
  xt::xtensor<int, 2> flags({n_stations, n_correlations_out}, false);

  const xt::xtensor<float, 3> stats_real = xt::real(stats);
  const xt::xtensor<float, 3> stats_imag = xt::imag(stats);

  for (size_t cor = 0; cor < n_correlations_out; ++cor) {
    // Process XX and YY, skip XY and YX, thus:
    // when n_correlations == 4:
    //  - n_correlations_out: 2
    //  - cor_in: [0, 3]
    //  - cor_out: [0, 1]
    // otherwise, n_correlations != 4:
    //  - n_correlations_out: 1
    //  - cor_in: [0]
    //  - cor_out: [0]
    const size_t cor_in = cor * 3;
    const size_t cor_out = cor;

    // xt::median with axis argument does not support strided views, a
    // copy of the stats for the current correlation is therefore required
    const xt::xtensor<float, 2> stats_real_cor =
        xt::view(stats_real, xt::all(), xt::all(), cor_in);
    const xt::xtensor<float, 2> stats_imag_cor =
        xt::view(stats_imag, xt::all(), xt::all(), cor_in);

    const float scale_real = xt::median(stats_real_cor);
    const float scale_imag = xt::median(stats_imag_cor);

    xt::xtensor<float, 1> median_real =
        xt::median(stats_real_cor, 1) / scale_real;
    xt::xtensor<float, 1> median_imag =
        xt::median(stats_imag_cor, 1) / scale_imag;

    xt::view(flags, xt::all(), cor_out) =
        FindOutliers(sigma, max_iterations, std::move(median_real)) |
        FindOutliers(sigma, max_iterations, std::move(median_imag));
  }

  return flags;
}

void Flagger::ComputeStats(const dp3::base::DPBuffer::DataType& data,
                           common::BaselineOrder baseline_order) {
  compute_statistics_timer_.start();
  xt::xtensor<std::complex<float>, 2> stats_baseline_stddev =
      ComputeStatsStdDev(data);
  xt::xtensor<std::complex<float>, 2> stats_baseline_sum_square =
      ComputeStatsSumSquare(data);
  stats_antenna_stddev_ =
      Flagger::GroupStats(n_stations_, n_antennas_per_station_,
                          stats_baseline_stddev, baseline_order);
  stats_antenna_sum_square_ =
      Flagger::GroupStats(n_stations_, n_antennas_per_station_,
                          stats_baseline_sum_square, baseline_order);
  compute_statistics_timer_.stop();
}

void Flagger::AssertStatsComputed() const {
  assert(stats_antenna_stddev_.shape(0) == n_stations_);
  assert(stats_antenna_stddev_.shape(1) == n_antennas_per_station_);
  assert(stats_antenna_stddev_.shape(2) == n_correlations_);
  assert(stats_antenna_stddev_.shape() == stats_antenna_sum_square_.shape());
}

xt::xtensor<int, 1> Flagger::FindBadAntennas(float sigma, int max_iterations) {
  AssertStatsComputed();

  find_bad_antennas_timer_.start();

  const xt::xtensor<int, 2> flags_stddev = Flagger::ComputeAntennaFlags(
      sigma, max_iterations, stats_antenna_stddev_);
  assert(flags_stddev.shape(0) == n_stations_);
  assert(flags_stddev.shape(1) == n_antennas_per_station_);

  const xt::xtensor<int, 2> flags_sum_square = Flagger::ComputeAntennaFlags(
      sigma, max_iterations, stats_antenna_sum_square_);

  const xt::xtensor<int, 1> flagged_antennas =
      xt::flatten(flags_stddev) & xt::flatten(flags_sum_square);

  find_bad_antennas_timer_.stop();

  return flagged_antennas;
}

xt::xtensor<int, 1> Flagger::FindBadStations(float sigma, int max_iterations) {
  AssertStatsComputed();

  find_bad_stations_timer_.start();

  const xt::xtensor<int, 2> flags_stddev = Flagger::ComputeStationFlags(
      sigma, max_iterations, stats_antenna_stddev_);
  assert(flags_stddev.shape(0) == n_stations_);
  assert(flags_stddev.shape(1) == (n_correlations_ == 4 ? 2 : 1));

  const xt::xtensor<int, 2> flags_sum_square = Flagger::ComputeStationFlags(
      sigma, max_iterations, stats_antenna_sum_square_);

  const xt::xtensor<size_t, 1> station_indices = xt::flatten_indices(
      xt::where(xt::prod(flags_stddev & flags_sum_square, {1})));

  xt::xtensor<int, 1> antenna_flags({n_antennas_}, false);

  for (size_t station : station_indices) {
    const size_t first_antenna = station * n_antennas_per_station_;
    const size_t last_antenna = first_antenna + n_antennas_per_station_;
    xt::view(antenna_flags, xt::range(first_antenna, last_antenna)) = true;
  }

  find_bad_stations_timer_.stop();

  return antenna_flags;
}

}  // namespace antennaflagger
}  // namespace dp3
