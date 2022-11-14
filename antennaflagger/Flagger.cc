// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <set>

#include "Flagger.h"
#include "./statistics.h"

namespace detail {
inline size_t DataIndex(size_t n_channels, size_t n_correlations, size_t bl,
                        size_t chan, size_t cor) {
  // data[n_baselines][n_channels][n_correlations]
  return (bl * n_channels + chan) * n_correlations + cor;
}

/*
 * Helper functions for AartfaacFlagger::find_bad_antennas
 */
template <typename InType, typename OutType>
void ComputeStatsPerBaseline(size_t n_baselines, size_t n_channels,
                             size_t n_correlations,
                             const std::vector<std::complex<InType>>& data,
                             std::vector<std::complex<OutType>>* stats_std,
                             std::vector<std::complex<OutType>>* stats_sump2) {
  //  -  in: data[n_baselines][n_channels][n_correlations]
  //  - out: stats[n_correlations][n_baselines]
  // Note that the baseline and correlation indices are swapped between
  // input and output. This is also different from the Python reference,
  // and is to avoid strided access later on where the stats are used.

  [[maybe_unused]] size_t n_in = n_baselines * n_channels * n_correlations;
  assert(data.size() == n_in);

  std::vector<OutType> result(n_baselines * n_correlations);
  size_t n_out = n_correlations * n_baselines;
  stats_std->resize(n_out);
  stats_sump2->resize(n_out);

  for (size_t bl = 0; bl < n_baselines; ++bl) {
    for (size_t cor = 0; cor < n_correlations; ++cor) {
      // Compute sum, sump2
      OutType sum_real = 0;
      OutType sum_imag = 0;
      OutType sump2_real = 0;
      OutType sump2_imag = 0;
      for (size_t chan = 0; chan < n_channels; ++chan) {
        size_t idx = DataIndex(n_channels, n_correlations, bl, chan, cor);
        std::complex<InType> value = data[idx];
        sum_real += value.real();
        sum_imag += value.imag();
        sump2_real += value.real() * value.real();
        sump2_imag += value.imag() * value.imag();
      }

      // Compute mean
      OutType mean_real = sum_real / n_channels;
      OutType mean_imag = sum_imag / n_channels;

      // Compute standard deviation
      OutType sum2_real = 0;
      OutType sum2_imag = 0;
      for (size_t chan = 0; chan < n_channels; ++chan) {
        size_t idx = DataIndex(n_channels, n_correlations, bl, chan, cor);
        std::complex<InType> value = data[idx];
        OutType diff_real = mean_real - value.real();
        OutType diff_imag = mean_imag - value.imag();
        sum2_real += diff_real * diff_real;
        sum2_imag += diff_imag * diff_imag;
      }
      OutType std_real = std::sqrt(sum2_real / n_channels);
      OutType std_imag = std::sqrt(sum2_imag / n_channels);

      // Store result
      size_t idx = cor * n_baselines + bl;
      (*stats_std)[idx] = std::complex<OutType>(std_real, std_imag);
      (*stats_sump2)[idx] = std::complex<OutType>(sump2_real, sump2_imag);
    }
  }
}  // end ComputeStatsPerBaseline

template <typename InType, typename OutType>
void ComputeStatsPerAntenna(
    size_t n_receivers, size_t n_baselines, size_t n_channels,
    size_t n_correlations, const std::vector<std::complex<InType>>& data,
    std::vector<std::complex<OutType>>* stats_std_grouped,
    std::vector<std::complex<OutType>>* stats_sump2_grouped) {
  // Compute stats for every baseline and correlation:
  //  -  in: data[n_baselines][n_channels][n_correlations]
  //  - out: stats_std[n_correlations][n_baselines]
  //  - out: stats_sump2[n_correlations][n_baselines]
  assert(data.size() == size_t(n_baselines * n_channels * n_correlations));
  std::vector<std::complex<OutType>> stats_std;
  std::vector<std::complex<OutType>> stats_sump2;
  ComputeStatsPerBaseline(n_baselines, n_channels, n_correlations, data,
                          &stats_std, &stats_sump2);
  assert(stats_std.size() == size_t(n_correlations * n_baselines));
  assert(stats_sump2.size() == size_t(n_correlations * n_baselines));

  // Group stats by antenna:
  // -  in: stats_std[n_correlations][n_baselines]
  // -  in: stats_sump2[n_correlations][n_baselines]
  // - out: stats_std_grouped[n_correlations][n_receivers]
  // - out: stats_sump2_grouped[n_correlations][n_receivers]
  stats_std_grouped->resize(n_receivers * n_correlations);
  stats_sump2_grouped->resize(n_receivers * n_correlations);

  for (size_t antenna = 0; antenna < n_receivers; ++antenna) {
    // Find all baselines with this antenna
    std::vector<size_t> indices = dp3::common::ComputeBaselineList(
        antenna, n_receivers, dp3::common::BaselineOrder::kRowMajor);

    // Sum statistics over these baselines
    for (size_t cor = 0; cor < n_correlations; ++cor) {
      for (const size_t& bl : indices) {
        size_t idx_src = cor * n_baselines + bl;
        size_t idx_dst = cor * n_receivers + antenna;
        (*stats_std_grouped)[idx_dst] += stats_std[idx_src];
        (*stats_sump2_grouped)[idx_dst] += stats_sump2[idx_src];
      }
    }
  }
}  // end ComputeStatsPerAntenna

/*
 * Helper functions for AaatfaacFlagger::find_bad_stations
 */
template <typename T>
void FindBadAntennas(size_t n_stations, size_t n_receivers_per_station,
                     size_t n_baselines, size_t n_correlations, float sigma,
                     int maxiters, const std::vector<std::complex<T>>& stats,
                     std::vector<size_t>* all_bad_stations) {
  int n_receivers = n_stations * n_receivers_per_station;
  // -  in: stats[n_correlations][n_receivers]
  // - out: all_bad_stations[?]
  assert(n_correlations == static_cast<int>(4));
  assert(stats.size() == static_cast<size_t>(n_correlations * n_receivers));

  // Find bad stations by looking at the statistics (per correlation) for all
  // antennas that belong to one station. The results are stored in a list of
  // lists such that stations can be processed independently.
  std::vector<std::vector<int>> bad_stations(n_stations);

  // Check XY, YX and YY
  // TODO: find out why XX is not checked and/or make this configurable
  for (size_t cor = 1; cor < n_correlations; ++cor) {
    // The check sum is the list of statistics for the current correlation
    //  - in: T2 check_sum[n_correlations][n_receivers]
    const std::complex<T>* check_sum = &stats[cor * n_receivers];

    for (size_t station = 0; station < n_stations; ++station) {
      // Make a list of all antennas for the current station,
      // and select the corresponding values from check_sum.
      std::vector<int> station_indices(n_receivers_per_station);
      std::vector<T> check_sum_station(n_receivers_per_station * 2);
      for (size_t i = 0; i < n_receivers_per_station; ++i) {
        const size_t station_index = station * n_receivers_per_station + i;
        station_indices[i] = station_index;
        T real =
            std::log(static_cast<T>(1.0) + check_sum[station_index].real());
        T imag =
            std::log(static_cast<T>(1.0) + check_sum[station_index].imag());
        check_sum_station[0 * n_receivers_per_station + i] = real;
        check_sum_station[1 * n_receivers_per_station + i] = imag;
      }

      // Find outliers
      const size_t n_flags = n_receivers_per_station * 2;
      std::vector<bool> flags =
          FindOutliers(sigma, maxiters, check_sum_station);
      assert(flags.size() == size_t(n_flags));

      // Count number of outliers
      int n_outliers = 0;
      for (size_t i = 0; i < n_flags; ++i) {
        n_outliers += flags[i];
      }

      // Add the outliers to the list of bad stations
      const size_t n = bad_stations[station].size();
      bad_stations[station].resize(n + n_outliers);
      size_t j = 0;
      for (size_t i = 0; i < n_receivers_per_station * 2; ++i) {
        if (flags[i]) {
          size_t station_index =
              station * n_receivers_per_station + (i % n_receivers_per_station);
          bad_stations[station][n + j++] = station_index;
        }
      }
    }
  }

  // Add the bad stations to the global list of bad stations
  for (size_t station = 0; station < n_stations; ++station) {
    all_bad_stations->insert(all_bad_stations->end(),
                             bad_stations[station].begin(),
                             bad_stations[station].end());
  }
}  // end FindBadAntennas

template <typename T>
std::vector<std::complex<T>> ComputeMeans(size_t n_stations,
                                          size_t n_receivers_per_station,
                                          const std::complex<T>* check_sum) {
  const size_t n_receivers = n_stations * n_receivers_per_station;

  // -  in: T2 check_sum[n_receivers]
  // - out: T2 means[n_stations]
  std::vector<std::complex<T>> means(n_stations);

  // Split the real and imaginary parts of the check sum
  std::vector<T> check_sum_real(n_receivers);
  std::vector<T> check_sum_imag(n_receivers);
  for (size_t antenna = 0; antenna < n_receivers; ++antenna) {
    check_sum_real[antenna] = check_sum[antenna].real();
    check_sum_imag[antenna] = check_sum[antenna].imag();
  }

  // Compute the median for both real and imaginary parts
  T median_real = ComputeMedian(check_sum_real);
  T median_imag = ComputeMedian(check_sum_imag);

  for (size_t station = 0; station < n_stations; ++station) {
    // Make a list of all antennas for the current station,
    // and select the corresponding values from check_sum
    std::vector<T> check_sum_station_real(n_receivers_per_station);
    std::vector<T> check_sum_station_imag(n_receivers_per_station);
    for (size_t i = 0; i < n_receivers_per_station; ++i) {
      size_t station_index = station * n_receivers_per_station + i;
      check_sum_station_real[i] = check_sum_real[station_index];
      check_sum_station_imag[i] = check_sum_imag[station_index];
    }

    // Compute the median for both the real and imaginary parts
    // of the check sum for the current station
    T median_station_real = ComputeMedian(check_sum_station_real);
    T median_station_imag = ComputeMedian(check_sum_station_imag);

    // Add the mean value for the current station
    means[station] = {median_station_real / median_real,
                      median_station_imag / median_imag};
  }

  return means;
}  // end ComputeMeans

};  // namespace detail

/*
 * Member functions
 */
namespace dp3 {
namespace antennaflagger {

void Flagger::ReportRuntime() const {
  std::cout << "statistics: " << computeStatisticsWatch_.getElapsed()
            << ", bad_antennas: " << findbadAntennasWatch_.getElapsed()
            << ", bad_stations: " << findBadStationsWatch_.getElapsed() << '\n';
}

void Flagger::ComputeStats(const std::vector<std::complex<double>>& data) {
  computeStatisticsWatch_.start();
  detail::ComputeStatsPerAntenna(n_receivers_, n_baselines_, n_channels_,
                                 n_correlations_, data, &stats_std_,
                                 &stats_sump2_);
  assert(stats_std_.size() == size_t(n_correlations_ * n_receivers_));
  assert(stats_sump2_.size() == size_t(n_correlations_ * n_receivers_));
  computeStatisticsWatch_.stop();
}

void Flagger::ComputeStats(const std::vector<std::complex<float>>& data) {
  computeStatisticsWatch_.start();
  detail::ComputeStatsPerAntenna(n_receivers_, n_baselines_, n_channels_,
                                 n_correlations_, data, &stats_std_,
                                 &stats_sump2_);
  assert(stats_std_.size() == size_t(n_correlations_ * n_receivers_));
  assert(stats_sump2_.size() == size_t(n_correlations_ * n_receivers_));
  computeStatisticsWatch_.stop();
}

void Flagger::AssertStatsComputed() const {
  if (!(stats_std_.size() == (n_correlations_ * n_receivers_)) ||
      !(stats_sump2_.size() == (n_correlations_ * n_receivers_))) {
    throw std::runtime_error("Call Flagger::compute_stats() first!");
  }
}

std::vector<size_t> Flagger::FindBadAntennas(float sigma, int maxiters) {
  AssertStatsComputed();

  findbadAntennasWatch_.start();

  std::vector<size_t> flagged_antennas;
  detail::FindBadAntennas(n_stations_, n_receivers_per_station, n_baselines_,
                          n_correlations_, sigma, maxiters, stats_sump2_,
                          &flagged_antennas);
  detail::FindBadAntennas(n_stations_, n_receivers_per_station, n_baselines_,
                          n_correlations_, sigma, maxiters, stats_std_,
                          &flagged_antennas);

  // Count the number of times an antenna is flagged
  std::vector<size_t> flagged_count(n_receivers_);
  for (size_t antenna : flagged_antennas) {
    flagged_count[antenna]++;
  }

  // Make a list of antennas that are flagged using both of the statistics that
  // are used above.
  std::vector<size_t> flagged_antennas2;
  for (size_t antenna = 0; antenna < n_receivers_; ++antenna) {
    if (flagged_count[antenna] > size_t(2)) {
      flagged_antennas2.push_back(antenna);
    }
  }

  findbadAntennasWatch_.stop();

  return flagged_antennas2;
}

std::vector<size_t> Flagger::FindBadStations(float sigma, int maxiters) {
  AssertStatsComputed();

  findBadStationsWatch_.start();

  // Compute means for all the check sums seperately for
  // the real and imaginary part of the check sum.
  // -  in: StatsType2 check_sums[n_check_sums][n_receivers]
  // - out: StatsType all_means[check_sums.size() * 2][n_stations]
  std::array<const std::complex<StatsType>*, 4> check_sums{
      &stats_std_[0 * n_receivers_], &stats_std_[3 * n_receivers_],
      &stats_sump2_[0 * n_receivers_], &stats_sump2_[3 * n_receivers_]};
  std::vector<std::vector<StatsType>> all_means(check_sums.size() * 2);

  for (std::vector<StatsType>& means : all_means) {
    means.resize(n_stations_);
  }

  for (size_t i = 0; i < check_sums.size(); ++i) {
    std::vector<std::complex<StatsType>> means = detail::ComputeMeans(
        n_stations_, n_receivers_per_station, check_sums[i]);
    for (size_t station = 0; station < n_stations_; ++station) {
      all_means[0 * check_sums.size() + i][station] = means[station].real();
      all_means[1 * check_sums.size() + i][station] = means[station].imag();
    }
  }

  // Flag all stations for which all of the means are an outlier
  std::vector<bool> station_flags(n_stations_, {false});
  for (size_t i = 0; i < all_means.size(); ++i) {
    const std::vector<StatsType>& means = all_means[i];

    // Find outliers
    std::vector<bool> flags = FindOutliers(sigma, maxiters, means);
    assert(flags.size() == size_t(n_stations_));

    for (size_t station = 0; station < n_stations_; ++station) {
      station_flags[station] =
          (flags[station] && (i == 0 || station_flags[station]));
    }
  }

  // Make a list of all the antenna indices for the flagged stations
  std::vector<size_t> flagged_antennas;
  for (size_t station = 0; station < n_stations_; ++station) {
    if (station_flags[station]) {
      // Make space for the new indices
      const size_t n = flagged_antennas.size();
      flagged_antennas.resize(n + n_receivers_per_station);

      // Add the new indices
      for (size_t i = 0; i < n_receivers_per_station; ++i) {
        const size_t station_index = station * n_receivers_per_station + i;
        flagged_antennas[n + i] = station_index;
      }
    }
  }

  findBadStationsWatch_.stop();

  return flagged_antennas;
}

}  // namespace antennaflagger
}  // namespace dp3
