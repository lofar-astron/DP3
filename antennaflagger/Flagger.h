// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_ANTENNAFLAGGER_H
#define DP3_ANTENNAFLAGGER_H

#include <vector>
#include <complex>

#include "../common/baseline_indices.h"
#include "../common/Timer.h"

namespace dp3 {
namespace antennaflagger {

/*
 * AntennaFlagger uses statistics to find bad antennas for a visibility data
 * set. It is based on the following Python reference code:
 * https://github.com/transientskp/A12_pipeline/blob/master/pyscripts/calc_antflags.py
 * and comprises four main steps:
 *  1) Compute standard deviation and sump2 for every baseline and correlation
 *     seperately. These statistics are integrated over frequency (channels)
 *  2) Group the statistics per antenna by integrating over all baselines where
 *     a particular antenna is involved. These statistics are stored as class
 *     members.
 *  3) Query these statistics to find outlier antennas
 *  4) Query these statistics to find outlier stations
 */
class Flagger {
 public:
  Flagger(size_t n_stations, size_t n_receivers_per_station, size_t n_channels,
          size_t n_correlations)
      : n_stations_(n_stations),
        n_receivers_per_station(n_receivers_per_station),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        n_receivers_(n_stations * n_receivers_per_station),
        n_baselines_(dp3::common::ComputeNBaselines(n_receivers_)) {}

  void ComputeStats(const std::vector<std::complex<double>>& data);
  void ComputeStats(const std::vector<std::complex<float>>& data);

  std::vector<size_t> FindBadAntennas(float sigma, int maxiters);
  std::vector<size_t> FindBadStations(float sigma, int maxiters);

  void ReportRuntime() const;

 private:
  // Helper function
  void AssertStatsComputed() const;

  // Constants
  const size_t n_stations_;
  const size_t n_receivers_per_station;
  const size_t n_channels_;
  const size_t n_correlations_;

  // Derived constants
  const size_t n_receivers_;
  const size_t n_baselines_;

  // Timings
  common::NSTimer computeStatisticsWatch_;
  common::NSTimer findbadAntennasWatch_;
  common::NSTimer findBadStationsWatch_;

  // Data
  using StatsType = float;
  std::vector<std::complex<StatsType>> stats_std_;
  std::vector<std::complex<StatsType>> stats_sump2_;
};

}  // namespace antennaflagger
}  // namespace dp3

#endif  // DP3_ANTENNAFLAGGER_H
