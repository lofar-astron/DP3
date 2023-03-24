// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_ANTENNAFLAGGER_FLAGGER_H_
#define DP3_ANTENNAFLAGGER_FLAGGER_H_

#include <vector>
#include <complex>

#include <aocommon/xt/span.h>

#include "../common/baseline_indices.h"
#include "../common/Timer.h"

namespace dp3 {
namespace antennaflagger {

/*
 * AntennaFlagger uses statistics to find bad antennas for a visibility data
 * set. It is based on the following Python reference code:
 * https://github.com/transientskp/A12_pipeline/blob/master/pyscripts/calc_antflags.py
 * and comprises four main steps:
 *  1) Compute standard deviation and sum of squares for every baseline and
 *     correlation seperately. These statistics are integrated over frequency
 *     (channels).
 *  2) Group the statistics per antenna by integrating over all baselines where
 *     a particular antenna is involved. These statistics are stored as class
 *     members.
 *  3) Query these statistics to find outlier antennas.
 *  4) Query these statistics to find outlier stations.
 */
class Flagger {
 public:
  using DataSpan = aocommon::xt::Span<std::complex<float>, 3>;

  Flagger(size_t n_stations, size_t n_antennas_per_station, size_t n_channels,
          size_t n_correlations)
      : n_stations_(n_stations),
        n_antennas_per_station_(n_antennas_per_station),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        n_antennas_(n_stations * n_antennas_per_station),
        n_baselines_(common::ComputeNBaselines(n_antennas_)) {}

  /**
   * This helper function computes statistics over the input visibility data.
   * These statistics are used by FindBadAntennas and FindBadStations, which
   * perform the actual flagging. Therefore, this function should always be
   * called first.
   *
   * @param data Three-dimensional input data of the form (baseline,
   * channels, polarizations).
   * @param baseline_order Baselines are assumed to be ordered according to this
   * ordering. See common/baseline_indices.h for more details on the ordering.
   */
  void ComputeStats(const DataSpan& data, common::BaselineOrder baseline_order);

  /**
   * Identify antennas that are outliers compared to the other antennas
   * belonging to the same station. Note that this is only usefull for
   * telescopes such as AARTFAAC where n_antennas_per_station > 1.
   *
   * @return Flags as a tensor of booleans of length n_stations *
   * n_antennas_per_station.
   */
  xt::xtensor<bool, 1> FindBadAntennas(float sigma, int max_iterations);

  /**
   * Identify stations that are outliers compared to the other stations.
   * In case of n_antennas_per_station > 1, a per-station statistic over all
   * antennas in a station is used.
   *
   * @return Flags as a tensor of booleans of length n_stations *
   * n_antennas_per_station.
   */
  xt::xtensor<bool, 1> FindBadStations(float sigma, int max_iterations);

  /**
   * Compute standard deviation for the real and imaginary component seperately
   * over the second axis of the input.
   *
   * @return Statistics as a tensor of complex numbers.
   */
  static xt::xtensor<std::complex<float>, 2> ComputeStatsStdDev(
      const DataSpan& data);

  /**
   * Compute sum of squares for the real and imaginary component seperately over
   * the second axis of the input.
   *
   * @return Statistics as a tensor of complex numbers.
   */
  static xt::xtensor<std::complex<float>, 2> ComputeStatsSumSquare(
      const DataSpan& data);

  /**
   * Compute the per antenna statistics by combining the
   * statistics of all baselines an antenna contributed to. The input is in the
   * form (n_stations, n_antennas_per_station, n_correlations),
   *
   * @return Statistics in the form (n_stations, n_antennas_per_station,
   * n_correlations).
   */
  static xt::xtensor<std::complex<float>, 3> GroupStats(
      size_t n_stations, size_t n_antennas_per_station,
      const xt::xtensor<std::complex<float>, 2>& stats_baseline,
      common::BaselineOrder baseline_order);

  /**
   * Per station, look at the statistics for all its antennas to find outliers.
   * The input is assumed to be in the form (n_stations, n_antennas_per_station,
   * n_correlations).
   *
   * @return Flags as a tensor of booleans in the form (n_stations,
   * n_antennas_per_station).
   */
  static xt::xtensor<bool, 2> ComputeAntennaFlags(
      float sigma, int max_iterations,
      const xt::xtensor<std::complex<float>, 3>& stats);

  /**
   * Per correlation (XX in case of n_correlations == 1, XX and YY in case of
   * n_correlations == 4), take the median of the statistics for all antennas in
   * a station and compare these to find outliers. The input is assumed to be in
   * the form (n_stations, n_antennas_per_station, n_correlations).
   *
   * @return Flags as a tensor of booleans in the form (n_stations, 1 || 2).
   */
  static xt::xtensor<bool, 2> ComputeStationFlags(
      float sigma, int max_iterations,
      const xt::xtensor<std::complex<float>, 3>& stats);

 private:
  // Helper function
  void AssertStatsComputed() const;

  // Constants
  const size_t n_stations_;
  const size_t n_antennas_per_station_;
  const size_t n_channels_;
  const size_t n_correlations_;

  // Derived constants
  const size_t n_antennas_;
  const size_t n_baselines_;

  // Timings
  common::NSTimer compute_statistics_timer_;
  common::NSTimer find_bad_antennas_timer_;
  common::NSTimer find_bad_stations_timer_;

  // Data
  xt::xtensor<std::complex<float>, 3> stats_antenna_stddev_;
  xt::xtensor<std::complex<float>, 3> stats_antenna_sum_square_;
};

}  // namespace antennaflagger
}  // namespace dp3

#endif  // DP3_ANTENNAFLAGGER_FLAGGER_H_
