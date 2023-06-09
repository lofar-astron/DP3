// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Flagger.h"

#include <aocommon/xt/span.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include "common/baseline_indices.h"

namespace {

const dp3::common::BaselineOrder kBaselineOrder =
    dp3::common::BaselineOrder::kRowMajor;

void CheckStatsStdDev(size_t n_baselines, size_t n_correlations,
                      const xt::xtensor<bool, 1>& baseline_flags,
                      xt::xtensor<std::complex<float>, 2> stats) {
  for (size_t bl = 0; bl < n_baselines; bl++) {
    auto baseline_stats = xt::view(stats, bl, xt::all());
    if (baseline_flags[bl]) {
      // The standard deviation should be zero for all flagged baselines.
      BOOST_CHECK(xt::allclose(baseline_stats, 0.0f));
    } else {
      BOOST_CHECK(xt::all(!xt::isclose(baseline_stats, 0.0f)));
    }
  }
}

void CheckStatsSumSquare(size_t n_baselines, size_t n_correlations,
                         const xt::xtensor<bool, 1>& baseline_flags,
                         xt::xtensor<std::complex<float>, 2> stats) {
  size_t first_flagged_baseline = 0;
  for (size_t bl = 0; bl < n_baselines; bl++) {
    if (baseline_flags[bl]) {
      first_flagged_baseline = bl;
      break;
    }
  }

  for (size_t bl = 0; bl < n_baselines; bl++) {
    auto baseline_stats = xt::view(stats, bl, xt::all());
    if (baseline_flags[bl]) {
      // The sum of squares should be the same for all flagged baselines.
      auto first_baseline_stats =
          xt::view(stats, first_flagged_baseline, xt::all());
      BOOST_CHECK(xt::allclose(baseline_stats, first_baseline_stats));
    } else {
      BOOST_CHECK(xt::all(!xt::isclose(baseline_stats, 0.0f)));
    }
  }
}

xt::xtensor<std::complex<float>, 3> CreateBrokenAntennaData(
    size_t n_baselines, size_t n_antennas, size_t n_channels,
    size_t n_correlations, size_t antenna1) {
  std::array<size_t, 3> data_shape{n_baselines, n_channels, n_correlations};
  xt::xtensor<std::complex<float>, 3> data(data_shape);

  // Generate input data
  xt::random::seed(0);
  xt::real(data) = xt::random::randn<float>(data.shape(), 0, 1);
  xt::imag(data) = xt::random::randn<float>(data.shape(), 0, 1);

  // Set broken antenna
  for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
    const size_t bl = dp3::common::ComputeBaselineIndex(
        antenna1, antenna2, n_antennas, kBaselineOrder);
    xt::view(data, bl, xt::all(), xt::all()) = 42.0f;
  }

  return data;
}

xt::xtensor<std::complex<float>, 3> CreateBrokenStationData(
    size_t n_baselines, size_t n_stations, size_t n_antennas_per_station,
    size_t n_channels, size_t n_correlations, size_t station1) {
  const size_t n_antennas = n_stations * n_antennas_per_station;

  std::array<size_t, 3> data_shape{n_baselines, n_channels, n_correlations};
  xt::xtensor<std::complex<float>, 3> data(data_shape);

  // Generate input data
  xt::random::seed(0);
  xt::real(data) = xt::random::randn<float>(data.shape(), 0, 1);
  xt::imag(data) = xt::random::randn<float>(data.shape(), 0, 1);

  // Set broken station
  for (size_t i = 0; i < n_antennas_per_station; ++i) {
    const size_t antenna1 = station1 * n_antennas_per_station + i;
    for (size_t antenna2 = 0; antenna2 < n_antennas; ++antenna2) {
      const size_t bl = dp3::common::ComputeBaselineIndex(
          antenna1, antenna2, n_antennas, kBaselineOrder);
      xt::view(data, bl, xt::all(), xt::all()) = 42.0f;
    }
  }

  return data;
}

xt::xtensor<bool, 1> ComputeBaselineSelectionForAntenna(size_t n_baselines,
                                                        size_t n_antennas,
                                                        size_t antenna) {
  xt::xtensor<bool, 1> selection({n_baselines}, false);
  const std::vector<size_t> indices =
      dp3::common::ComputeBaselineList(antenna, n_antennas, kBaselineOrder);
  xt::view(selection, xt::keep(indices)) = true;
  return selection;
}

xt::xtensor<bool, 1> ComputeBaselineSelectionForStation(
    size_t n_baselines, size_t n_stations, size_t n_antennas_per_station,
    size_t station) {
  const size_t n_antennas = n_stations * n_antennas_per_station;
  xt::xtensor<bool, 1> selection({n_baselines}, false);
  for (size_t i = 0; i < n_antennas_per_station; ++i) {
    const size_t antenna = station * n_antennas_per_station + i;
    const std::vector<size_t> indices =
        dp3::common::ComputeBaselineList(antenna, n_antennas, kBaselineOrder);
    xt::view(selection, xt::keep(indices)) = true;
  }
  return selection;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(antennaflagger)

BOOST_AUTO_TEST_CASE(compute_antenna_flags) {
  const unsigned int kNStations = 2;
  const unsigned int kNAntennasPerStation = 48;
  const unsigned int kNChannels = 8;
  const unsigned int kNCorrelations = 4;

  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  // Init
  const size_t antenna1 = 42;
  xt::xtensor<std::complex<float>, 3> data = CreateBrokenAntennaData(
      kNBaselines, kNAntennas, kNChannels, kNCorrelations, antenna1);
  const dp3::antennaflagger::Flagger::DataSpan data_span =
      aocommon::xt::CreateSpan(data);
  const xt::xtensor<bool, 1> baseline_flags =
      ComputeBaselineSelectionForAntenna(kNBaselines, kNAntennas, antenna1);

  // ComputeStats
  const xt::xtensor<std::complex<float>, 2> stats_baseline_stddev =
      dp3::antennaflagger::Flagger::ComputeStatsStdDev(data_span);
  CheckStatsStdDev(kNBaselines, kNCorrelations, baseline_flags,
                   stats_baseline_stddev);
  const xt::xtensor<std::complex<float>, 2> stats_baseline_sum_square =
      dp3::antennaflagger::Flagger::ComputeStatsSumSquare(data_span);
  CheckStatsSumSquare(kNBaselines, kNCorrelations, baseline_flags,
                      stats_baseline_sum_square);

  // GroupStats
  xt::xtensor<std::complex<float>, 3> stats_antenna_stddev =
      dp3::antennaflagger::Flagger::GroupStats(kNStations, kNAntennasPerStation,
                                               stats_baseline_stddev,
                                               kBaselineOrder);
  xt::xtensor<std::complex<float>, 3> stats_antenna_sum_square =
      dp3::antennaflagger::Flagger::GroupStats(kNStations, kNAntennasPerStation,
                                               stats_baseline_sum_square,
                                               kBaselineOrder);

  // ComputeAntennaFlags
  const float kSigma = 4;
  const int kMaxIterations = 2;
  const xt::xtensor<bool, 2> flags_stddev =
      dp3::antennaflagger::Flagger::ComputeAntennaFlags(kSigma, kMaxIterations,
                                                        stats_antenna_stddev);
  const xt::xtensor<bool, 2> flags_sum_square =
      dp3::antennaflagger::Flagger::ComputeAntennaFlags(
          kSigma, kMaxIterations, stats_antenna_sum_square);

  // Check flags
  xt::xtensor<bool, 2> expected_flags({kNStations, kNAntennasPerStation},
                                      false);
  expected_flags(0, antenna1) = true;
  BOOST_CHECK_EQUAL(flags_stddev, expected_flags);
  BOOST_CHECK_EQUAL(flags_sum_square, expected_flags);
}

BOOST_AUTO_TEST_CASE(find_bad_antennas) {
  const unsigned int kNStations = 2;
  const unsigned int kNAntennasPerStation = 48;
  const unsigned int kNChannels = 8;
  const unsigned int kNCorrelations = 4;

  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  // Init
  const size_t antenna1 = 42;
  xt::xtensor<std::complex<float>, 3> data = CreateBrokenAntennaData(
      kNBaselines, kNAntennas, kNChannels, kNCorrelations, antenna1);
  const dp3::antennaflagger::Flagger::DataSpan data_span =
      aocommon::xt::CreateSpan(data);

  // Flagger
  dp3::antennaflagger::Flagger flagger(kNStations, kNAntennasPerStation,
                                       kNChannels, kNCorrelations);

  flagger.ComputeStats(data_span, kBaselineOrder);

  const float kSigma = 4;
  const int kMaxIterations = 2;

  const xt::xtensor<bool, 1> antenna_flags =
      flagger.FindBadAntennas(kSigma, kMaxIterations);

  // Check flags
  BOOST_REQUIRE_EQUAL(antenna_flags.size(), kNAntennas);
  xt::xtensor<bool, 1> expected_flags({kNAntennas}, false);
  expected_flags(antenna1) = true;
  BOOST_CHECK_EQUAL(antenna_flags, expected_flags);
}

BOOST_AUTO_TEST_CASE(compute_station_flags) {
  const unsigned int kNStations = 12;
  const unsigned int kNAntennasPerStation = 1;
  const unsigned int kNChannels = 8;
  const unsigned int kNCorrelations = 4;

  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  // Init
  const size_t station1 = 6;
  xt::xtensor<std::complex<float>, 3> data =
      CreateBrokenStationData(kNBaselines, kNStations, kNAntennasPerStation,
                              kNChannels, kNCorrelations, station1);
  const dp3::antennaflagger::Flagger::DataSpan data_span =
      aocommon::xt::CreateSpan(data);

  const xt::xtensor<bool, 1> baseline_flags =
      ComputeBaselineSelectionForStation(kNBaselines, kNStations,
                                         kNAntennasPerStation, station1);

  // ComputeStats
  const xt::xtensor<std::complex<float>, 2> stats_baseline_stddev =
      dp3::antennaflagger::Flagger::ComputeStatsStdDev(data_span);
  CheckStatsStdDev(kNBaselines, kNCorrelations, baseline_flags,
                   stats_baseline_stddev);
  const xt::xtensor<std::complex<float>, 2> stats_baseline_sum_square =
      dp3::antennaflagger::Flagger::ComputeStatsSumSquare(data_span);
  CheckStatsSumSquare(kNBaselines, kNCorrelations, baseline_flags,
                      stats_baseline_sum_square);

  // GroupStats
  const xt::xtensor<std::complex<float>, 3> stats_antenna_stddev =
      dp3::antennaflagger::Flagger::GroupStats(kNStations, kNAntennasPerStation,
                                               stats_baseline_stddev,
                                               kBaselineOrder);
  const xt::xtensor<std::complex<float>, 3> stats_antenna_sum_square =
      dp3::antennaflagger::Flagger::GroupStats(kNStations, kNAntennasPerStation,
                                               stats_baseline_sum_square,
                                               kBaselineOrder);

  // ComputeStationFlags
  const float kSigma = 3;
  const int kMaxIterations = 2;

  const xt::xtensor<bool, 2> flags_stddev =
      dp3::antennaflagger::Flagger::ComputeStationFlags(kSigma, kMaxIterations,
                                                        stats_antenna_stddev);
  BOOST_REQUIRE_EQUAL(flags_stddev.shape(0), kNStations);
  BOOST_REQUIRE_EQUAL(flags_stddev.shape(1), 2);  // XX and XY only

  const xt::xtensor<bool, 2> flags_sum_square =
      dp3::antennaflagger::Flagger::ComputeStationFlags(
          kSigma, kMaxIterations, stats_antenna_sum_square);
  BOOST_REQUIRE_EQUAL(flags_sum_square.shape(0), kNStations);
  BOOST_REQUIRE_EQUAL(flags_sum_square.shape(1), 2);  // XX and XY only

  // Check flags
  xt::xtensor<bool, 2> expected_flags({kNStations, 2}, false);
  xt::view(expected_flags, station1, xt::all()) = true;
  BOOST_CHECK_EQUAL(flags_stddev, expected_flags);
  BOOST_CHECK_EQUAL(flags_sum_square, expected_flags);
}

BOOST_AUTO_TEST_CASE(find_bad_stations) {
  const unsigned int kNStations = 12;
  const unsigned int kNAntennasPerStation = 48;
  const unsigned int kNChannels = 8;
  const unsigned int kNCorrelations = 4;

  const unsigned int kNAntennas = kNStations * kNAntennasPerStation;
  const unsigned int kNBaselines = dp3::common::ComputeNBaselines(kNAntennas);

  // Init
  const size_t station = 6;
  xt::xtensor<std::complex<float>, 3> data =
      CreateBrokenStationData(kNBaselines, kNStations, kNAntennasPerStation,
                              kNChannels, kNCorrelations, station);
  const dp3::antennaflagger::Flagger::DataSpan data_span =
      aocommon::xt::CreateSpan(data);

  // Flagger
  dp3::antennaflagger::Flagger flagger(kNStations, kNAntennasPerStation,
                                       kNChannels, kNCorrelations);

  flagger.ComputeStats(data_span, kBaselineOrder);

  const float kSigma = 3;
  const int kMaxIterations = 2;

  const xt::xtensor<bool, 1> antenna_flags =
      flagger.FindBadStations(kSigma, kMaxIterations);

  // Check flags
  BOOST_REQUIRE_EQUAL(xt::sum(antenna_flags)(), kNAntennasPerStation);
  xt::xtensor<bool, 1> expected_flags({kNAntennas}, false);
  const size_t first_antenna = station * kNAntennasPerStation;
  const size_t last_antenna = first_antenna + kNAntennasPerStation;
  xt::view(expected_flags, xt::range<size_t>(first_antenna, last_antenna)) =
      true;
  BOOST_CHECK_EQUAL(antenna_flags, expected_flags);
}

BOOST_AUTO_TEST_SUITE_END()
