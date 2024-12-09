// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/SolveData.h"

#include <dp3/base/BdaBuffer.h>
#include <dp3/base/DPBuffer.h>
#include "../../gain_solvers/BdaSolverBuffer.h"
#include "../../gain_solvers/SolverTools.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iterator>
#include <random>

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using ChannelBlockData = dp3::ddecal::FullSolveData::ChannelBlockData;

namespace {

const size_t kNPolarizations = 4;
const size_t kNChannels = 7;
const size_t kNChannelBlocks = 2;
const size_t kNBaselines = 3;
const std::array<size_t, 3> kShape{kNBaselines, kNChannels, kNPolarizations};
const size_t kNAntennas = 3;
const std::vector<int> kAntennas1 = {0, 0, 0};
const std::vector<int> kAntennas2 = {1, 2, 0};

void FillRandomData(std::vector<std::complex<float>>& data, size_t size) {
  // Reusing these variables ensures that each call generates different data.
  static std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  static std::mt19937 mt(0);

  data.clear();
  data.reserve(size);
  std::generate_n(std::back_inserter(data), size, []() {
    return std::complex<float>(uniform_data(mt), uniform_data(mt));
  });
}

void FillRegularData(DPBuffer::DataType& data) {
  // Reusing these variables ensures that each call generates different data.
  static std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  static std::mt19937 mt(0);

  std::generate_n(data.begin(), data.size(), []() {
    return std::complex<float>(uniform_data(mt), uniform_data(mt));
  });
}

void FillBdaBuffer(BdaBuffer& buffer, size_t avg_channels, size_t all_channels,
                   double unit_interval) {
  std::vector<std::complex<float>> data;
  const std::vector<float> weights(all_channels * kNPolarizations, 1.0f);

  for (size_t time_index = 0; time_index != 2; ++time_index) {
    const double row_time = 2.0 * time_index * unit_interval;
    // Add averaged rows for baselines 0 (cross-correlation 0 x 1)
    FillRandomData(data, avg_channels * kNPolarizations);
    BOOST_REQUIRE(buffer.AddRow(
        row_time + unit_interval, unit_interval * 2.0, unit_interval * 2.0, 0,
        avg_channels, kNPolarizations, data.data(), nullptr, weights.data()));

    // Add averaged rows for baselines 2 (auto-correlation 0 x 0)
    FillRandomData(data, avg_channels * kNPolarizations);
    BOOST_REQUIRE(buffer.AddRow(
        row_time + unit_interval, unit_interval * 2.0, unit_interval * 2.0, 2,
        avg_channels, kNPolarizations, data.data(), nullptr, weights.data()));

    // Add non-averaged rows for baseline 1 (cross-correlation 0 x 2)
    for (int j = 0; j < 2; ++j) {
      FillRandomData(data, all_channels * kNPolarizations);
      BOOST_REQUIRE(buffer.AddRow(
          row_time + (0.5 + j) * unit_interval, unit_interval, unit_interval, 1,
          all_channels, kNPolarizations, data.data(), nullptr, weights.data()));
    }
  }
}

}  // namespace

// SolveData should only copy data, so the tolerance is zero for these tests.
BOOST_AUTO_TEST_SUITE(solve_data, *boost::unit_test::tolerance(0.0f))

BOOST_AUTO_TEST_CASE(regular) {
  // Test that all data from several DPBuffers ends up in a SolveData structure.
  // The test uses three baselines. The third baseline contains
  // auto-correlations, which SolveData should ignore.
  const size_t kNTimes = 2;
  const size_t kNSolveDataBaselines = kNBaselines - 1;
  const size_t kNDirections = 1;
  const std::string kDirectionName = "direction_foo";
  const std::vector<size_t> kNSolutionsPerDirection(kNDirections, 1);
  const std::vector<size_t> kExpectedChannelBlockSizes{3, 4};

  std::vector<std::unique_ptr<DPBuffer>> unweighted_buffers;
  std::vector<DPBuffer> weighted_buffers(kNTimes);

  for (size_t time = 0; time < kNTimes; ++time) {
    unweighted_buffers.emplace_back(std::make_unique<DPBuffer>(time, 1.0));
    unweighted_buffers.back()->GetData().resize(kShape);
    FillRegularData(unweighted_buffers.back()->GetData(""));

    unweighted_buffers.back()->AddData(kDirectionName);
    FillRegularData(unweighted_buffers.back()->GetData(kDirectionName));

    unweighted_buffers.back()->GetWeights().resize(kShape);
    unweighted_buffers.back()->GetWeights().fill(1.0f);

    unweighted_buffers.back()->GetFlags().resize(kShape);
    unweighted_buffers.back()->GetFlags().fill(false);
  }

  dp3::ddecal::AssignAndWeight(unweighted_buffers, {kDirectionName},
                               weighted_buffers, false, false);

  const dp3::ddecal::SolveData data(
      weighted_buffers, {kDirectionName}, kNChannelBlocks, kNAntennas,
      kNSolutionsPerDirection, kAntennas1, kAntennas2);
  BOOST_TEST_REQUIRE(data.NChannelBlocks() == kNChannelBlocks);

  for (size_t ch_block = 0; ch_block < kNChannelBlocks; ++ch_block) {
    const ChannelBlockData& cb_data = data.ChannelBlock(ch_block);
    BOOST_TEST_REQUIRE(cb_data.NDirections() == kNDirections);
    // SolveData does not include baselines with auto-correlations,
    // so it should not contain the last baseline.
    BOOST_TEST_REQUIRE(cb_data.NVisibilities() ==
                       kNTimes * kNSolveDataBaselines *
                           kExpectedChannelBlockSizes[ch_block]);
    for (size_t ant = 0; ant < kNAntennas; ++ant) {
      BOOST_TEST(cb_data.NAntennaVisibilities(ant) ==
                 kNTimes * ((ant == 0) ? 2 : 1) *
                     kExpectedChannelBlockSizes[ch_block]);
    }

    const size_t first_channel = ch_block * kNChannels / kNChannelBlocks;
    for (size_t v = 0; v < cb_data.NVisibilities(); ++v) {
      const size_t time =
          v / (kNSolveDataBaselines * kExpectedChannelBlockSizes[ch_block]);
      const size_t baseline =
          (v / kExpectedChannelBlockSizes[ch_block]) % kNSolveDataBaselines;
      const size_t channel =
          first_channel + v % kExpectedChannelBlockSizes[ch_block];

      BOOST_TEST(cb_data.Antenna1Index(v) == size_t(kAntennas1[baseline]));
      BOOST_TEST(cb_data.Antenna2Index(v) == size_t(kAntennas2[baseline]));

      const std::complex<float>* const expected_data =
          &weighted_buffers[time].GetData("")(baseline, channel, 0);
      const aocommon::MC2x2F& data = cb_data.Visibility(v);
      for (size_t pol = 0; pol < kNPolarizations; ++pol) {
        BOOST_TEST(expected_data[pol] == data.Get(pol));
      }

      for (size_t direction = 0; direction < kNDirections; ++direction) {
        const std::complex<float>* const expected_model_data =
            &weighted_buffers[time].GetData(kDirectionName)(baseline, channel,
                                                            0);
        const aocommon::MC2x2F& model_data =
            cb_data.ModelVisibility(direction, v);

        BOOST_TEST(cb_data.SolutionIndex(direction, v) == direction);
        BOOST_TEST(cb_data.NSolutionVisibilities(direction, 0) ==
                   kNTimes * kNSolveDataBaselines *
                       kExpectedChannelBlockSizes[ch_block]);

        for (size_t pol = 0; pol < kNPolarizations; ++pol) {
          BOOST_TEST(expected_model_data[pol] == model_data.Get(pol));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(regular_with_dd_intervals) {
  // Test that all data from several DPBuffers ends up in a SolveData structure.
  // The test uses three baselines. The third baseline contains
  // auto-correlations, which SolveData should ignore.
  const size_t kNTimes = 2;
  const size_t kNSolveDataBaselines = kNBaselines - 1;
  const std::vector<std::string> kDirectionNames{"foo_direction", "bar_dir"};
  const std::vector<size_t> kNSolutionsPerDirection{1, 2};
  const std::vector<size_t> kExpectedChannelBlockSizes{3, 4};

  std::vector<std::unique_ptr<DPBuffer>> unweighted_buffers;
  std::vector<DPBuffer> weighted_buffers(kNTimes);

  for (size_t time = 0; time < kNTimes; ++time) {
    unweighted_buffers.emplace_back(std::make_unique<DPBuffer>(time, 1.0));
    unweighted_buffers.back()->GetData().resize(kShape);
    FillRegularData(unweighted_buffers.back()->GetData(""));

    for (const std::string& name : kDirectionNames) {
      unweighted_buffers.back()->AddData(name);
      FillRegularData(unweighted_buffers.back()->GetData(name));
    }

    unweighted_buffers.back()->GetWeights().resize(kShape);
    unweighted_buffers.back()->GetWeights().fill(1.0f);

    unweighted_buffers.back()->GetFlags().resize(kShape);
    unweighted_buffers.back()->GetFlags().fill(false);
  }

  dp3::ddecal::AssignAndWeight(unweighted_buffers, kDirectionNames,
                               weighted_buffers, false, false);

  const dp3::ddecal::SolveData data(
      weighted_buffers, kDirectionNames, kNChannelBlocks, kNAntennas,
      kNSolutionsPerDirection, kAntennas1, kAntennas2);
  BOOST_TEST_REQUIRE(data.NChannelBlocks() == kNChannelBlocks);

  for (size_t ch_block = 0; ch_block < kNChannelBlocks; ++ch_block) {
    const ChannelBlockData& cb_data = data.ChannelBlock(ch_block);
    BOOST_TEST_REQUIRE(cb_data.NDirections() == kDirectionNames.size());
    // SolveData does not include baselines with auto-correlations,
    // so it should not contain the last baseline.
    BOOST_TEST_REQUIRE(cb_data.NVisibilities() ==
                       kNTimes * kNSolveDataBaselines *
                           kExpectedChannelBlockSizes[ch_block]);
    for (size_t ant = 0; ant < kNAntennas; ++ant) {
      BOOST_TEST(cb_data.NAntennaVisibilities(ant) ==
                 kNTimes * ((ant == 0) ? 2 : 1) *
                     kExpectedChannelBlockSizes[ch_block]);
    }

    BOOST_TEST(cb_data.NSolutionsForDirection(0) == 1u);
    BOOST_TEST(cb_data.NSolutionsForDirection(1) == 2u);
    BOOST_TEST(cb_data.SolutionIndex(0, 0) == 0);
    BOOST_TEST(cb_data.SolutionIndex(1, 0) == 1);
    BOOST_TEST(cb_data.SolutionIndex(1, cb_data.NVisibilities() - 1) == 2);
    BOOST_TEST(cb_data.NSolutionVisibilities(0, 0) ==
               kNTimes * kNSolveDataBaselines *
                   kExpectedChannelBlockSizes[ch_block]);
    BOOST_TEST(cb_data.NSolutionVisibilities(1, 1) ==
               kNTimes * kNSolveDataBaselines *
                   kExpectedChannelBlockSizes[ch_block] / 2);
    BOOST_TEST(cb_data.NSolutionVisibilities(1, 2) ==
               kNTimes * kNSolveDataBaselines *
                   kExpectedChannelBlockSizes[ch_block] / 2);
    for (size_t v = 1; v < cb_data.NVisibilities(); ++v) {
      BOOST_TEST(cb_data.SolutionIndex(0, v) == 0);
      BOOST_TEST(cb_data.SolutionIndex(1, v) >=
                 cb_data.SolutionIndex(1, v - 1));
    }
  }
}

BOOST_AUTO_TEST_CASE(bda) {
  // The BDA data from the SolverTester is too complex for a simple test.
  // -> Use a simpler BdaBuffer with three baselines:
  // - An 'averaged' baseline: Single BDA row with kNAveragedChannels channels.
  // - A 'normal' baseline: Two BDA rows with successive time values and
  //   kNAllChannels channels each.
  // - An auto-correlation baseline: Same layout as the averaged baseline.
  //   SolveData should skip all auto-correlations.
  const std::vector<size_t> kNRowsPerBaseline{2, 4, 2};
  const size_t kNAllChannels = 7;
  const size_t kNAveragedChannels = 4;
  const size_t kBdaBufferSize =
      (2 * kNAveragedChannels + 2 * kNAllChannels) * kNPolarizations * 2;
  const size_t kNDirections = 2;
  const double kUnitTimeInterval = 2.0;
  const double kSolverInterval = kUnitTimeInterval * 4;

  // Number of visibilities for a baseline in one row
  // Outer index: channel block index; inner index: baseline index.
  const std::vector<std::vector<size_t>> kNVisibilitiesPerBaseline{{2, 3},
                                                                   {2, 4}};

  BdaBuffer::Fields bda_fields(true);

  auto bda_data_buffer =
      std::make_unique<BdaBuffer>(kBdaBufferSize, bda_fields);
  bda_fields.weights = false;
  FillBdaBuffer(*bda_data_buffer, kNAveragedChannels, kNAllChannels,
                kUnitTimeInterval);

  std::vector<std::unique_ptr<BdaBuffer>> bda_model_buffers;
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    bda_model_buffers.push_back(
        std::make_unique<BdaBuffer>(kBdaBufferSize, bda_fields));
    FillBdaBuffer(*bda_model_buffers.back(), kNAveragedChannels, kNAllChannels,
                  kUnitTimeInterval);
  }

  dp3::ddecal::BdaSolverBuffer solver_buffer(kNDirections, 0.0, kSolverInterval,
                                             kNBaselines);
  solver_buffer.AppendAndWeight(std::move(bda_data_buffer),
                                std::move(bda_model_buffers));
  BOOST_TEST_REQUIRE(solver_buffer.GetDataRows().size() ==
                     kNRowsPerBaseline[0] + kNRowsPerBaseline[1] +
                         kNRowsPerBaseline[2]);

  const std::vector<size_t> solutions_per_direction{1, 2};
  const dp3::ddecal::SolveData solve_data(solver_buffer, kNChannelBlocks,
                                          kNAntennas, solutions_per_direction,
                                          kAntennas1, kAntennas2, true);
  BOOST_TEST_REQUIRE(solve_data.NChannelBlocks() == kNChannelBlocks);

  // There should be two baselines:
  // - A baseline with 4 visibilities per row, of which 2 in channelblock 0.
  //   This baseline has 2 timesteps, which yields 4 visibilities in total.
  // - A baseline with 7 visibilities per row, of which 3 in channelblock 0.
  //   This baseline has 4 timesteps, which yields 12 visibilities in total.
  // Hence, we expect 16 visibilities:
  BOOST_TEST(solve_data.ChannelBlock(0).NSolutionVisibilities(0, 0) == 16);
  // Same as for solution 0, but this interval was split in two (see
  // solutions_per_direction):
  BOOST_TEST(solve_data.ChannelBlock(0).NSolutionVisibilities(1, 1) == 8);
  BOOST_TEST(solve_data.ChannelBlock(0).NSolutionVisibilities(1, 2) == 8);

  // There should be two baselines:
  // - Baseline 0 has 4 visibilities per row, of which 2 in channelblock 1.
  //   This baseline has 2 timesteps, which yields 4 visibilities in total.
  // - Baseline 1 has 7 visibilities per row, of which 4 in channelblock 1.
  //   This baseline has 4 timesteps, which yields 16 visibilities in total.
  // Hence, we expect 20 visibilities in a full interval.
  BOOST_TEST(solve_data.ChannelBlock(1).NSolutionVisibilities(0, 0) == 20);
  BOOST_TEST(solve_data.ChannelBlock(1).NSolutionVisibilities(1, 1) == 10);
  BOOST_TEST(solve_data.ChannelBlock(1).NSolutionVisibilities(1, 2) == 10);

  for (size_t ch_block = 0; ch_block < kNChannelBlocks; ++ch_block) {
    const ChannelBlockData& cb_data = solve_data.ChannelBlock(ch_block);
    BOOST_TEST_REQUIRE(cb_data.NDirections() == kNDirections);

    // SolveData does not include baselines with auto-correlations,
    // so it should not contain the last baseline, with index 2.
    const size_t n_visibilities = kNVisibilitiesPerBaseline[ch_block][0] +
                                  2 * kNVisibilitiesPerBaseline[ch_block][1];
    BOOST_TEST(cb_data.NVisibilities() == 2 * n_visibilities);
    BOOST_TEST(cb_data.NAntennaVisibilities(0) == 2 * n_visibilities);
    BOOST_TEST(cb_data.NAntennaVisibilities(1) ==
               2 * kNVisibilitiesPerBaseline[ch_block][0]);
    BOOST_TEST(cb_data.NAntennaVisibilities(2) ==
               4 * kNVisibilitiesPerBaseline[ch_block][1]);

    for (size_t v = 0; v < 2 * n_visibilities; ++v) {
      const aocommon::MC2x2F& data = cb_data.Visibility(v);
      const aocommon::MC2x2F& model_data = cb_data.ModelVisibility(0, v);

      const size_t baseline =
          ((v % n_visibilities) < kNVisibilitiesPerBaseline[ch_block][0]) ? 0
                                                                          : 1;
      BOOST_TEST(cb_data.Antenna1Index(v) == 0);
      BOOST_TEST(cb_data.Antenna2Index(v) == ((baseline == 0) ? 1 : 2));
      BOOST_TEST(cb_data.SolutionIndex(0, v) == 0);
      const size_t dir1_solution_index = v < n_visibilities ? 1 : 2;
      BOOST_TEST(cb_data.SolutionIndex(1, v) == dir1_solution_index);

      // Determine the BDA row and the data offset in that row.
      constexpr size_t kRowsBlock1[] = {0, 0, 2, 2, 2, 3, 3, 3,
                                        4, 4, 6, 6, 6, 7, 7, 7};
      constexpr size_t kRowsBlock2[] = {0, 0, 2, 2, 2, 2, 3, 3, 3, 3,
                                        4, 4, 6, 6, 6, 6, 7, 7, 7, 7};
      constexpr size_t kOffsets1[] = {0, 1, 0, 1, 2, 0, 1, 2,
                                      0, 1, 0, 1, 2, 0, 1, 2};
      constexpr size_t kOffsets2[] = {0, 1, 0, 1, 2, 3, 0, 1, 2, 3,
                                      0, 1, 0, 1, 2, 3, 0, 1, 2, 3};
      const size_t row = ch_block == 0 ? kRowsBlock1[v] : kRowsBlock2[v];
      size_t offset = ch_block == 0 ? kOffsets1[v] : kOffsets2[v];
      const size_t n_channels =
          (baseline == 0) ? kNAveragedChannels : kNAllChannels;
      const size_t first_channel = ch_block * n_channels / kNChannelBlocks;
      offset = (offset + first_channel) * kNPolarizations;

      const std::complex<float>* const data_ptr =
          solver_buffer.GetDataRows()[row]->data + offset;
      const std::complex<float>* const model_data_ptr =
          solver_buffer.GetModelDataRows(0)[row]->data + offset;
      for (size_t p = 0; p < kNPolarizations; ++p) {
        BOOST_TEST(data.Get(p) == data_ptr[p]);
        BOOST_TEST(model_data.Get(p) == model_data_ptr[p]);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
