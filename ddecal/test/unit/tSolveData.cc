// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../gain_solvers/SolveData.h"

#include "SolverTester.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_unique.hpp>

#include <algorithm>
#include <iterator>
#include <random>

using dp3::base::BDABuffer;
using dp3::ddecal::test::SolverTester;
using ChannelBlockData = dp3::ddecal::SolveData::ChannelBlockData;

namespace {

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

void FillBdaBuffer(BDABuffer& buffer, size_t avg_channels,
                   size_t all_channels) {
  const size_t kNPolarizations = 4;
  std::vector<std::complex<float>> data;
  const std::vector<float> weights(all_channels * kNPolarizations, 1.0f);

  // Add an averaged row for baseline 0.
  FillRandomData(data, avg_channels * kNPolarizations);
  BOOST_REQUIRE(buffer.AddRow(1.0, 2.0, 2.0, 0, avg_channels, kNPolarizations,
                              data.data(), nullptr, weights.data()));

  // Add non-averaged rows for baseline 1.
  for (int j = 0; j < 2; ++j) {
    FillRandomData(data, all_channels * kNPolarizations);
    BOOST_REQUIRE(buffer.AddRow(0.5 + j, 1.0, 1.0, 1, all_channels,
                                kNPolarizations, data.data(), nullptr,
                                weights.data()));
  }
}

}  // namespace

// SolveData should only copy data, so the tolerance is zero for these tests.
BOOST_AUTO_TEST_SUITE(solve_data, *boost::unit_test::tolerance(0.0f))

BOOST_FIXTURE_TEST_CASE(regular, SolverTester) {
  // Test that all data from a SolverBuffer ends up in a SolveData structure.
  // A 'complete' SolverBuffer, with kNRegularTimes time steps, is not needed.
  const size_t kNTimes = 2;
  const dp3::ddecal::SolverBuffer& solver_buffer = FillData(kNTimes);
  const dp3::ddecal::SolveData data(solver_buffer, kNChannelBlocks,
                                    kNDirections, kNAntennas, Antennas1(),
                                    Antennas2());
  BOOST_REQUIRE(data.NChannelBlocks() == kNChannelBlocks);

  const std::vector<size_t> kExpectedChannelBlockSizes{2, 3, 2, 3};

  for (size_t ch_block = 0; ch_block < kNChannelBlocks; ++ch_block) {
    const ChannelBlockData& cb_data = data.ChannelBlock(ch_block);
    BOOST_REQUIRE(cb_data.NDirections() == kNDirections);
    // SolveData does not include baselines with auto-correlations. However,
    // the SolverTester creates data without auto-correlations, so the number
    // of baselines remains equal.
    BOOST_TEST(cb_data.NVisibilities() ==
               kNTimes * kNBaselines * kExpectedChannelBlockSizes[ch_block]);
    for (size_t ant = 0; ant < kNAntennas; ++ant) {
      BOOST_TEST(cb_data.NAntennaVisibilities(ant) ==
                 kNTimes * (kNAntennas - 1) *
                     kExpectedChannelBlockSizes[ch_block]);
    }

    const size_t first_channel = ch_block * kNChannels / kNChannelBlocks;
    for (size_t v = 0; v < cb_data.NVisibilities(); ++v) {
      const size_t time =
          v / (kNBaselines * kExpectedChannelBlockSizes[ch_block]);
      const size_t baseline =
          (v / kExpectedChannelBlockSizes[ch_block]) % kNBaselines;
      const size_t channel =
          first_channel + v % kExpectedChannelBlockSizes[ch_block];

      BOOST_TEST(cb_data.Antenna1Index(v) == size_t(Antennas1()[baseline]));
      BOOST_TEST(cb_data.Antenna2Index(v) == size_t(Antennas2()[baseline]));

      const std::complex<float>* const expected_data =
          solver_buffer.DataPointer(time, baseline, channel);
      const aocommon::MC2x2F& data = cb_data.Visibility(v);
      for (size_t pol = 0; pol < kNPolarizations; ++pol) {
        BOOST_TEST(expected_data[pol] == data[pol]);
      }

      for (size_t direction = 0; direction < kNDirections; ++direction) {
        const std::complex<float>* const expected_model_data =
            solver_buffer.ModelDataPointer(time, direction, baseline, channel);
        const aocommon::MC2x2F& model_data =
            cb_data.ModelVisibility(direction, v);
        const aocommon::MC2x2F& model_vector_data =
            cb_data.ModelVisibilityVector(direction)[v];

        for (size_t pol = 0; pol < kNPolarizations; ++pol) {
          BOOST_TEST(expected_model_data[pol] == model_data[pol]);
          BOOST_TEST(expected_model_data[pol] == model_vector_data[pol]);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(bda) {
  // The BDA data from the SolverTester is too complex for a simple test.
  // -> Use a simpler BDABuffer with two baselines:
  // - An 'averaged' baseline: Single BDA row with kAveragedChannels channels.
  // - A 'normal' baseline: Two BDA rows with kAllChannels channels each.
  const size_t kNPolarizations = 4;
  const size_t kNChannelBlocks = 2;
  const size_t kAllChannels = 7;
  const size_t kAveragedChannels = 4;
  const size_t kBdaBufferSize =
      (kAveragedChannels + 2 * kAllChannels) * kNPolarizations;
  const size_t kNDirections = 1;
  const size_t kNBaselines = 2;
  const size_t kNAntennas = 3;
  const std::vector<int> kAntennas1 = {0, 0};
  const std::vector<int> kAntennas2 = {1, 2};

  // Outer index: channel block index; inner index: baseline index.
  const std::vector<std::vector<size_t>> kBaselineVisibilities{{2, 3}, {2, 4}};

  BDABuffer::Fields bda_fields(false);
  bda_fields.data = true;
  bda_fields.weights = true;
  auto bda_data_buffer =
      boost::make_unique<BDABuffer>(kBdaBufferSize, bda_fields);
  bda_fields.weights = false;
  std::vector<std::unique_ptr<BDABuffer>> bda_model_buffers;
  bda_model_buffers.push_back(
      boost::make_unique<BDABuffer>(kBdaBufferSize, bda_fields));

  FillBdaBuffer(*bda_data_buffer, kAveragedChannels, kAllChannels);
  FillBdaBuffer(*bda_model_buffers.back(), kAveragedChannels, kAllChannels);

  dp3::ddecal::BdaSolverBuffer solver_buffer(kNDirections, -1.0, 10.0,
                                             kNBaselines);
  solver_buffer.AppendAndWeight(std::move(bda_data_buffer),
                                std::move(bda_model_buffers));
  BOOST_REQUIRE(solver_buffer.GetDataRows().size() == size_t(1 + 2));

  const dp3::ddecal::SolveData solve_data(solver_buffer, kNChannelBlocks,
                                          kNDirections, kNAntennas, kAntennas1,
                                          kAntennas2);
  BOOST_REQUIRE(solve_data.NChannelBlocks() == kNChannelBlocks);

  for (size_t ch_block = 0; ch_block < kNChannelBlocks; ++ch_block) {
    const ChannelBlockData& cb_data = solve_data.ChannelBlock(ch_block);
    BOOST_REQUIRE(cb_data.NDirections() == kNDirections);

    const size_t n_visibilities = kBaselineVisibilities[ch_block][0] +
                                  2 * kBaselineVisibilities[ch_block][1];
    BOOST_TEST(cb_data.NVisibilities() == n_visibilities);
    BOOST_TEST(cb_data.NAntennaVisibilities(0) == n_visibilities);
    BOOST_TEST(cb_data.NAntennaVisibilities(1) ==
               kBaselineVisibilities[ch_block][0]);
    BOOST_TEST(cb_data.NAntennaVisibilities(2) ==
               2 * kBaselineVisibilities[ch_block][1]);

    for (size_t v = 0; v < n_visibilities; ++v) {
      const aocommon::MC2x2F& data = cb_data.Visibility(v);
      const aocommon::MC2x2F& model_data = cb_data.ModelVisibility(0, v);
      const aocommon::MC2x2F& model_vector_data =
          cb_data.ModelVisibilityVector(0)[v];

      const size_t baseline = (v < kBaselineVisibilities[ch_block][0]) ? 0 : 1;
      BOOST_TEST(cb_data.Antenna1Index(v) == 0);
      BOOST_TEST(cb_data.Antenna2Index(v) == ((baseline == 0) ? 1 : 2));

      // Determine the BDA row and the data offset in that row.
      size_t row = 0;
      size_t offset = v;
      if (baseline == 1) {
        ++row;
        offset -= kBaselineVisibilities[ch_block][0];
        if (offset >= kBaselineVisibilities[ch_block][1]) {
          ++row;
          offset -= kBaselineVisibilities[ch_block][1];
        }
      }
      const size_t n_channels =
          (baseline == 0) ? kAveragedChannels : kAllChannels;
      const size_t first_channel = ch_block * n_channels / kNChannelBlocks;
      offset = (offset + first_channel) * kNPolarizations;

      const std::complex<float>* const data_ptr =
          solver_buffer.GetDataRows()[row]->data + offset;
      const std::complex<float>* const model_data_ptr =
          solver_buffer.GetModelDataRows(0)[row]->data + offset;
      for (size_t p = 0; p < kNPolarizations; ++p) {
        BOOST_TEST(data[p] == data_ptr[p]);
        BOOST_TEST(model_data[p] == model_data_ptr[p]);
        BOOST_TEST(model_vector_data[p] == model_data_ptr[p]);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
