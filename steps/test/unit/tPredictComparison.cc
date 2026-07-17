// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef USE_FAST_PREDICT

#include "steps/FastPredict.h"
#include "steps/OnePredict.h"

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "common/ParameterSet.h"
#include "common/test/unit/fixtures/fDirectory.h"
#include "steps/test/unit/tPredict.h"
#include "steps/ResultStep.h"

using aocommon::Logger;

namespace {

constexpr unsigned int kNCorr = 4;
constexpr unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
constexpr double kStartTime = 0.0;
constexpr double kInterval = 1.0;
constexpr std::size_t kNBaselines = 3;

static std::unique_ptr<dp3::base::DPBuffer> CreateBuffer(
    const double time, const double interval, std::size_t n_baselines,
    const std::vector<std::size_t>& channel_counts, const float base_value,
    const float weight = 1.0) {
  const std::array<std::size_t, 3> kShape{n_baselines, channel_counts.size(),
                                          kNCorr};

  auto buffer = std::make_unique<dp3::base::DPBuffer>(time, interval);
  buffer->GetData().resize(kShape);
  buffer->GetWeights().resize(kShape);
  buffer->GetFlags().resize(kShape);
  buffer->GetUvw().resize({n_baselines, 3});

  buffer->GetFlags().fill(false);
  buffer->GetWeights().fill(weight);

  for (std::size_t baseline = 0; baseline < n_baselines; ++baseline) {
    const float baseline_value = (baseline * 100.0f) + (base_value / weight);

    std::size_t channel = 0;
    float channel_value = baseline_value;
    for (std::size_t channel_count : channel_counts) {
      const float value = channel_value + 5.0f * (channel_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        buffer->GetData()(baseline, channel, corr) = value + corr;
        buffer->GetWeights()(baseline, channel, corr) *= channel_count;
      }
      ++channel;
      channel_value += channel_count * 10.0f;
    }
    buffer->GetUvw()(baseline, 0) = baseline_value + 0.0;
    buffer->GetUvw()(baseline, 1) = baseline_value + 1.0;
    buffer->GetUvw()(baseline, 2) = baseline_value + 2.0;
  }

  return buffer;
}

template <typename TPredict>
void SetInfo(const std::shared_ptr<TPredict>& predict) {
  dp3::base::DPInfo info(kNCorr, kNChan);
  info.setTimes(0.5, 9.5, 1.0);

  const std::vector<int> kAnt1{0, 0, 1};
  const std::vector<int> kAnt2{1, 2, 2};
  const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
  const std::vector<double> kAntDiam(3, 1.0);
  const std::vector<casacore::MPosition> kAntPos(3);
  info.setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1, kAnt2);

  std::vector<double> chan_freqs(kNChan, 10.0e6);
  std::vector<double> chan_widths(kNChan, 3.0e6);

  info.setChannels(std::move(chan_freqs), std::move(chan_widths));
  predict->setInfo(info);
}

double ComputeMaxAbsDiff(const dp3::base::DPBuffer::DataType& lhs,
                         const dp3::base::DPBuffer::DataType& rhs) {
  const auto& lhs_shape = lhs.shape();
  const auto& rhs_shape = rhs.shape();
  BOOST_REQUIRE_EQUAL(lhs_shape.size(), rhs_shape.size());
  for (size_t dim = 0; dim < lhs_shape.size(); ++dim) {
    BOOST_REQUIRE_EQUAL(lhs_shape[dim], rhs_shape[dim]);
  }

  double max_abs_diff = 0.0;
  for (size_t bl = 0; bl < lhs.shape()[0]; ++bl) {
    for (size_t ch = 0; ch < lhs.shape()[1]; ++ch) {
      for (size_t cr = 0; cr < lhs.shape()[2]; ++cr) {
        const double abs_diff =
            std::abs(static_cast<std::complex<double>>(lhs(bl, ch, cr)) -
                     static_cast<std::complex<double>>(rhs(bl, ch, cr)));
        max_abs_diff = std::max(max_abs_diff, abs_diff);
      }
    }
  }

  return max_abs_diff;
}

}  // namespace

BOOST_FIXTURE_TEST_SUITE(predictcomparison, dp3::common::test::FixtureDirectory)

BOOST_AUTO_TEST_CASE(gaussian_only_fast_vs_onepredict_difference) {
  dp3::common::ParameterSet one_parset;
  one_parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);

  dp3::common::ParameterSet fast_parset;
  fast_parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  fast_parset.add("usefastpredict", "True");

  const std::vector<std::string> kSourcePatterns{"0010.1+3045"};

  auto one_predict =
      std::make_shared<dp3::steps::OnePredict>(one_parset, "", kSourcePatterns);
  auto fast_predict = std::make_shared<dp3::steps::FastPredict>(
      fast_parset, "", kSourcePatterns);

  auto one_result = std::make_shared<dp3::steps::ResultStep>();
  auto fast_result = std::make_shared<dp3::steps::ResultStep>();

  one_predict->setNextStep(one_result);
  fast_predict->setNextStep(fast_result);

  SetInfo(one_predict);
  SetInfo(fast_predict);

  one_predict->process(CreateBuffer(kStartTime * kInterval, kInterval,
                                    kNBaselines, kChannelCounts, 0.0f));
  fast_predict->process(CreateBuffer(kStartTime * kInterval, kInterval,
                                     kNBaselines, kChannelCounts, 0.0f));

  std::unique_ptr<dp3::base::DPBuffer> one_buffer = one_result->take();
  std::unique_ptr<dp3::base::DPBuffer> fast_buffer = fast_result->take();

  const double max_abs_diff =
      ComputeMaxAbsDiff(one_buffer->GetData(), fast_buffer->GetData());
  Logger::Info << "Gaussian-only OnePredict vs FastPredict max abs diff: "
               << max_abs_diff << std::endl;

  // Highlight current behavior: FastPredict and OnePredict diverge on this
  // Gaussian patch by more than pure floating-point epsilon.
  BOOST_TEST(max_abs_diff > 1.0e-7);
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_FAST_PREDICT
