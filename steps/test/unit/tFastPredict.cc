// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef USE_FAST_PREDICT

#include "steps/FastPredict.h"

#include <regex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "base/DP3.h"

#include "common/ParameterSet.h"
#include "steps/ApplyCal.h"
#include "steps/NullStep.h"

#include "tPredict.h"
#include "H5ParmFixture.h"

using dp3::steps::FastPredict;
using dp3::steps::Step;

namespace {

// Constants copy pasted from steps/test/unit/tIDGPredict.cc
constexpr unsigned int kNCorr = 4;
constexpr unsigned int kNChan = 5;
const std::vector<std::size_t> kChannelCounts(kNChan, 1);
constexpr double kStartTime = 0.0;
constexpr double kInterval = 1.0;
constexpr std::size_t kNBaselines = 3;

class FastPredictFixture {
 public:
  FastPredictFixture() : predict_() {
    dp3::common::ParameterSet parset;
    parset.add("fixture.sourcedb", dp3::steps::test::kPredictSkymodel);
    predict_ = std::make_shared<FastPredict>(parset, "fixture.",
                                             std::vector<std::string>());
    predict_->setNextStep(std::make_shared<dp3::steps::NullStep>());
    SetInfo(predict_);
  }

  static void SetInfo(std::shared_ptr<FastPredict> predict) {
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

 protected:
  std::shared_ptr<FastPredict> predict_;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(fastpredict)

BOOST_FIXTURE_TEST_CASE(constructor, FastPredictFixture) {
  // Nothing to do: The fixture calls the constructor.
}

BOOST_FIXTURE_TEST_CASE(getfirstdirection, FastPredictFixture) {
  const dp3::base::Direction first_direction = predict_->GetFirstDirection();

  BOOST_CHECK_CLOSE(first_direction.ra,
                    dp3::steps::test::kExpectedFirstDirection.ra, 1.0e-3);
  BOOST_CHECK_CLOSE(first_direction.dec,
                    dp3::steps::test::kExpectedFirstDirection.dec, 1.0e-3);
}

BOOST_FIXTURE_TEST_CASE(fields_defaults, FastPredictFixture) {
  BOOST_TEST(predict_->getRequiredFields() == Step::kUvwField);
  BOOST_TEST(predict_->getProvidedFields() == Step::kDataField);
}

BOOST_DATA_TEST_CASE(fields_add_subtract,
                     boost::unit_test::data::make({"add", "subtract"}),
                     operation) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("operation", operation);
  parset.add("usefastpredict", "True");
  const FastPredict predict(parset, "", {});
  BOOST_TEST(predict.getRequiredFields() ==
             (Step::kDataField | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() == Step::kDataField);
}

BOOST_FIXTURE_TEST_CASE(fields_applycal, dp3::steps::test::H5ParmFixture) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("usefastpredict", "True");
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  const FastPredict predict(parset, "", {});

  // FastPredict uses ApplyCal which has a OneApplyCal sub-step as next step.
  const dp3::steps::ApplyCal apply_cal(parset, "applycal.", true);

  const dp3::common::Fields apply_cal_required =
      dp3::base::GetChainRequiredFields(
          std::make_shared<dp3::steps::ApplyCal>(apply_cal));
  // TODO(AST-1033) Determine ApplyCal provided fields using generic DP3
  // functions.
  const dp3::common::Fields apply_cal_provided =
      apply_cal.getNextStep()->getProvidedFields();
  BOOST_TEST(predict.getRequiredFields() ==
             (apply_cal_required | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() ==
             (apply_cal_provided | Step::kDataField));
}

BOOST_DATA_TEST_CASE_F(dp3::steps::test::H5ParmFixture,
                       fields_applycal_add_subtract,
                       boost::unit_test::data::make({"add", "subtract"}),
                       operation) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("usefastpredict", "True");
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  parset.add("operation", operation);
  const FastPredict predict(parset, "", {});

  // When operation is "add" or "subtract", FastPredict only combines the
  // required fields of its ApplyCal sub-step.
  const dp3::steps::ApplyCal apply_cal(parset, "applycal.", true);

  const dp3::common::Fields apply_cal_required =
      dp3::base::GetChainRequiredFields(
          std::make_shared<dp3::steps::ApplyCal>(apply_cal));

  BOOST_TEST(predict.getRequiredFields() ==
             (apply_cal_required | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() == Step::kDataField);
}

/**
 * Create a buffer with artificial data values.
 * @param time Start time for the buffer.
 * @param interval Interval duration for the buffer.
 * @param n_baselines Number of baselines in the buffer.
 * @param base_value Base value for the data values, for distinguishing
 buffers.
 *        For distinguishing baselines, this function adds baseline_nr *
 100.0.
 *        When the buffer represents averaged data, the base_value should be
 *        the total of the base values of the original buffers.
 *        This function divides the base_value by the supplied weight so the
 *        caller does not have to do that division.
 * @param channel_counts List for generating channel data.
 *        For input buffers, this list should contain a 1 for each channel.
 *        When generating expected output data, this list should contain the
 *        number of averaged input buffers for each output buffer.
 * @param weight Weight value for the data values in the buffer.
 *
 * @note The function has been copied from @ref
 steps/test/unit/tIDGPredict.cc.
 */
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
    // Base value for this baseline.
    const float baseline_value = (baseline * 100.0) + (base_value / weight);

    std::size_t channel = 0;
    float channel_value = baseline_value;  // Base value for a group of channels
    for (std::size_t channel_count : channel_counts) {
      // For each channel, increase channel_value by 10.0.
      // When channel_count == 1, 'value' should equal channel_value.
      // When channel_count > 1, 'value' should be the average for multiple
      // channels.
      const float value = channel_value + 5.0 * (channel_count - 1);
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        buffer->GetData()(baseline, channel, corr) = value + corr;
        buffer->GetWeights()(baseline, channel, corr) *= channel_count;
      }
      ++channel;
      channel_value += channel_count * 10.0;
    }
    buffer->GetUvw()(baseline, 0) = baseline_value + 0.0;
    buffer->GetUvw()(baseline, 1) = baseline_value + 1.0;
    buffer->GetUvw()(baseline, 2) = baseline_value + 2.0;
  }

  return buffer;
}

BOOST_AUTO_TEST_CASE(outputmodelname) {
  std::unique_ptr<dp3::base::DPBuffer> input_buffer = CreateBuffer(
      kStartTime * kInterval, kInterval, kNBaselines, kChannelCounts, 0.);
  std::string output_model_name = "a_model_name";

  // Predict visibilities to main data buffer, replacing the input visibilities.
  // Make step chain
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkymodel);
  parset.add("usefastpredict", "True");
  auto predict =
      std::make_shared<FastPredict>(parset, "", std::vector<std::string>());
  auto predict_result = std::make_shared<dp3::steps::ResultStep>();
  predict->setNextStep(predict_result);
  FastPredictFixture::SetInfo(predict);

  // Process and verify
  predict->process(std::make_unique<dp3::base::DPBuffer>(*input_buffer));
  std::unique_ptr<dp3::base::DPBuffer> result_main = predict_result->take();

  BOOST_CHECK(!(result_main->HasData(output_model_name)));
  BOOST_CHECK(!xt::allclose(result_main->GetData(), input_buffer->GetData()));

  // Predict visibilities to an extra data buffer in the output DPBuffer.
  // Make step chain (with extra parset pair)
  parset.add("outputmodelname", output_model_name);
  predict =
      std::make_shared<FastPredict>(parset, "", std::vector<std::string>());
  predict_result = std::make_shared<dp3::steps::ResultStep>();
  predict->setNextStep(predict_result);
  FastPredictFixture::SetInfo(predict);

  // Process and verify
  predict->process(std::make_unique<dp3::base::DPBuffer>(*input_buffer));
  std::unique_ptr<dp3::base::DPBuffer> result_extra = predict_result->take();

  // Verify main data buffer still contains the original visibilities
  BOOST_CHECK(xt::allclose(result_extra->GetData(), input_buffer->GetData()));
  // Verify predicted visibilities are present and identical
  BOOST_CHECK(xt::allclose(result_extra->GetData(output_model_name),
                           result_main->GetData()));
  // Verify fields differ from the 'fields_defaults' test
  BOOST_TEST(predict->getRequiredFields() == Step::kUvwField);
  BOOST_TEST(predict->getProvidedFields() == dp3::common::Fields());
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_FAST_PREDICT
