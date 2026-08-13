// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef USE_FAST_PREDICT

#include "steps/FastPredict.h"

#include <regex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <predict/BeamResponse.h>

#include "base/DP3.h"

#include "common/ParameterSet.h"
#include "steps/ApplyCal.h"
#include "steps/NullStep.h"
#include "steps/OnePredict.h"

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
// Realistic MJD-seconds timestamp (~2013) matching tNDPPP-generic.MS,
// required so EveryBeam beam calculations fall within the IERS table range.
constexpr double kBeamTestTime = 4.87128e+09;

dp3::common::ParameterSet MakeBeamParset() {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("usebeammodel", "True");
  parset.add("beammode", "array_factor");
  parset.add("beam_interval", "120");
  return parset;
}

struct BeamConsistencyFixture {
 public:
  BeamConsistencyFixture() {
    const dp3::common::ParameterSet predict_parset = MakeBeamParset();

    one_predict = std::make_shared<dp3::steps::OnePredict>(
        predict_parset, "", std::vector<std::string>());
    one_result = std::make_shared<dp3::steps::ResultStep>();
    one_predict->setNextStep(one_result);

    fast_predict = std::make_shared<FastPredict>(predict_parset, "",
                                                 std::vector<std::string>());
    fast_result = std::make_shared<dp3::steps::ResultStep>();
    fast_predict->setNextStep(fast_result);

    dp3::base::DPInfo info(kNCorr, kNChan, "HBA_DUAL_INNER");
    info.setMsName("tNDPPP-generic.MS");
    info.setTimes(kBeamTestTime, kBeamTestTime + 9.0, 1.0);
    const std::vector<int> ant1{0, 0, 1};
    const std::vector<int> ant2{1, 2, 2};
    const std::vector<std::string> ant_names{"CS001HBA0", "CS002HBA0",
                                             "CS002HBA1"};

    // Real LOFAR HBA ITRF positions are required for EveryBeam to produce a
    // time-varying array factor. Zero positions give a degenerate result.
    casacore::Vector<double> vals(3);
    std::vector<casacore::MPosition> ant_pos(3);
    vals[0] = 3828763;
    vals[1] = 442449;
    vals[2] = 5064923;
    ant_pos[0] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828746;
    vals[1] = 442592;
    vals[2] = 5064924;
    ant_pos[1] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828827;
    vals[1] = 442642;
    vals[2] = 5064875;
    ant_pos[2] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    info.setAntennas(ant_names, std::vector<double>(3, 70.0), ant_pos, ant1,
                     ant2);
    info.setChannels(std::vector<double>(kNChan, 120.0e6),
                     std::vector<double>(kNChan, 3.0e6));
    one_predict->setInfo(info);
    fast_predict->setInfo(info);
  }

  std::shared_ptr<dp3::steps::OnePredict> one_predict;
  std::shared_ptr<FastPredict> fast_predict;
  std::shared_ptr<dp3::steps::ResultStep> one_result;
  std::shared_ptr<dp3::steps::ResultStep> fast_result;
};

class FastPredictFixture {
 public:
  FastPredictFixture() : predict_() {
    dp3::common::ParameterSet parset;
    parset.add("fixture.sourcedb", dp3::steps::test::kPredictSkyModel);
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
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("operation", operation);
  parset.add("usefastpredict", "True");
  const FastPredict predict(parset, "", {});
  BOOST_TEST(predict.getRequiredFields() ==
             (Step::kDataField | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() == Step::kDataField);
}

BOOST_FIXTURE_TEST_CASE(fields_applycal, dp3::steps::test::H5ParmFixture) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("usefastpredict", "True");
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  const FastPredict predict(parset, "", {});

  // FastPredict uses ApplyCal which has a OneApplyCal sub-step as next step.
  const std::shared_ptr<dp3::steps::ApplyCal> apply_cal =
      std::make_shared<dp3::steps::ApplyCal>(parset, "applycal.", true);

  const dp3::common::Fields apply_cal_required =
      dp3::base::GetChainRequiredFields(apply_cal);
  // TODO(AST-1033) Determine ApplyCal provided fields using generic DP3
  // functions.
  const dp3::common::Fields apply_cal_provided =
      apply_cal->getNextStep()->getProvidedFields();
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
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("usefastpredict", "True");
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  parset.add("operation", operation);
  const FastPredict predict(parset, "", {});

  // When operation is "add" or "subtract", FastPredict only combines the
  // required fields of its ApplyCal sub-step.
  const std::shared_ptr<dp3::steps::ApplyCal> apply_cal =
      std::make_shared<dp3::steps::ApplyCal>(parset, "applycal.", true);

  const dp3::common::Fields apply_cal_required =
      dp3::base::GetChainRequiredFields(apply_cal);

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
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
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

BOOST_AUTO_TEST_CASE(full_beam_sparse_station_ids) {
  const std::vector<predict::Baseline> baselines{{15, 16}};
  const xt::xtensor<double, 1> frequencies = {120.0e6, 121.0e6};

  // For better locality, the baseline and polarization dimensions are reversed
  // in the ApplyBeamToDataAndAdd method.
  Buffer4D direction_buffer({baselines.size(), 4, 2, frequencies.size()}, 0.0f);
  Buffer4D model_data({4, baselines.size(), 2, frequencies.size()}, 0.0f);

  // Use a diagonal 2x2 visibility of 1+0i for all channels.
  for (size_t ch = 0; ch < frequencies.size(); ++ch) {
    direction_buffer(0, 0, 0, ch) = 1.0f;
    direction_buffer(0, 3, 0, ch) = 1.0f;
  }

  // Beam values are produced in compact (unique-station) indexing order.
  xt::xtensor<float, 4> beam_values({2, 4, 2, frequencies.size()}, 0.0f);
  for (size_t ch = 0; ch < frequencies.size(); ++ch) {
    beam_values(0, 0, 0, ch) = 2.0f;
    beam_values(0, 3, 0, ch) = 2.0f;
    beam_values(1, 0, 0, ch) = 3.0f;
    beam_values(1, 3, 0, ch) = 3.0f;
  }

  predict::BeamResponsePlan beam_plan;
  beam_plan.SetTime(0.0);
  beam_plan.SetFieldId(0);
  beam_plan.SetBeamMode(everybeam::BeamMode::kFull);
  beam_plan.ApplyBeamToDataAndAdd(baselines, frequencies, direction_buffer,
                                  model_data, beam_values);

  for (size_t ch = 0; ch < frequencies.size(); ++ch) {
    BOOST_CHECK_CLOSE(model_data(0, 0, 0, ch), 6.0f, 1.0e-6);
    BOOST_CHECK_SMALL(model_data(1, 0, 0, ch), 1.0e-6f);
    BOOST_CHECK_SMALL(model_data(2, 0, 0, ch), 1.0e-6f);
    BOOST_CHECK_CLOSE(model_data(3, 0, 0, ch), 6.0f, 1.0e-6);
    BOOST_CHECK_SMALL(model_data(0, 0, 1, ch), 1.0e-6f);
    BOOST_CHECK_SMALL(model_data(1, 0, 1, ch), 1.0e-6f);
    BOOST_CHECK_SMALL(model_data(2, 0, 1, ch), 1.0e-6f);
    BOOST_CHECK_SMALL(model_data(3, 0, 1, ch), 1.0e-6f);
  }
}

BOOST_FIXTURE_TEST_CASE(beam_time_consistency_with_array_factor,
                        BeamConsistencyFixture) {
  using BufferPtr = std::unique_ptr<dp3::base::DPBuffer>;
  BufferPtr input =
      CreateBuffer(kBeamTestTime, kInterval, kNBaselines, kChannelCounts, 1.0f);

  one_predict->process(std::make_unique<dp3::base::DPBuffer>(*input));
  fast_predict->process(std::make_unique<dp3::base::DPBuffer>(*input));

  BufferPtr one_output = one_result->take();
  BufferPtr fast_output = fast_result->take();

  const dp3::base::DPBuffer::DataType& one_data = one_output->GetData();
  const dp3::base::DPBuffer::DataType& fast_data = fast_output->GetData();

  for (std::size_t i = 0; i < one_data.size(); ++i) {
    BOOST_CHECK_CLOSE(fast_data.data()[i], one_data.data()[i], 1.0e-4f);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_FAST_PREDICT
