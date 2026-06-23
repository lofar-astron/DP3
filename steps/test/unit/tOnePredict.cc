// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/OnePredict.h"

#include <regex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "base/DP3.h"

#include "common/ParameterSet.h"
#include "steps/ApplyCal.h"
#include "steps/NullStep.h"

#include "tPredict.h"
#include "H5ParmFixture.h"
#include "steps/test/unit/mock/MockTelescope.h"

using dp3::steps::ApplyCal;
using dp3::steps::OnePredict;
using dp3::steps::Step;

namespace {

// Constants copy pasted from steps/test/unit/tIDGPredict.cc
constexpr unsigned int kNCorr = 4;
constexpr double kTime = 0.0;
constexpr double kInterval = 1.0;
constexpr std::size_t kNBaselines = 3;
const std::vector<double> kChannelFrequencies = {
    1.1e7, 1.2e7, 1.3e7, 1.4e7, 1.5e7,
};
constexpr double kChannelWidth = 1.0e6;
const std::size_t kNChan = kChannelFrequencies.size();

class OnePredictFixture {
 public:
  OnePredictFixture() : predict_() {
    dp3::common::ParameterSet parset;
    parset.add("fixture.sourcedb", dp3::steps::test::kPredictSkyModel);
    predict_ = std::make_shared<OnePredict>(parset, "fixture.",
                                            std::vector<std::string>());
    predict_->setNextStep(std::make_shared<dp3::steps::NullStep>());
    predict_->updateInfo(MakeInfo());
  }

  static dp3::base::DPInfo MakeInfo() {
    dp3::base::DPInfo info(kNCorr, kNChan);
    info.setTimes(0.5, 9.5, 1.0);

    const std::vector<int> kAnt1{0, 0, 1};
    const std::vector<int> kAnt2{1, 2, 2};
    // For the reusebeammodel test, the antenna names should be real names
    // from the MS used in the test.
    const std::vector<std::string> kAntNames{"CS001HBA0", "CS002HBA0",
                                             "CS002HBA1"};
    const std::vector<double> kAntDiam(3, 1.0);
    const std::vector<casacore::MPosition> kAntPos(3);
    info.setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1, kAnt2);

    info.setChannels(std::vector<double>(kChannelFrequencies),
                     std::vector<double>(kNChan, kChannelWidth));
    return info;
  }

 protected:
  std::shared_ptr<OnePredict> predict_;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(onepredict)

BOOST_FIXTURE_TEST_CASE(constructor, OnePredictFixture) {
  // Nothing to do: The fixture calls the constructor.
}

BOOST_FIXTURE_TEST_CASE(getfirstdirection, OnePredictFixture) {
  const dp3::base::Direction first_direction = predict_->GetFirstDirection();

  BOOST_CHECK_CLOSE(first_direction.ra,
                    dp3::steps::test::kExpectedFirstDirection.ra, 1.0e-3);
  BOOST_CHECK_CLOSE(first_direction.dec,
                    dp3::steps::test::kExpectedFirstDirection.dec, 1.0e-3);
}

BOOST_FIXTURE_TEST_CASE(fields_defaults, OnePredictFixture) {
  BOOST_TEST(predict_->getRequiredFields() == Step::kUvwField);
  BOOST_TEST(predict_->getProvidedFields() == Step::kDataField);
}

BOOST_DATA_TEST_CASE(fields_add_subtract,
                     boost::unit_test::data::make({"add", "subtract"}),
                     operation) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("operation", operation);
  const OnePredict predict(parset, "", {});
  BOOST_TEST(predict.getRequiredFields() ==
             (Step::kDataField | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() == Step::kDataField);
}

BOOST_FIXTURE_TEST_CASE(fields_applycal, dp3::steps::test::H5ParmFixture) {
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  const OnePredict predict(parset, "", {});

  // OnePredict uses ApplyCal which has a OneApplyCal sub-step as next step.
  auto apply_cal = std::make_shared<ApplyCal>(parset, "applycal.", true);

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
  parset.add("applycal.parmdb", kParmDb);
  parset.add("applycal.correction", kSoltabName);
  parset.add("operation", operation);
  const OnePredict predict(parset, "", {});

  // When operation is "add" or "subtract", OnePredict only combines the
  // required fields of its ApplyCal sub-step.
  auto apply_cal = std::make_shared<ApplyCal>(parset, "applycal.", true);

  const dp3::common::Fields apply_cal_required =
      dp3::base::GetChainRequiredFields(apply_cal);

  BOOST_TEST(predict.getRequiredFields() ==
             (apply_cal_required | Step::kUvwField));
  BOOST_TEST(predict.getProvidedFields() == Step::kDataField);
}

/**
 * Create a buffer with artificial data values.
 * @note This function was originally copied from @ref
 * steps/test/unit/tIDGPredict.cc.
 */
static std::unique_ptr<dp3::base::DPBuffer> CreateBuffer() {
  const std::array<std::size_t, 3> kShape{kNBaselines, kNChan, kNCorr};

  auto buffer = std::make_unique<dp3::base::DPBuffer>(kTime, kInterval);
  buffer->GetData().resize(kShape);
  buffer->GetWeights().resize(kShape);
  buffer->GetFlags().resize(kShape);
  buffer->GetUvw().resize({kNBaselines, 3});

  buffer->GetFlags().fill(false);
  buffer->GetWeights().fill(1.0);

  for (std::size_t baseline = 0; baseline < kNBaselines; ++baseline) {
    // Base value for this baseline.
    const float baseline_value = baseline * 100.0;

    float channel_value = baseline_value;  // Base value for a group of channels
    for (std::size_t channel = 0; channel < kNChan; ++channel) {
      for (unsigned int corr = 0; corr < kNCorr; ++corr) {
        buffer->GetData()(baseline, channel, corr) = channel_value + corr;
      }
      // For each channel, increase channel_value by 10.0.
      channel_value += 10.0;
    }
    buffer->GetUvw()(baseline, 0) = baseline_value + 0.0;
    buffer->GetUvw()(baseline, 1) = baseline_value + 1.0;
    buffer->GetUvw()(baseline, 2) = baseline_value + 2.0;
  }

  return buffer;
}

BOOST_FIXTURE_TEST_CASE(showTimings, OnePredictFixture) {
  dp3::common::NSTimer timer;
  {
    const dp3::common::NSTimer::StartStop scoped_timer(timer);
    predict_->process(CreateBuffer());
  }

  std::stringstream sstr;
  // Ensure the test doesn't depend on the system's locale settings.
  sstr.imbue(std::locale::classic());
  predict_->showTimings(sstr, timer.getElapsed());
  const std::string output = sstr.str();

  {
    // The output percentage is between "  0.x" and "100.x".
    // Percentages above the 100% aren't validated.
    const std::regex regex{
        R"(  (1[0-9]| [ 0-9])[0-9]\.[0-9]% \([ 0-9]{5} [m ]s\) OnePredict fixture.\n)"
        R"(          (1[0-9]| [ 0-9])[0-9]\.[0-9]% \([ 0-9]{5} [m ]s\) of it spent in predict\n)"
        R"(          (1[0-9]| [ 0-9])[0-9]\.[0-9]% \([ 0-9]{5} [m ]s\) of it spent in apply beam\n)"};
    BOOST_CHECK(std::regex_match(output.begin(), output.end(), regex));
  }
  {
    // At the moment no beam is applied, so the percentages are fixed.
    // TODO Add an additional test to test with a beam applied.
    const std::regex regex{
        R"(  (1[0-9]| [ 0-9])[0-9]\.[0-9]% \([ 0-9]{5} [m ]s\) OnePredict fixture.\n)"
        R"(          100\.0% \([ 0-9]{5} [m ]s\) of it spent in predict\n)"
        R"(            0\.0% \(    0 ms\) of it spent in apply beam\n)"};
    BOOST_CHECK(std::regex_match(output.begin(), output.end(), regex));
  }
}

BOOST_AUTO_TEST_CASE(outputmodelname) {
  std::unique_ptr<dp3::base::DPBuffer> input_buffer = CreateBuffer();
  std::string output_model_name = "a_model_name";

  // Predict visibilities to main data buffer, replacing the input visibilities.
  // Make step chain
  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  auto predict =
      std::make_shared<OnePredict>(parset, "", std::vector<std::string>());
  auto predict_result = std::make_shared<dp3::steps::ResultStep>();
  predict->setNextStep(predict_result);
  predict->updateInfo(OnePredictFixture::MakeInfo());

  // Process and verify
  predict->process(std::make_unique<dp3::base::DPBuffer>(*input_buffer));
  std::unique_ptr<dp3::base::DPBuffer> result_main = predict_result->take();

  BOOST_CHECK(!(result_main->HasData(output_model_name)));
  BOOST_CHECK(!xt::allclose(result_main->GetData(), input_buffer->GetData()));

  // Predict visibilities to an extra data buffer in the output DPBuffer.
  // Make step chain (with extra parset pair)
  parset.add("outputmodelname", output_model_name);
  predict =
      std::make_shared<OnePredict>(parset, "", std::vector<std::string>());
  predict_result = std::make_shared<dp3::steps::ResultStep>();
  predict->setNextStep(predict_result);
  predict->updateInfo(OnePredictFixture::MakeInfo());

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

BOOST_AUTO_TEST_CASE(reusebeammodel) {
  const everybeam::vector3r_t kExpectedDirection{
      0.46738011434338467, -0.71981911052266068, 0.5132409539804188};

  dp3::common::ParameterSet parset;
  parset.add("sourcedb", dp3::steps::test::kPredictSkyModel);
  // Only use the first source in this test, which disables computing sources
  // in parallel. It also allows using MockTelescope / MockPointResponse, which
  // expect one MockPointResponse::Response call for each frequency.
  parset.add("sources", "[0002.2+3139]");
  parset.add("usebeammodel", "true");
  parset.add("reusebeammodel", "true");
  parset.add("elementmodel", "something_invalid");
  parset.add("usechannelfreq", "also_invalid");
  parset.add("coefficients_path", "/invalid/path");

  OnePredict predict(parset, "", {});
  const std::vector<std::string> kExpectedUnusedKeys = {
      "coefficients_path", "elementmodel", "usechannelfreq"};
  const std::vector<std::string> unused_keys = parset.unusedKeys();
  BOOST_CHECK_EQUAL_COLLECTIONS(unused_keys.begin(), unused_keys.end(),
                                kExpectedUnusedKeys.begin(),
                                kExpectedUnusedKeys.end());

  predict.setNextStep(std::make_shared<dp3::steps::NullStep>());

  dp3::base::DPInfo info = OnePredictFixture::MakeInfo();
  info.SetTelescope(std::make_shared<dp3::test::MockTelescope>(
      kChannelFrequencies, kExpectedDirection));
  predict.updateInfo(info);
  BOOST_CHECK(predict.getInfoOut().HasTelescope());
  BOOST_CHECK(&predict.getInfoOut().GetTelescope() == &info.GetTelescope());

  // OnePredict::process should call MockPointResponse::Response for each
  // frequency in kChannelFrequencies.
  predict.process(CreateBuffer());
}

BOOST_AUTO_TEST_SUITE_END()
