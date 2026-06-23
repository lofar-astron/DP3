// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/ApplyBeam.h"

#include <boost/test/unit_test.hpp>
#include <EveryBeam/load.h>
#include <xtensor/misc/xcomplex.hpp>
#include <xtensor/views/xview.hpp>

#include "mock/ThrowStep.h"
#include "common/ParameterSet.h"
#include "steps/ResultStep.h"
#include "steps/test/unit/mock/MockTelescope.h"

using dp3::steps::ApplyBeam;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(apply_beam)

namespace {
const double kTime{4.87128e+09};
const double kExposure{10.0139};
const std::vector<double> kFrequencies{1.34288e+08, 1.34312e+08};
}  // namespace

// Create and populate a buffer with the bare minimum values required to run
// ApplyBeam successfully
std::unique_ptr<dp3::base::DPBuffer> CreateBuffer(
    const xt::xtensor<std::complex<float>, 3>& data,
    const xt::xtensor<float, 3>& weights) {
  std::unique_ptr<dp3::base::DPBuffer> buffer =
      std::make_unique<dp3::base::DPBuffer>(kTime, kExposure);

  buffer->GetData().resize(data.shape());
  buffer->GetWeights().resize(weights.shape());

  buffer->GetData() = data;
  buffer->GetWeights() = weights;

  return buffer;
}

// Create and populate an info object with the bare minimum values required to
// run ApplyBeam successfully
dp3::base::DPInfo MakeInfo() {
  dp3::base::DPInfo info(4, 8, "HBA_DUAL_INNER");
  info.setMsName("tNDPPP-generic.MS");
  info.setChannels(std::vector<double>(kFrequencies), {24414.1, 24414.1},
                   {24414.1, 24414.1}, {24414.1, 24414.1}, 1.34375e+08, 0);

  info.setTimes(3.0, 3.0, 5.0);

  const std::vector<int> antenna_1{0};
  const std::vector<int> antenna_2{0};

  const std::vector<std::string> antenna_names{"CS001HBA0", "CS002HBA0"};

  std::vector<casacore::MPosition> antenna_position(2);
  casacore::Vector<double> vals(3);
  vals[0] = 3828763;
  vals[1] = 442449;
  vals[2] = 5064923;
  antenna_position[0] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vals[0] = 3828746;
  vals[1] = 442592;
  vals[2] = 5064924;
  antenna_position[1] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);

  const std::vector<double> antenna_diameter(2.0, 70.0);
  info.setAntennas(antenna_names, antenna_diameter, antenna_position, antenna_1,
                   antenna_2);
  return info;
}

// Perform a very small/minimal run of the ArrayBeam step and compare against
// expected result to help catch any regressions
BOOST_AUTO_TEST_CASE(test_step_basic) {
  const xt::xtensor<std::complex<float>, 3> kTestData = {
      {{{0.000105909, -6.15896e-06},
        {-5.75999e-05, -0.000169081},
        {6.09389e-05, -0.000144006},
        {7.58708e-05, -0.000118892}},
       {{5.04177e-05, -0.000176506},
        {0.000159816, 3.01513e-06},
        {8.10528e-05, 9.37562e-05},
        {0.000141897, -7.32e-05}}},
      {{{0.000170118, 5.01791e-05},
        {6.19825e-05, 3.1506e-05},
        {0.000152869, -8.22876e-05},
        {3.44475e-05, -0.000179837}},
       {{0.000143976, -0.000107417},
        {0.000130452, -3.40554e-05},
        {0.000100048, 7.11484e-05},
        {-0.000127001, -6.57001e-05}}},
      {{{6.60433e-05, -9.32823e-05},
        {6.31116e-05, 0.000100894},
        {-9.42722e-05, -0.000188316},
        {-5.7855e-05, -0.000132138}},
       {{0.000108384, 8.16866e-05},
        {-0.00012894, -0.000180018},
        {0.000106929, -2.96774e-05},
        {0.000101331, -1.57673e-05}}}};
  const xt::xtensor<float, 3> kTestWeights = {
      {{7.429210e+08, 1.792340e+09, 5.485560e+08, 1.411610e+09},
       {4.384190e+08, 5.486580e+08, 5.480330e+08, 4.262130e+08}},
      {{1.745190e+09, 8.292680e+08, 1.791950e+09, 1.608370e+09},
       {4.232110e+08, 4.369040e+08, 6.651330e+08, 1.709300e+09}},
      {{9.355210e+08, 1.963980e+09, 6.532890e+08, 5.810230e+08},
       {5.264440e+08, 1.566430e+09, 1.715470e+09, 1.069660e+09}}};
  const xt::xtensor<std::complex<float>, 3> kExpectedData = {
      {{{886.133, -2354.56},
        {771.412, -602.496},
        {-367.09, -843.605},
        {883.182, 710.44}},
       {{2249.85, -563.42},
        {199.536, 273.074},
        {926.346, -564.403},
        {-137.817, -1638.45}}},
      {{{0.000170118, 5.01791e-05},
        {6.19825e-05, 3.1506e-05},
        {0.000152869, -8.22876e-05},
        {3.44475e-05, -0.000179837}},
       {{0.000143976, -0.000107417},
        {0.000130452, -3.40554e-05},
        {0.000100048, 7.11484e-05},
        {-0.000127001, -6.57001e-05}}},
      {{{6.60433e-05, -9.32823e-05},
        {6.31116e-05, 0.000100894},
        {-9.42722e-05, -0.000188316},
        {-5.7855e-05, -0.000132138}},
       {{0.000108384, 8.16866e-05},
        {-0.00012894, -0.000180018},
        {0.000106929, -2.96774e-05},
        {0.000101331, -1.57673e-05}}}};
  const xt::xtensor<float, 3> kExpectedWeights = {
      {{7.42921e+08, 1.79234e+09, 5.48556e+08, 1.41161e+09},
       {4.38419e+08, 5.48658e+08, 5.48033e+08, 4.26213e+08}},
      {{1.74519e+09, 8.29268e+08, 1.79195e+09, 1.60837e+09},
       {4.23211e+08, 4.36904e+08, 6.65133e+08, 1.7093e+09}},
      {{9.35521e+08, 1.96398e+09, 6.53289e+08, 5.81023e+08},
       {5.26444e+08, 1.56643e+09, 1.71547e+09, 1.06966e+09}}};

  auto apply_beam =
      std::make_unique<ApplyBeam>(dp3::common::ParameterSet(), "", false);
  auto result = std::make_shared<dp3::steps::ResultStep>();
  apply_beam->setNextStep(result);
  apply_beam->updateInfo(MakeInfo());

  auto buffer = CreateBuffer(kTestData, kTestWeights);
  apply_beam->process(std::move(buffer));

  std::unique_ptr<dp3::base::DPBuffer> result_buffer = result->take();
  BOOST_TEST_REQUIRE(result_buffer);
  BOOST_TEST(result_buffer->GetTime() == kTime);
  BOOST_TEST(result_buffer->GetExposure() == kExposure);
  BOOST_TEST(xt::allclose(result_buffer->GetData(), kExpectedData));
  BOOST_TEST(result_buffer->GetWeights() == kExpectedWeights,
             boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(reusebeammodel) {
  const std::array<size_t, 3> shape = {3, 2, 4};

  const xt::xtensor<std::complex<float>, 3> kTestData(
      shape, std::complex<float>(1.0, 1.0));
  const xt::xtensor<float, 3> kWeights(shape, 1.0f);

  // The expected data is 0.25 for the first baseline, since
  // MockPointResponse returns a beam value of 2.0 and the input data is 1.0.
  // The other two baselines are unchanged. Since MakeInfo() only creates a
  // single baseline, ApplyBeam should not touch the extra data in the buffer.
  xt::xtensor<std::complex<float>, 3> expected_data = kTestData;
  xt::view(expected_data, 0, xt::all(), xt::all()) =
      std::complex<float>(0.25, 0.25);

  // Expected direction for MockPointResponse::Response() calls.
  const everybeam::vector3r_t kExpectedDirection{
      -0.0011852071523323795, 0.00053077365480819693, 0.99999915678131113};

  dp3::common::ParameterSet parset;
  parset.add("reusebeammodel", "true");
  parset.add("elementmodel", "something_invalid");
  parset.add("usechannelfreq", "also_invalid");
  parset.add("coefficients_path", "/invalid/path");

  ApplyBeam apply_beam(parset, "");
  const std::vector<std::string> kExpectedUnusedKeys = {
      "coefficients_path", "elementmodel", "usechannelfreq"};
  const std::vector<std::string> unused_keys = parset.unusedKeys();
  BOOST_CHECK_EQUAL_COLLECTIONS(unused_keys.begin(), unused_keys.end(),
                                kExpectedUnusedKeys.begin(),
                                kExpectedUnusedKeys.end());

  auto result = std::make_shared<dp3::steps::ResultStep>();
  apply_beam.setNextStep(result);

  dp3::base::DPInfo info = MakeInfo();
  info.SetTelescope(std::make_shared<dp3::test::MockTelescope>(
      kFrequencies, kExpectedDirection));
  apply_beam.updateInfo(info);
  BOOST_CHECK(apply_beam.getInfoOut().HasTelescope());
  BOOST_CHECK(&apply_beam.getInfoOut().GetTelescope() == &info.GetTelescope());

  apply_beam.process(CreateBuffer(kTestData, kWeights));
  std::unique_ptr<dp3::base::DPBuffer> result_buffer = result->take();
  BOOST_TEST_REQUIRE(result_buffer);
  BOOST_TEST(result_buffer->GetData() == expected_data,
             boost::test_tools::per_element());
  BOOST_TEST(result_buffer->GetWeights() == kWeights,
             boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(reusebeammodel_no_telescope) {
  dp3::common::ParameterSet parset;
  parset.add("reusebeammodel", "true");
  ApplyBeam apply_beam(parset, "");
  BOOST_CHECK_THROW(apply_beam.updateInfo(MakeInfo()), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(fields_defaults) {
  dp3::common::ParameterSet parset;
  const ApplyBeam apply_beam(parset, "");
  BOOST_TEST(apply_beam.getRequiredFields() == Step::kDataField);
  BOOST_TEST(apply_beam.getProvidedFields() == Step::kDataField);
}

BOOST_AUTO_TEST_CASE(fields_updateweights) {
  dp3::common::ParameterSet parset;
  parset.add("updateweights", "true");
  const ApplyBeam updates_weights(parset, "");

  const dp3::common::Fields kExpectedFields =
      Step::kDataField | Step::kWeightsField;
  BOOST_TEST(updates_weights.getRequiredFields() == kExpectedFields);
  BOOST_TEST(updates_weights.getProvidedFields() == kExpectedFields);
}

BOOST_AUTO_TEST_SUITE_END()
