// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSBDAReader.h"

#include <boost/test/unit_test.hpp>

#include "mock/MockStep.h"
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"

using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::MSBDAReader;

namespace {
const std::string kPrefix = "";
const dp3::common::ParameterSet kParset;

/// Runs the test in a separate directory containing tNDPPP_bda_tmp.MS.
class BdaMsFixture : public dp3::common::test::FixtureDirectory {
 public:
  BdaMsFixture() : FixtureDirectory() {
    ExtractResource("tNDPPP_bda.in_MS.tgz");
  }
};

}  // namespace

BOOST_FIXTURE_TEST_SUITE(msbdareader, BdaMsFixture)

BOOST_AUTO_TEST_CASE(constructor) {
  const casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  MSBDAReader reader(ms, kParset, kPrefix);
  BOOST_TEST(reader.table().isRootTable());
  BOOST_TEST(reader.table().isSameRoot(ms));

  BOOST_TEST((reader.outputs() == dp3::steps::Step::MsType::kBda));
}

BOOST_AUTO_TEST_CASE(set_info) {
  casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  MSBDAReader reader(ms, kParset, kPrefix);

  reader.setInfo(DPInfo());
  const double start_time =
      4472025725;  // conversion of start time for tNDPPP_bda_tmp.MS
                   // (2000/08/03 13h22m05.000) into seconds

  const DPInfo& info = reader.getInfoOut();
  BOOST_TEST(info.spectralWindow() == 0U);
  BOOST_TEST(info.nchan() == 16U);
  BOOST_TEST(info.ncorr() == 4U);
  // With BDA we approximate this amount of buffers to be streamed.
  BOOST_TEST(info.ntime() == reader.table().nrow() / info.nbaselines());
  BOOST_TEST(info.channelsAreRegular());
  BOOST_TEST(info.hasBDAChannels());
  BOOST_TEST(info.nbaselines() == 6U);
  BOOST_TEST(info.startchan() == 0U);
  BOOST_TEST(info.startTime() == start_time);
  BOOST_TEST(info.timeInterval() == 30U);
}

BOOST_AUTO_TEST_CASE(process, *boost::unit_test::tolerance(0.0001) *
                                  boost::unit_test::tolerance(0.0001f)) {
  casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  MSBDAReader reader(ms, kParset, kPrefix);
  reader.setFieldsToRead(dp3::steps::Step::kDataField |
                         dp3::steps::Step::kWeightsField);
  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  reader.updateInfo(DPInfo());
  // Call process function with null (DPBuffer) pointer to ensure it doesn't do
  // anything with the input.
  reader.process(std::unique_ptr<dp3::base::DPBuffer>());
  reader.finish();

  auto kExpectedData = std::complex<float>(2.75794, 0.899097);
  double kExpectedUVW[] = {-0.125726, -766.917, 397.793};
  auto kExpectedWeights = std::vector<float>(
      reader.getInfoOut().ncorr() * reader.getInfoOut().nchan(), 1.0);

  BOOST_TEST(mock_step->FinishCount() == std::size_t(1));
  BOOST_TEST(mock_step->TotalRowCount() == std::size_t(6));

  BOOST_REQUIRE(!mock_step->GetBdaBuffers().empty());
  const dp3::base::BdaBuffer& buffer = *mock_step->GetBdaBuffers().front();
  BOOST_REQUIRE(!buffer.GetRows().empty());
  BOOST_TEST(buffer.GetData(0)->imag() == kExpectedData.imag());
  BOOST_TEST(buffer.GetData(0)->real() == kExpectedData.real());
  BOOST_TEST(buffer.GetRows()[0].uvw == kExpectedUVW);
  for (std::size_t i = 0; i < kExpectedWeights.size(); ++i) {
    BOOST_TEST(buffer.GetWeights(0)[i] == kExpectedWeights[i]);
  }
}

BOOST_AUTO_TEST_CASE(process_no_fields_to_read) {
  casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  MSBDAReader reader(ms, kParset, kPrefix);
  BOOST_TEST(reader.getFieldsToRead() == dp3::common::Fields());
  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);
  reader.updateInfo(DPInfo());

  reader.process(std::unique_ptr<dp3::base::BdaBuffer>());
  reader.finish();

  BOOST_TEST_REQUIRE(mock_step->GetBdaBuffers().size() == 1);
  BOOST_TEST(!mock_step->GetBdaBuffers().front()->GetData());
  BOOST_TEST(!mock_step->GetBdaBuffers().front()->GetWeights());
  BOOST_TEST(!mock_step->GetBdaBuffers().front()->GetFlags());
}

BOOST_AUTO_TEST_CASE(show) {
  // Check that the show methods do not throw any error.
  const casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  const MSBDAReader reader(ms, kParset, kPrefix);

  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(reader.show(nullout));
  BOOST_CHECK_NO_THROW(reader.showCounts(nullout));
  BOOST_CHECK_NO_THROW(reader.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_SUITE_END()
