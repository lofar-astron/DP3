// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mock/MockStep.h"
#include "../../MSBDAReader.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::MSBDAReader;

namespace {
const std::string kPrefix = "";
const dp3::common::ParameterSet kParset;
}  // namespace

BOOST_AUTO_TEST_SUITE(msbdareader)

BOOST_AUTO_TEST_CASE(constructor) {
  const casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  MSBDAReader reader(ms, kParset, kPrefix);
  BOOST_TEST(reader.table().isRootTable());
  BOOST_TEST(reader.table().isSameRoot(ms));

  BOOST_TEST((reader.outputs() == dp3::steps::Step::MsType::kBda));
}

BOOST_AUTO_TEST_CASE(set_info) {
  casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  DPInfo info;
  MSBDAReader reader(ms, kParset, kPrefix);

  reader.setInfo(info);
  const double start_time =
      4472025725;  // conversion of start time for tNDPPP_bda_tmp.MS
                   // (2000/08/03 13h22m05.000) into seconds

  info = reader.getInfo();
  BOOST_TEST(reader.spectralWindow() == 0U);
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
  DPInfo info;
  MSBDAReader reader(ms, kParset, kPrefix);
  reader.setReadVisData(true);
  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);
  reader.setInfo(info);
  dp3::base::DPBuffer buf;

  reader.process(buf);
  reader.finish();

  auto kExpectedData = std::complex<float>(2.75794, 0.899097);
  double kExpectedUVW[] = {-0.125726, -766.917, 397.793};
  auto kExpectedWeights = std::vector<float>(
      reader.getInfo().ncorr() * reader.getInfo().nchan(), 1);

  BOOST_TEST(mock_step->FinishCount() == std::size_t(1));
  BOOST_TEST(mock_step->TotalRowCount() == std::size_t(6));

  BOOST_REQUIRE(!mock_step->GetBdaBuffers().empty());
  auto rows = mock_step->GetBdaBuffers().front()->GetRows();
  BOOST_REQUIRE(!rows.empty());
  BOOST_TEST(rows[0].data->imag() == kExpectedData.imag());
  BOOST_TEST(rows[0].data->real() == kExpectedData.real());
  BOOST_TEST(rows[0].uvw == kExpectedUVW);
  for (unsigned i = 0; i < kExpectedWeights.size(); ++i) {
    BOOST_TEST(*(rows[0].weights + i) == kExpectedWeights[i]);
  }

  // Check that the print methods do not throw any error
  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(reader.show(nullout));
  BOOST_CHECK_NO_THROW(reader.showCounts(nullout));
  BOOST_CHECK_NO_THROW(reader.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_CASE(process_nan) {
  casacore::MeasurementSet ms("tNDPPP_bda_tmp.MS");
  DPInfo info;
  MSBDAReader reader(ms, kParset, kPrefix);
  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);
  reader.setInfo(info);
  dp3::base::DPBuffer buf;

  reader.setReadVisData(false);
  reader.process(buf);
  reader.finish();

  std::complex<float>* data = mock_step->GetBdaBuffers()[0]->GetRows()[0].data;
  BOOST_TEST(std::isnan(data->imag()));
  BOOST_TEST(std::isnan(data->real()));
}

BOOST_AUTO_TEST_SUITE_END()
