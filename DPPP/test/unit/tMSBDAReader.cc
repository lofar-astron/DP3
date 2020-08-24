// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "mock/MockStep.h"
#include "../../MSBDAReader.h"
#include "../../DPInfo.h"
#include "../../../Common/ParameterSet.h"

#include <boost/test/unit_test.hpp>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"

using DP3::ParameterSet;
using DP3::DPPP::DPInfo;
using DP3::DPPP::MSBDAReader;

namespace {
const std::string prefix = "";
}  // namespace

BOOST_AUTO_TEST_SUITE(dpinfo)

BOOST_AUTO_TEST_CASE(flow) {
  DPInfo info;
  ParameterSet parset;
  MSBDAReader reader("does_not_exist.MS", parset, prefix);

  BOOST_TEST(reader.outputs() == DP3::DPPP::DPStep::MSType::BDA);
}

BOOST_AUTO_TEST_CASE(empty_input) {
  DPInfo info;
  ParameterSet parset;
  MSBDAReader reader("does_not_exist.MS", parset, prefix);
  reader.setReadVisData(true);

  BOOST_CHECK_THROW(reader.setInfo(info), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(wrong_input_format) {
  DPInfo info;
  ParameterSet parset;
  MSBDAReader reader("tNDPPP_tmp.MS", parset, prefix);

  BOOST_CHECK_THROW(reader.setInfo(info), std::domain_error);
}

BOOST_AUTO_TEST_CASE(set_info) {
  const std::string kMSName = "tNDPPP_bda_tmp.MS";
  DPInfo info;
  ParameterSet parset;
  MSBDAReader reader(kMSName, parset, prefix);

  reader.setInfo(info);

  info = reader.getInfo();
  BOOST_TEST(reader.table().tableName().compare(
                 reader.table().tableName().length() - kMSName.length(),
                 kMSName.length(), kMSName) == 0);
  BOOST_TEST(reader.spectralWindow() == -1U);
  BOOST_TEST(info.nchan() == 16U);
  BOOST_TEST(info.ncorr() == 4U);
  BOOST_TEST(info.ntime() == 1U);
  BOOST_TEST(info.channelsAreRegular());
  BOOST_TEST(info.hasBDAChannels());
  BOOST_TEST(info.nbaselines() == 6U);
  BOOST_TEST(info.startchan() == 0U);
  BOOST_TEST(info.startTime() == 0U);
  BOOST_TEST(info.timeInterval() == 30U);
}

BOOST_AUTO_TEST_CASE(process, *boost::unit_test::tolerance(0.001) *
                                  boost::unit_test::tolerance(0.0001f)) {
  DPInfo info;
  ParameterSet parset;
  parset.add("memoryperc", "100");
  parset.add("memory_max", "1");
  MSBDAReader reader("tNDPPP_bda_tmp.MS", parset, prefix);
  auto mock_step = std::make_shared<DP3::DPPP::MockStep>();
  reader.setNextStep(mock_step);
  reader.setInfo(info);
  DP3::DPPP::DPBuffer buf;

  reader.process(buf);
  reader.finish();

  auto kExpectedData = std::complex<float>(2.75794, 0.899097);
  double kExpectedUVW[] = {-0.125726, -766.917, 397.793};
  auto kExpectedWeights = std::vector<float>(
      reader.getInfo().ncorr() * reader.getInfo().nchan(), 1);

  mock_step->CheckFinishCount(1);
  BOOST_TEST(mock_step->GetBdaBuffers().size() == 1U);
  BOOST_TEST(mock_step->GetBdaBuffers()[0]->GetRows().size() == 120U);
  BOOST_TEST(mock_step->GetBdaBuffers()[0]->GetRows()[0].data->imag() ==
             kExpectedData.imag());
  BOOST_TEST(mock_step->GetBdaBuffers()[0]->GetRows()[0].data->real() ==
             kExpectedData.real());
  BOOST_TEST(mock_step->GetBdaBuffers()[0]->GetRows()[0].uvw == kExpectedUVW);
  for (unsigned i = 0; i < kExpectedWeights.size(); ++i) {
    BOOST_TEST(*(mock_step->GetBdaBuffers()[0]->GetRows()[0].weights + i) ==
               kExpectedWeights[i]);
  }

  // Check that the print methods do not throw any error
  std::ostream nullout(nullptr);
  BOOST_CHECK_NO_THROW(reader.show(nullout));
  BOOST_CHECK_NO_THROW(reader.showCounts(nullout));
  BOOST_CHECK_NO_THROW(reader.showTimings(nullout, 10));
}

BOOST_AUTO_TEST_SUITE_END()
