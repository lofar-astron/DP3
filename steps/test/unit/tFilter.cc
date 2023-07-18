// tFilter.cc: Test program for class Filter
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Filter.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>

#include <boost/test/unit_test.hpp>

#include <xtensor/xio.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xview.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Filter;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(filter)

BOOST_AUTO_TEST_CASE(fields_default) {
  const ParameterSet parset;
  const Filter filter(parset, "");
  BOOST_CHECK_EQUAL(filter.getProvidedFields(), dp3::common::Fields());
  BOOST_CHECK_EQUAL(filter.getRequiredFields(), dp3::common::Fields());
}

BOOST_AUTO_TEST_CASE(fields_channel_selection) {
  // Test fields when enabling channel filtering.
  ParameterSet parset;
  parset.add("startchan", "4");
  parset.add("nchan", "2");
  const Filter channel_only(parset, "");
  const dp3::common::Fields kExpectedFields =
      Step::kDataField | Step::kFlagsField | Step::kWeightsField;
  BOOST_CHECK_EQUAL(channel_only.getRequiredFields(), kExpectedFields);
  BOOST_CHECK_EQUAL(channel_only.getProvidedFields(), kExpectedFields);

  ParameterSet parset_remove_ant = parset;
  parset_remove_ant.add("remove", "true");
  const Filter remove_ant(parset_remove_ant, "");
  const dp3::common::Fields kExpectedFieldsRemove =
      kExpectedFields | Step::kUvwField;
  BOOST_CHECK_EQUAL(remove_ant.getRequiredFields(), kExpectedFieldsRemove);
  BOOST_CHECK_EQUAL(remove_ant.getProvidedFields(), kExpectedFieldsRemove);
}

BOOST_AUTO_TEST_CASE(fields_baseline_selection) {
  ParameterSet parset;
  parset.add("baseline", "42");
  const Filter baseline_only(parset, "");

  const dp3::common::Fields kExpectedFields =
      Step::kDataField | Step::kFlagsField | Step::kWeightsField |
      Step::kUvwField;
  BOOST_CHECK_EQUAL(baseline_only.getRequiredFields(), kExpectedFields);
  BOOST_CHECK_EQUAL(baseline_only.getProvidedFields(), kExpectedFields);

  parset.add("startchan", "4");
  parset.add("nchan", "2");
  const Filter baseline_and_channel(parset, "");
  BOOST_CHECK_EQUAL(baseline_and_channel.getRequiredFields(), kExpectedFields);
  BOOST_CHECK_EQUAL(baseline_and_channel.getProvidedFields(), kExpectedFields);
}

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(std::size_t ntime, std::size_t nbl, std::size_t nchan,
            std::size_t ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {
    info() = DPInfo(ncorr, nchan);
    info().setTimes(100.0, 100.0 + (ntime - 1) * 5.0, 5.0);
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(nbl);
    std::vector<int> ant2(nbl);
    int st1 = 0;
    int st2 = 0;
    for (std::size_t i = 0; i < nbl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antNames{"rs01.s01", "rs02.s01", "cs01.s01",
                                      "cs01.s02"};
    // Define their positions (more or less WSRT RT0-3).
    std::vector<casacore::MPosition> antPos(4);
    casacore::Vector<double> vals(3);
    vals[0] = 3828763;
    vals[1] = 442449;
    vals[2] = 5064923;
    antPos[0] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828746;
    vals[1] = 442592;
    vals[2] = 5064924;
    antPos[1] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828729;
    vals[1] = 442735;
    vals[2] = 5064925;
    antPos[2] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    vals[0] = 3828713;
    vals[1] = 442878;
    vals[2] = 5064926;
    antPos[3] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
    std::vector<double> antDiam(4, 70.0);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 100000.0);
    std::vector<double> chanFreqs;
    for (std::size_t i = 0; i < nchan; i++) {
      chanFreqs.push_back(1050000.0 + i * 100000.0);
    }
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
  }

 private:
  bool process(std::unique_ptr<DPBuffer>) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    auto buffer = std::make_unique<DPBuffer>();
    buffer->setTime(itsCount * 5 + 2);
    buffer->setExposure(0.1 * (itsCount + 1));
    std::array<size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};
    buffer->GetData().resize(data_shape);
    for (std::size_t i = 0; i < buffer->GetData().size(); ++i) {
      buffer->GetData().data()[i] = std::complex<float>(
          i + itsCount * 10.0,
          static_cast<int>(i) - 1000 + static_cast<int>(itsCount) * 6.0);
    }
    buffer->GetWeights().resize(data_shape);
    const float last_weight = 0.5f + (itsNCorr * itsNChan * itsNBl * 0.01f);
    buffer->GetWeights() =
        xt::arange<float>(0.5, last_weight, 0.01).reshape(data_shape);
    buffer->ResizeFlags(data_shape);
    buffer->GetFlags().fill(itsFlag);
    // Set every third channel for every fourth baseline to the opposite flag to
    // make incorrect channel or baseline removal detectable.
    xt::strided_view(buffer->GetFlags(),
                     {xt::range(0, itsNBl, 4), xt::range(0, itsNChan, 3),
                      xt::all()}) = !itsFlag;
    buffer->ResizeUvw(itsNBl);
    buffer->GetUvw() =
        xt::arange(itsCount * 100.0, (itsNBl * 3.0) + (itsCount * 100.0), 1.0)
            .reshape({itsNBl, 3});
    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  std::size_t itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of averaging TestInput.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(std::size_t ntime, std::size_t nbl, std::size_t nchan,
             std::size_t ncorr, std::size_t nblout, std::size_t stchan,
             std::size_t nchanOut, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsNBlOut(nblout),
        itsStChan(stchan),
        itsNChanOut(nchanOut),
        itsFlag(flag) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Fill expected results to match input
    std::array<size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> data(data_shape);
    for (std::size_t i = 0; i < data.size(); ++i) {
      data.data()[i] = std::complex<float>(
          i + itsCount * 10.0,
          static_cast<int>(i) - 1000 + static_cast<int>(itsCount) * 6.0);
    }
    const float last_weight = 0.5f + (itsNCorr * itsNChan * itsNBl * 0.01f);
    xt::xtensor<float, 3> weights =
        xt::arange<float>(0.5, last_weight, 0.01).reshape(data_shape);
    xt::xtensor<bool, 3> flags(data_shape);
    flags.fill(itsFlag);
    // Set every third channel for every fourth baseline to the opposite flag to
    // make incorrect channel or baseline removal detectable.
    xt::strided_view(flags, {xt::range(0, itsNBl, 4), xt::range(0, itsNChan, 3),
                             xt::all()}) = !itsFlag;
    xt::xtensor<double, 2> uvw =
        xt::arange(itsCount * 100.0, (itsNBl * 3.0) + (itsCount * 100.0), 1.0)
            .reshape({itsNBl, 3});
    BOOST_TEST(buffer->GetData() ==
                   xt::view(data, xt::range(0, itsNBlOut),
                            xt::range(itsStChan, itsStChan + itsNChanOut),
                            xt::range(0, itsNCorr)),
               boost::test_tools::per_element());
    BOOST_TEST(buffer->GetFlags() ==
                   xt::view(flags, xt::range(0, itsNBlOut),
                            xt::range(itsStChan, itsStChan + itsNChanOut),
                            xt::range(0, itsNCorr)),
               boost::test_tools::per_element());
    BOOST_TEST(buffer->GetWeights() ==
                   xt::view(weights, xt::range(0, itsNBlOut),
                            xt::range(itsStChan, itsStChan + itsNChanOut),
                            xt::range(0, itsNCorr)),
               boost::test_tools::per_element());
    BOOST_TEST(
        buffer->GetUvw() == xt::view(uvw, xt::range(0, itsNBlOut), xt::all()),
        boost::test_tools::per_element());
    BOOST_CHECK_CLOSE(buffer->getTime(), itsCount * 5.0 + 2, 1.0e-3);
    BOOST_CHECK_CLOSE(buffer->getExposure(), 0.1 * (itsCount + 1), 1.0e-3);
    ++itsCount;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(itsNChan, info.origNChan());
    BOOST_CHECK_EQUAL(itsNChanOut, info.nchan());
    BOOST_CHECK_EQUAL(itsNBlOut, info.nbaselines());
    BOOST_CHECK_EQUAL(itsNTime, info.ntime());
    BOOST_CHECK_EQUAL(5.0, info.timeInterval());
    BOOST_CHECK_EQUAL(1, info.nchanAvg());
    BOOST_CHECK_EQUAL(1, info.ntimeAvg());
  }

  std::size_t itsCount, itsNTime, itsNBl, itsNChan, itsNCorr, itsNBlOut,
      itsStChan, itsNChanOut;
  bool itsFlag;
};

void TestChannelsOnly(std::size_t ntime, std::size_t nbl, std::size_t nchan,
                      std::size_t ncorr, std::size_t startchan,
                      std::size_t nchanout, bool flag) {
  auto in = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("startchan", std::to_string(startchan));
  parset.add("nchan", std::to_string(nchanout) + "+nchan-nchan");
  auto filter = std::make_shared<Filter>(parset, "");
  auto out = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, nbl,
                                          startchan, nchanout, flag);
  dp3::steps::test::Execute({in, filter, out});
}

void TestChannelsAndBaselines(std::size_t ntime, std::size_t nbl,
                              std::size_t nchan, std::size_t ncorr,
                              std::size_t startchan, std::size_t nchanout,
                              bool flag) {
  BOOST_CHECK(nbl <=
              4);  // otherwise baseline selection removes more than the first

  auto in = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("startchan", std::to_string(startchan) + "+nchan-nchan");
  parset.add("nchan", std::to_string(nchanout));
  // This removes the first baseline.
  parset.add("baseline", "[[rs01.s01,rs*]]");
  auto filter = std::make_shared<Filter>(parset, "");
  auto out = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, 2,
                                          startchan, nchanout, flag);
  dp3::steps::test::Execute({in, filter, out});
}

BOOST_AUTO_TEST_CASE(filter_channels_1) {
  TestChannelsOnly(10, 3, 32, 4, 2, 24, false);
}

BOOST_AUTO_TEST_CASE(filter_channels_2) {
  TestChannelsOnly(10, 10, 30, 1, 3, 3, true);
}

BOOST_AUTO_TEST_CASE(filter_channels_3) {
  TestChannelsOnly(10, 10, 1, 4, 0, 1, true);
}

BOOST_AUTO_TEST_CASE(channels_and_baselines) {
  TestChannelsAndBaselines(10, 4, 32, 4, 2, 24, false);
}

BOOST_AUTO_TEST_SUITE_END()
