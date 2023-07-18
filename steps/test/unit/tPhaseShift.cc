// tPhaseShift.cc: Test program for class PhaseShift
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../PhaseShift.h"

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::PhaseShift;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(phaseshift)

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
    info().setTimes(0.0, (ntime - 1) * 10.0, 10.0);
    casacore::MDirection phaseCenter(casacore::Quantity(45, "deg"),
                                     casacore::Quantity(30, "deg"),
                                     casacore::MDirection::J2000);
    info().setArrayInformation(casacore::MPosition(), phaseCenter, phaseCenter,
                               phaseCenter);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 100000.0);
    std::vector<double> chanFreqs;
    for (std::size_t i = 0; i < nchan; i++) {
      chanFreqs.push_back(1050000.0 + i * 100000.0);
    }
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
    // Fill the baseline stations.
    // Determine nr of stations using:  na*(na+1)/2 = nbl
    // If many baselines, divide into groups of 6 to test if
    // PhaseShift disentangles it correctly.
    std::size_t nant = std::size_t(-0.5 + sqrt(0.25 + 2 * nbl));
    if (nant * (nant + 1) / 2 < nbl) ++nant;
    int grpszant = 3;
    std::size_t grpszbl = grpszant * (grpszant + 1) / 2;
    if (nbl > grpszbl) {
      nant = grpszant * (nbl + grpszbl - 1) / grpszbl;
    } else {
      grpszant = nant;
      grpszbl = nbl;
    }
    vector<int> ant1(nbl);
    vector<int> ant2(nbl);
    int st1 = 0;
    int st2 = 0;
    int lastant = grpszant;
    for (std::size_t i = 0; i < nbl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (i % grpszbl == grpszbl - 1) {
        st1 = lastant;
        st2 = lastant;
        lastant += grpszant;
      } else {
        if (++st2 == lastant) {
          st2 = ++st1;
        }
      }
    }
    vector<string> antNames(nant);
    vector<casacore::MPosition> antPos(nant);
    vector<double> antDiam(nant, 70.0);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    std::array<size_t, 2> stat_uvw_shape{itsNBl, 3};
    itsStatUVW.resize(stat_uvw_shape);
    for (std::size_t i = 0; i < nant; ++i) {
      itsStatUVW(i, 0) = 0.01 + i * 0.02;
      itsStatUVW(i, 1) = 0.05 + i * 0.03;
      itsStatUVW(i, 2) = 0.015 + i * 0.025;
    }
  }

  void fillUVW(xt::xtensor<double, 2>& uvw, int count) {
    for (std::size_t i = 0; i < itsNBl; ++i) {
      uvw(i, 0) = (itsStatUVW(getInfo().getAnt2()[i], 0) + count * 0.002 -
                   (itsStatUVW(getInfo().getAnt1()[i], 0) + count * 0.002));
      uvw(i, 1) = (itsStatUVW(getInfo().getAnt2()[i], 1) + count * 0.004 -
                   (itsStatUVW(getInfo().getAnt1()[i], 1) + count * 0.004));
      uvw(i, 2) = (itsStatUVW(getInfo().getAnt2()[i], 2) + count * 0.006 -
                   (itsStatUVW(getInfo().getAnt1()[i], 2) + count * 0.006));
    }
  }

 private:
  bool process(std::unique_ptr<DPBuffer>) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    std::array<size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};
    auto buffer = std::make_unique<DPBuffer>();
    buffer->setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buffer->GetData().resize(data_shape);
    auto& data = buffer->GetData();
    for (std::size_t i = 0; i < data.size(); ++i) {
      data.data()[i] = std::complex<float>(
          i + itsCount * 10,
          static_cast<int>(i) - 1000 + static_cast<int>(itsCount) * 6);
    }
    buffer->GetWeights().resize(data_shape);
    buffer->GetWeights().fill(1.0);
    buffer->GetFlags().resize(data_shape);
    buffer->GetFlags().fill(itsFlag);
    std::array<size_t, 2> uvw_shape{itsNBl, 3};
    xt::xtensor<double, 2> uvw(uvw_shape);
    fillUVW(uvw, itsCount);
    buffer->ResizeUvw(itsNBl);
    buffer->GetUvw() = uvw;
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
  xt::xtensor<double, 2> itsStatUVW;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(TestInput& input, std::size_t ntime, std::size_t nbl,
             std::size_t nchan, std::size_t ncorr, bool flag)
      : itsInput(&input),
        itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buf) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    std::array<size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> expected_result(data_shape);
    for (std::size_t i = 0; i < expected_result.size(); ++i) {
      expected_result.data()[i] = std::complex<float>(
          i + itsCount * 10,
          static_cast<int>(i) - 1000 + static_cast<int>(itsCount) * 6);
    }
    std::array<size_t, 2> uvw_shape{itsNBl, 3};
    xt::xtensor<double, 2> uvw(uvw_shape);
    itsInput->fillUVW(uvw, itsCount);
    // Check the expected result against the actual result.
    BOOST_CHECK(xt::allclose(buf->GetData(), expected_result, 1.0e-7));
    BOOST_CHECK(xt::all(xt::equal(buf->GetFlags(), itsFlag)));
    BOOST_CHECK_CLOSE_FRACTION(buf->getTime(), 2.0 + 5 * itsCount, 1.0e-6);
    BOOST_CHECK(xt::allclose(buf->GetUvw(), uvw, 1.0e-7));
    ++itsCount;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    casacore::MVDirection dir = infoIn.phaseCenter().getValue();
    BOOST_CHECK_CLOSE_FRACTION(dir.getLong("deg").getValue(), 45.0, 1.0e-6);
    BOOST_CHECK_CLOSE_FRACTION(dir.getLat("deg").getValue(), 30.0, 1.0e-6);
  }

  TestInput* itsInput;
  std::size_t itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput1 : public dp3::steps::test::ThrowStep {
 public:
  TestOutput1(TestInput& input, std::size_t ntime, std::size_t nbl,
              std::size_t nchan, std::size_t ncorr, bool flag)
      : itsInput(&input),
        itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buf) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    std::array<size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> expected_result(data_shape);
    for (std::size_t i = 0; i < expected_result.size(); ++i) {
      expected_result.data()[i] = std::complex<float>(
          i + itsCount * 10,
          static_cast<int>(i) - 1000 + static_cast<int>(itsCount) * 6);
    }
    std::array<size_t, 2> uvw_shape{itsNBl, 3};
    xt::xtensor<double, 2> uvw(uvw_shape);
    itsInput->fillUVW(uvw, itsCount);
    // Check the expected result against the actual result.
    BOOST_CHECK(!xt::allclose(buf->GetData(), expected_result));
    BOOST_CHECK(!xt::all(xt::equal(buf->GetData(), expected_result)));
    BOOST_CHECK(xt::all(xt::equal(buf->GetFlags(), itsFlag)));
    BOOST_CHECK_CLOSE_FRACTION(buf->getTime(), 2. + 5 * itsCount, 1.0e-5);
    BOOST_CHECK(!xt::allclose(buf->GetUvw(), uvw));
    ++itsCount;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    casacore::MVDirection dir = infoIn.phaseCenter().getValue();
    BOOST_CHECK_CLOSE_FRACTION(dir.getLong("deg").getValue(), 50.0, 1.0e-5);
    BOOST_CHECK_CLOSE_FRACTION(dir.getLat("deg").getValue(), 35.0, 1.0e-5);
  }

  TestInput* itsInput;
  std::size_t itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Test with a shift to the original center.
void test1(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  auto in = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  // Keep phase center the same to be able to check if data are correct.
  parset.add("phasecenter", "[45deg, 30deg]");
  auto phase_shift = std::make_shared<PhaseShift>(parset, "");
  auto out = std::make_shared<TestOutput>(*in, ntime, nbl, nchan, ncorr, flag);
  dp3::steps::test::Execute({in, phase_shift, out});
}

// Test with a shift to another and then to the original phase center.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  auto in = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  // First shift to another center, then back to original.
  ParameterSet parset1;
  parset1.add("phasecenter", "[50deg, 35deg]");
  ParameterSet parset2;
  parset2.add("phasecenter", "[]");
  auto phase_shift1 = std::make_shared<PhaseShift>(parset1, "");
  auto out1 =
      std::make_shared<TestOutput1>(*in, ntime, nbl, nchan, ncorr, flag);
  auto phase_shift2 = std::make_shared<PhaseShift>(parset2, "");
  auto out2 = std::make_shared<TestOutput>(*in, ntime, nbl, nchan, ncorr, flag);
  dp3::steps::test::Execute({in, phase_shift1, out1, phase_shift2, out2});
}

BOOST_AUTO_TEST_CASE(test1a) { test1(10, 3, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test1b) { test1(10, 10, 30, 1, true); }

BOOST_AUTO_TEST_CASE(test2a) { test2(10, 6, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test2b) { test2(10, 6, 30, 1, true); }

BOOST_AUTO_TEST_SUITE_END()
