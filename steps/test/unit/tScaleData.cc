// tScaleData.cc: Test program for class ScaleData
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../ScaleData.h"
#include "../../InputStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::ScaleData;
using std::vector;

namespace {
const double kBaseFreq = 10.5;  // MHz
}

BOOST_AUTO_TEST_SUITE(scaledata)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// It can be used with different nr of times, channels, etc.
class TestInput : public InputStep {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {
    info().init(ncorr, 0, nchan, ntime, 0., 5., string(), string());
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(nbl);
    vector<int> ant2(nbl);
    for (int i = 0; i < nbl; ++i) {
      ant1[i] = i / 4;
      ant2[i] = i % 4;
    }
    vector<string> antNames{"rs01.s01", "rs02.s01", "cs01.s01", "cs01.s02"};
    // Define their positions (more or less WSRT RT0-3).
    vector<casacore::MPosition> antPos(4);
    vector<double> vals(3);
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
    vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 1e6);
    std::vector<double> chanFreqs;
    for (int i = 0; i < nchan; i++) {
      chanFreqs.push_back((kBaseFreq + i) * 1e6);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
  }

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 10 + itsCount * 6);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 0.5f, 0.01f);
    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime(itsCount * 30 + 4472025740.0);
    buf.setData(data);
    buf.setWeights(weights);
    buf.setUVW(uvw);
    casacore::Cube<bool> flags(data.shape());
    flags = false;
    buf.setFlags(flags);
    casacore::Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = false;
    buf.setFullResFlags(fullResFlags);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of TestInput run by TestScaling.
class TestOutput : public dp3::steps::Step {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  void addData(casacore::Cube<casacore::Complex>& to,
               const casacore::Cube<casacore::Complex>& from, int bl) {
    to += from(casacore::IPosition(3, 0, 0, bl),
               casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl));
  }
  void addConjData(casacore::Cube<casacore::Complex>& to,
                   const casacore::Cube<casacore::Complex>& from, int bl) {
    to +=
        conj(from(casacore::IPosition(3, 0, 0, bl),
                  casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl)));
  }
  virtual bool process(const DPBuffer& buf) {
    // Fill data and scale as needed.
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    casacore::Complex* dataPtr = data.data();
    int cnt = 0;
    for (int i = 0; i < itsNBl; ++i) {
      for (int j = 0; j < itsNChan; ++j) {
        double freq = kBaseFreq + j;
        double coeff1 = 2 + 0.5 * freq;
        double coeff2 = 3 + 2 * freq + 1 * freq * freq;
        // The first antenna uses coeff1, the others use coeff2.
        double sc1 = info().getAnt1()[i] == 0 ? coeff1 : coeff2;
        double sc2 = info().getAnt2()[i] == 0 ? coeff1 : coeff2;
        double scale = sqrt(sc1 * sc2);
        for (int k = 0; k < itsNCorr; ++k) {
          *dataPtr =
              casacore::Complex(cnt + itsCount * 10, cnt - 10 + itsCount * 6) *
              float(scale);
          ++dataPtr;
          ++cnt;
        }
      }
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 0.5f, 0.01f);
    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    BOOST_CHECK(allEQ(buf.getData(), data));
    BOOST_CHECK_EQUAL(buf.getFlags().shape(),
                      casacore::IPosition(3, itsNCorr, itsNChan, itsNBl));
    BOOST_CHECK(allEQ(buf.getFlags(), false));
    BOOST_CHECK(allEQ(buf.getWeights(), weights));
    BOOST_CHECK(allEQ(buf.getUVW(), uvw));
    BOOST_CHECK_EQUAL(buf.getFullResFlags().shape(),
                      casacore::IPosition(3, itsNChan, 1, itsNBl));
    BOOST_CHECK(allEQ(buf.getFullResFlags(), false));
    itsCount++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.nbaselines()), itsNBl);
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

void TestScaling(int ntime, int nbl, int nchan, int ncorr) {
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr);
  ParameterSet parset;
  parset.add("stations", "[rs01.s01, *]");
  parset.add("coeffs", "[[2,0.5],[3,2,1]]");
  parset.add("scalesize", "false");
  auto step2 = std::make_shared<ScaleData>(step1.get(), parset, "",
                                           ScaleData::MsType::kRegular);
  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_CASE(test_scale_data1) { TestScaling(2, 4, 4, 1); }

BOOST_AUTO_TEST_CASE(test_scale_data2) { TestScaling(10, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(test_scale_data3) { TestScaling(10, 12, 16, 2); }

BOOST_AUTO_TEST_SUITE_END()
