// tScaleData.cc: Test program for class ScaleData
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <complex>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include "../../ScaleData.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::ScaleData;
using std::vector;

namespace {
const double kBaseFreq = 10.5;  // MHz
}

BOOST_AUTO_TEST_SUITE(scaledata)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(std::size_t ntime, std::size_t nbl, std::size_t nchan,
            std::size_t ncorr)
      : count_(0),
        n_time_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        n_correlations_(ncorr) {
    GetWritableInfoOut() = DPInfo(ncorr, nchan);
    GetWritableInfoOut().setTimes(0.0, (ntime - 1) * 5.0, 5.0);
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(nbl);
    vector<int> ant2(nbl);
    for (std::size_t i = 0; i < nbl; ++i) {
      ant1[i] = i / 4;
      ant2[i] = i % 4;
    }
    vector<std::string> antNames{"rs01.s01", "rs02.s01", "cs01.s01",
                                 "cs01.s02"};
    // Define their positions (more or less WSRT RT0-3).
    vector<casacore::MPosition> antPos(4);
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
    vector<double> antDiam(4, 70.);
    GetWritableInfoOut().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 1e6);
    std::vector<double> chanFreqs;
    for (std::size_t i = 0; i < nchan; i++) {
      chanFreqs.push_back((kBaseFreq + i) * 1e6);
    }
    GetWritableInfoOut().setChannels(std::move(chanFreqs),
                                     std::move(chanWidth));
  }

 private:
  bool process(std::unique_ptr<dp3::base::DPBuffer> buffer) override {
    // Stop when all times are done.
    if (count_ == n_time_) {
      return false;
    }

    buffer->SetTime(count_ * 30 + 4472025740.0);
    const std::array shape{n_baselines_, n_channels_, n_correlations_};
    buffer->GetData().resize(shape);
    std::generate(buffer->GetData().begin(), buffer->GetData().end(),
                  [&, i = -1]() mutable {
                    ++i;
                    return std::complex<float>(i + count_ * 10,
                                               i - 10 + count_ * 6);
                  });
    buffer->GetWeights().resize(shape);
    std::generate(buffer->GetWeights().begin(), buffer->GetWeights().end(),
                  [result = 0.5f]() mutable {
                    result += 0.01f;
                    return result;
                  });

    buffer->GetFlags().resize(shape);
    buffer->GetFlags().fill(false);

    buffer->GetUvw().resize({n_baselines_, 3});
    for (std::size_t i = 0; i < n_baselines_; ++i) {
      buffer->GetUvw()(i, 0) = 1 + count_ + i;
      buffer->GetUvw()(i, 1) = 2 + count_ + i;
      buffer->GetUvw()(i, 2) = 3 + count_ + i;
    }

    getNextStep()->process(std::move(buffer));
    ++count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  int count_;
  int n_time_;
  std::size_t n_baselines_;
  std::size_t n_channels_;
  std::size_t n_correlations_;
};

// Class to check result of TestInput run by TestScaling.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(std::size_t ntime, std::size_t nbl, std::size_t nchan,
             std::size_t ncorr)
      : count_(0),
        n_time_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        n_correlations_(ncorr) {}

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
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Fill data and scale as needed.
    const std::array shape{n_baselines_, n_channels_, n_correlations_};
    xt::xtensor<std::complex<float>, 3> data{shape};

    int cnt = 0;
    for (std::size_t i = 0; i < n_baselines_; ++i) {
      for (std::size_t j = 0; j < n_channels_; ++j) {
        double freq = kBaseFreq + j;
        double coeff1 = 2 + 0.5 * freq;
        double coeff2 = 3 + 2 * freq + 1 * freq * freq;
        // The first antenna uses coeff1, the others use coeff2.
        double sc1 = getInfoOut().getAnt1()[i] == 0 ? coeff1 : coeff2;
        double sc2 = getInfoOut().getAnt2()[i] == 0 ? coeff1 : coeff2;
        float scale = sqrt(sc1 * sc2);
        for (std::size_t k = 0; k < n_correlations_; ++k) {
          data(i, j, k) =
              std::complex<float>(cnt + count_ * 10, cnt - 10 + count_ * 6) *
              scale;
          ++cnt;
        }
      }
    }
    xt::xtensor<float, 3> weights{shape};
    std::generate(weights.begin(), weights.end(), [result = 0.5f]() mutable {
      result += 0.01f;
      return result;
    });

    xt::xtensor<double, 2> uvw({n_baselines_, 3});
    for (std::size_t i = 0; i < n_baselines_; ++i) {
      uvw(i, 0) = 1 + count_ + i;
      uvw(i, 1) = 2 + count_ + i;
      uvw(i, 2) = 3 + count_ + i;
    }
    BOOST_CHECK(xt::allclose(buffer->GetData(), data));

    BOOST_CHECK_EQUAL_COLLECTIONS(buffer->GetFlags().shape().begin(),
                                  buffer->GetFlags().shape().end(),
                                  shape.begin(), shape.end());
    BOOST_CHECK(xt::all(xt::equal(buffer->GetFlags(), false)));

    BOOST_CHECK(xt::allclose(buffer->GetWeights(), weights));
    BOOST_CHECK(xt::allclose(buffer->GetUvw(), uvw));
    count_++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    Step::updateInfo(infoIn);
    BOOST_CHECK_EQUAL(infoIn.origNChan(), n_channels_);
    BOOST_CHECK_EQUAL(infoIn.nchan(), n_channels_);
    BOOST_CHECK_EQUAL(infoIn.ntime(), n_time_);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(infoIn.nbaselines(), n_baselines_);
  }

  int count_;
  int n_time_;
  std::size_t n_baselines_;
  std::size_t n_channels_;
  std::size_t n_correlations_;
};

void TestScaling(int ntime, int nbl, int nchan, int ncorr) {
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr);
  ParameterSet parset;
  parset.add("stations", "[rs01.s01, *]");
  parset.add("coeffs", "[[2,0.5],[3,2,1]]");
  parset.add("scalesize", "false");
  auto step2 =
      std::make_shared<ScaleData>(parset, "", ScaleData::MsType::kRegular);
  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_CASE(test_scale_data1) { TestScaling(2, 4, 4, 1); }

BOOST_AUTO_TEST_CASE(test_scale_data2) { TestScaling(10, 16, 32, 4); }

BOOST_AUTO_TEST_CASE(test_scale_data3) { TestScaling(10, 12, 16, 2); }

BOOST_AUTO_TEST_SUITE_END()
