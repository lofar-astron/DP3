// tStationAdder.cc: Test program for class StationAdder
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include "../../StationAdder.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::StationAdder;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(stationadder)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(std::size_t ntime, std::size_t nbl, std::size_t nchan,
            std::size_t ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {
    info() = DPInfo(ncorr, nchan);
    info().setTimes(0.0, (ntime - 1) * 5.0, 5.0);
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(nbl);
    vector<int> ant2(nbl);
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
    vector<string> antNames{"rs01.s01", "rs02.s01", "cs01.s01", "cs01.s02"};
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
    vector<double> antDiam(4, 70.0);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 1000000.0);
    std::vector<double> chanFreqs;
    for (std::size_t i = 0; i < nchan; i++) {
      chanFreqs.push_back(10500000.0 + i * 1000000.0);
    }
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
  }

 private:
  bool process(std::unique_ptr<DPBuffer>) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }

    const std::array<std::size_t, 3> data_shape{itsNBl, itsNChan, itsNCorr};

    xt::xtensor<std::complex<float>, 3> data(data_shape);
    for (int i = 0; i < int(data.size()); ++i) {
      data[i] = std::complex<float>(i + itsCount * 10.0f,
                                    i - 10.0f + itsCount * 6.0f);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 0.5f, 0.01f);
    casacore::Matrix<double> uvw(3, itsNBl);
    for (std::size_t i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime(itsCount * 30 + 4472025740.0);
    buf.ResizeData(data_shape);
    buf.GetData() = data;
    buf.setWeights(weights);
    buf.setUVW(uvw);
    buf.ResizeFlags(data_shape);
    buf.GetFlags().fill(false);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  std::size_t itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of TestInput run by test1.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, std::size_t nbl, std::size_t nchan, std::size_t ncorr,
             bool sumauto)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsSumAuto(sumauto) {}

 private:
  void addData(xt::xtensor<std::complex<float>, 2>& to,
               const xt::xtensor<std::complex<float>, 3>& from, int bl) {
    to += xt::view(from, bl, xt::all(), xt::all());
  }
  void addConjData(xt::xtensor<std::complex<float>, 2>& to,
                   const xt::xtensor<std::complex<float>, 3>& from, int bl) {
    to += xt::conj(xt::view(from, bl, xt::all(), xt::all()));
  }
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    const std::array<std::size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> data(shape);

    for (int i = 0; i < int(data.size()); ++i) {
      data[i] = std::complex<float>(i + itsCount * 10, i - 10 + itsCount * 6);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 0.5f, 0.01f);
    xt::xtensor<std::complex<float>, 2> databl0(
        std::array<std::size_t, 2>{itsNChan, itsNCorr});
    xt::xtensor<std::complex<float>, 2> databl1(
        std::array<std::size_t, 2>{itsNChan, itsNCorr});
    // "{ns:[rs01.s01, rs02.s01, cs01.s02]}" was given resulting in 2 new
    // baselines (ns-ns and cs01.s01-ns).
    // Thus adding the baselines below.
    float weight = 0;
    if (itsSumAuto) {
      // add autocorr to form new autocorr
      addData(databl0, data, 0);
      addData(databl0, data, 5);
      addData(databl0, data, 15);
      weight = 3;
    } else {
      // add crosscorr to form new autocorr
      addData(databl0, data, 1);
      addData(databl0, data, 3);
      addData(databl0, data, 4);
      addData(databl0, data, 7);
      addData(databl0, data, 12);
      addData(databl0, data, 13);
      weight = 6;
    }
    addData(databl1, data, 8);
    addData(databl1, data, 9);
    addData(databl1, data, 11);
    addConjData(databl1, data, 2);
    addConjData(databl1, data, 6);
    addConjData(databl1, data, 14);
    xt::xtensor<double, 2> uvw(std::array<std::size_t, 2>{itsNBl, 3});
    for (std::size_t i = 0; i < itsNBl; ++i) {
      uvw(i, 0) = 1 + itsCount + i;
      uvw(i, 1) = 2 + itsCount + i;
      uvw(i, 2) = 3 + itsCount + i;
    }

    auto xt_weights = aocommon::xt::CreateSpan(
        weights.data(),
        std::array{static_cast<std::size_t>(weights.shape()[2]),
                   static_cast<std::size_t>(weights.shape()[1]),
                   static_cast<std::size_t>(weights.shape()[0])});

    BOOST_TEST(data == xt::view(buffer->GetData(), xt::range(0, itsNBl),
                                xt::all(), xt::all()));
    BOOST_TEST(buffer->GetFlags().shape() ==
                   (std::array<std::size_t, 3>{itsNBl + 2, itsNChan, itsNCorr}),
               boost::test_tools::per_element());
    BOOST_CHECK(!xt::any(buffer->GetFlags()));
    BOOST_TEST(xt_weights == xt::view(buffer->GetWeights(),
                                      xt::range(0, itsNBl), xt::all(),
                                      xt::all()));
    BOOST_TEST(uvw ==
               xt::view(buffer->GetUvw(), xt::range(0, itsNBl), xt::all()));

    // Now check data of new baselines.
    BOOST_CHECK(
        xt::allclose(xt::view(buffer->GetData(), itsNBl, xt::all(), xt::all()),
                     (databl0 / weight)));
    BOOST_CHECK(xt::allclose(
        xt::view(buffer->GetWeights(), itsNBl, xt::all(), xt::all()), weight));
    BOOST_CHECK(xt::allclose(
        xt::view(buffer->GetData(), itsNBl + 1, xt::all(), xt::all()),
        (databl1 / 6.0f)));
    BOOST_CHECK(xt::allclose(
        xt::view(buffer->GetWeights(), itsNBl + 1, xt::all(), xt::all()),
        6.0f));
    itsCount++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.nbaselines()), itsNBl + 2);
    BOOST_CHECK_EQUAL(int(infoIn.antennaNames().size()), 5);
    BOOST_CHECK_EQUAL(int(infoIn.antennaDiam().size()), 5);
    BOOST_CHECK_EQUAL(int(infoIn.antennaPos().size()), 5);
    BOOST_CHECK_EQUAL(infoIn.antennaNames()[4], "ns");
    casacore::Vector<double> pos1(infoIn.antennaPos()[4].getValue().getValue());
    BOOST_CHECK(
        casacore::near(pos1[0], (3828763.0 + 3828746.0 + 3828713.0) / 3));
    BOOST_CHECK(casacore::near(pos1[1], (442449.0 + 442592.0 + 442878.0) / 3));
    BOOST_CHECK(
        casacore::near(pos1[2], (5064923.0 + 5064924.0 + 5064926.0) / 3));
    // Check diam.
    double d1 = sqrt((pos1[0] - 3828763) * (pos1[0] - 3828763) +
                     (pos1[1] - 442449) * (pos1[1] - 442449) +
                     (pos1[2] - 5064923) * (pos1[2] - 5064923));
    double d2 = sqrt((pos1[0] - 3828746) * (pos1[0] - 3828746) +
                     (pos1[1] - 442592) * (pos1[1] - 442592) +
                     (pos1[2] - 5064924) * (pos1[2] - 5064924));
    double d3 = sqrt((pos1[0] - 3828713) * (pos1[0] - 3828713) +
                     (pos1[1] - 442878) * (pos1[1] - 442878) +
                     (pos1[2] - 5064926) * (pos1[2] - 5064926));
    BOOST_CHECK(casacore::near(infoIn.antennaDiam()[4],
                               70 + 2 * std::max(d1, std::max(d2, d3))));
  }

  int itsCount;
  std::size_t itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsSumAuto;
};

// Class to check result of flagged, unaveraged TestInput run by test2.
class TestOutput2 : public dp3::steps::test::ThrowStep {
 public:
  TestOutput2(int ntime, int nbl, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  void addData(casacore::Cube<casacore::Complex>& to,
               const casacore::Cube<casacore::Complex>& from,
               casacore::Cube<float>& tow, const casacore::Cube<float>& weights,
               int bl) {
    casacore::Cube<casacore::Complex> tmp = from.copy();
    casacore::Cube<casacore::Complex>::iterator tmpit = tmp.begin();
    casacore::Cube<float>::const_iterator weightit = weights.begin();
    for (; tmpit != tmp.end() && weightit != weights.end();
         tmpit++, weightit++) {
      *tmpit *= *weightit;
    }
    to += tmp(casacore::IPosition(3, 0, 0, bl),
              casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl));
    tow += weights(casacore::IPosition(3, 0, 0, bl),
                   casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl));
  }
  void addConjData(casacore::Cube<casacore::Complex>& to,
                   const casacore::Cube<casacore::Complex>& from,
                   casacore::Cube<float>& tow,
                   const casacore::Cube<float>& weights, int bl) {
    casacore::Cube<casacore::Complex> tmp = from.copy();
    casacore::Cube<casacore::Complex>::iterator tmpit = tmp.begin();
    casacore::Cube<float>::const_iterator weightit = weights.begin();
    for (; tmpit != tmp.end() && weightit != weights.end();
         tmpit++, weightit++) {
      *tmpit *= *weightit;
    }
    to +=
        conj(tmp(casacore::IPosition(3, 0, 0, bl),
                 casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl)));
    tow += weights(casacore::IPosition(3, 0, 0, bl),
                   casacore::IPosition(3, to.nrow() - 1, to.ncolumn() - 1, bl));
  }
  bool process(const DPBuffer& buf) override {
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 10 + itsCount * 6);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 0.5f, 0.01f);
    casacore::Cube<casacore::Complex> databl0(itsNCorr, itsNChan, 1);
    casacore::Cube<casacore::Complex> databl1(itsNCorr, itsNChan, 1);
    casacore::Cube<casacore::Complex> databl2(itsNCorr, itsNChan, 1);
    casacore::Cube<casacore::Complex> databl3(itsNCorr, itsNChan, 1);
    casacore::Cube<casacore::Complex> databl4(itsNCorr, itsNChan, 1);
    casacore::Cube<float> weightbl0(itsNCorr, itsNChan, 1, 0.);
    casacore::Cube<float> weightbl1(itsNCorr, itsNChan, 1, 0.);
    casacore::Cube<float> weightbl2(itsNCorr, itsNChan, 1, 0.);
    casacore::Cube<float> weightbl3(itsNCorr, itsNChan, 1, 0.);
    casacore::Cube<float> weightbl4(itsNCorr, itsNChan, 1, 0.);
    // "{ns1:[rs01.s01, rs02.s01], ns2:[cs01.s02, cs01.s01]}" was given.
    addData(databl0, data, weightbl0, weights, 8);
    addData(databl0, data, weightbl0, weights, 9);
    addData(databl1, data, weightbl1, weights, 12);
    addData(databl1, data, weightbl1, weights, 13);
    addData(databl2, data, weightbl2, weights, 2);
    addData(databl2, data, weightbl2, weights, 3);
    addData(databl3, data, weightbl3, weights, 6);
    addData(databl3, data, weightbl3, weights, 7);
    addConjData(databl0, data, weightbl0, weights, 2);
    addConjData(databl0, data, weightbl0, weights, 6);
    addConjData(databl1, data, weightbl1, weights, 3);
    addConjData(databl1, data, weightbl1, weights, 7);
    addConjData(databl2, data, weightbl2, weights, 8);
    addConjData(databl2, data, weightbl2, weights, 12);
    addConjData(databl3, data, weightbl3, weights, 9);
    addConjData(databl3, data, weightbl3, weights, 13);
    addConjData(databl4, databl0, weightbl4, weightbl0, 0);
    addConjData(databl4, databl1, weightbl4, weightbl1, 0);
    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    casacore::IPosition end(3, itsNCorr - 1, itsNChan - 1, itsNBl - 1);
    BOOST_CHECK(
        allEQ(buf.GetCasacoreData()(casacore::IPosition(3, 0), end), data));
    BOOST_CHECK_EQUAL(buf.GetCasacoreFlags().shape(),
                      casacore::IPosition(3, itsNCorr, itsNChan, itsNBl + 5));
    BOOST_CHECK(allEQ(buf.GetCasacoreFlags(), false));
    BOOST_CHECK(allEQ(buf.GetCasacoreWeights()(casacore::IPosition(3, 0), end),
                      weights));
    BOOST_CHECK(
        allEQ(buf.GetCasacoreUvw()(casacore::IPosition(2, 0),
                                   casacore::IPosition(2, 2, itsNBl - 1)),
              uvw));
    // Now check data of new baselines.
    end[2] = itsNBl;
    BOOST_CHECK(allNear(
        buf.GetCasacoreData()(casacore::IPosition(3, 0, 0, itsNBl), end),
        databl0, 1e-5));
    BOOST_CHECK(allNear(
        buf.GetCasacoreWeights()(casacore::IPosition(3, 0, 0, itsNBl), end),
        weightbl0, 1e-5));
    end[2] = itsNBl + 1;
    BOOST_CHECK(allNear(
        buf.GetCasacoreData()(casacore::IPosition(3, 0, 0, itsNBl + 1), end),
        databl1, 1e-5));
    BOOST_CHECK(allNear(
        buf.GetCasacoreWeights()(casacore::IPosition(3, 0, 0, itsNBl + 1), end),
        weightbl1, 1e-5));
    end[2] = itsNBl + 2;
    BOOST_CHECK(allNear(
        buf.GetCasacoreData()(casacore::IPosition(3, 0, 0, itsNBl + 2), end),
        databl2, 1e-5));
    BOOST_CHECK(allNear(
        buf.GetCasacoreWeights()(casacore::IPosition(3, 0, 0, itsNBl + 2), end),
        weightbl2, 1e-5));
    end[2] = itsNBl + 3;
    BOOST_CHECK(allNear(
        buf.GetCasacoreData()(casacore::IPosition(3, 0, 0, itsNBl + 3), end),
        databl3, 1e-5));
    BOOST_CHECK(allNear(
        buf.GetCasacoreWeights()(casacore::IPosition(3, 0, 0, itsNBl + 3), end),
        weightbl3, 1e-5));
    end[2] = itsNBl + 4;
    BOOST_CHECK(allNear(
        buf.GetCasacoreData()(casacore::IPosition(3, 0, 0, itsNBl + 4), end),
        databl4, 1e-5));
    BOOST_CHECK(allNear(
        buf.GetCasacoreWeights()(casacore::IPosition(3, 0, 0, itsNBl + 4), end),
        weightbl4, 1e-5));
    itsCount++;
    return true;
    BOOST_CHECK(allEQ(buf.GetCasacoreFlags(), false));
    itsCount++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.nbaselines()), itsNBl + 5);
    BOOST_CHECK_EQUAL(int(infoIn.antennaNames().size()), 6);
    BOOST_CHECK_EQUAL(infoIn.antennaNames()[4], "ns1");
    BOOST_CHECK_EQUAL(infoIn.antennaNames()[5], "ns2");
    casacore::Vector<double> pos1(infoIn.antennaPos()[4].getValue().getValue());
    BOOST_CHECK(casacore::near(pos1[0], (3828763. + 3828746.) / 2));
    BOOST_CHECK(casacore::near(pos1[1], (442449. + 442592.) / 2));
    BOOST_CHECK(casacore::near(pos1[2], (5064923. + 5064924.) / 2));
    casacore::Vector<double> pos2(infoIn.antennaPos()[5].getValue().getValue());
    BOOST_CHECK(casacore::near(pos2[0], (3828729. + 3828713.) / 2));
    BOOST_CHECK(casacore::near(pos2[1], (442735. + 442878.) / 2));
    BOOST_CHECK(casacore::near(pos2[2], (5064925. + 5064926.) / 2));
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Class to check result of TestInput run by test4.
class TestOutput4 : public dp3::steps::test::ThrowStep {
 public:
  TestOutput4(int ntime, int nbl, int nchan, int /*ncorr*/)
      : itsNTime(ntime), itsNBl(nbl), itsNChan(nchan) {}

 private:
  bool process(const DPBuffer&) override { return true; }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(int(infoIn.origNChan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.nchan()), itsNChan);
    BOOST_CHECK_EQUAL(int(infoIn.ntime()), itsNTime);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5);
    BOOST_CHECK_EQUAL(int(infoIn.nchanAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.ntimeAvg()), 1);
    BOOST_CHECK_EQUAL(int(infoIn.nbaselines()), itsNBl);
    BOOST_CHECK_EQUAL(int(infoIn.antennaNames().size()), 4);
    BOOST_CHECK_EQUAL(int(infoIn.antennaDiam().size()), 4);
    BOOST_CHECK_EQUAL(int(infoIn.antennaPos().size()), 4);
  }

  int itsNTime, itsNBl, itsNChan;
};

BOOST_DATA_TEST_CASE(test_add_three_stations,
                     boost::unit_test::data::make({true, false}), sumauto) {
  // Test must be done with with 16 baselines.
  const int kNTime = 10;
  const int kNBl = 16;
  const int kNChan = 32;
  const int kNCorr = 4;

  auto step1 = std::make_shared<TestInput>(kNTime, kNBl, kNChan, kNCorr);
  dp3::common::ParameterSet parset;
  parset.add("stations", "{ns:[rs01.s01, rs02.s01, cs01.s02]}");
  parset.add("autocorr", "true");
  if (!sumauto) {
    parset.add("sumauto", "false");
  }
  parset.add("average", "true");
  parset.add("useweights", "false");
  auto step2 = std::make_shared<StationAdder>(parset, "");
  auto step3 =
      std::make_shared<TestOutput>(kNTime, kNBl, kNChan, kNCorr, sumauto);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_CASE(test_add_two_groups_of_two_stations) {
  const int kNTime = 10;
  const int kNBl = 16;
  const int kNChan = 32;
  const int kNCorr = 4;

  auto step1 = std::make_shared<TestInput>(kNTime, kNBl, kNChan, kNCorr);
  dp3::common::ParameterSet parset;
  parset.add("stations",
             "{ns1:[rs01.s01, rs02.s01], ns2:[cs01.s02, cs01.s01]}");
  parset.add("autocorr", "false");
  parset.add("average", "false");
  auto step2 = std::make_shared<StationAdder>(parset, "");
  auto step3 = std::make_shared<TestOutput2>(kNTime, kNBl, kNChan, kNCorr);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_DATA_TEST_CASE(
    test_invalid_station,
    boost::unit_test::data::make(
        {// Unknown station.
         // Enabling this test crashes the test framework...
         // "{ns1:unknown, ns2:[cs01.s02, cs01.s01]}",

         // New station already used.
         "{ns1:[rs01.s01, rs02.s01], cs01.s02:[cs01.s02, cs01.s01]}",

         // Old station doubly used.
         "{ns1:[rs01.s01, rs02.s01], ns2:[rs01.s01, cs01.s01]}"}),
    stations) {
  auto step1 = std::make_shared<TestInput>(2, 8, 4, 4);
  dp3::common::ParameterSet parset;
  parset.add("stations", stations);
  parset.add("autocorr", "true");
  parset.add("average", "false");
  auto step2 = std::make_shared<StationAdder>(parset, "");
  auto step3 = std::make_shared<dp3::steps::test::ThrowStep>();
  BOOST_CHECK_THROW(dp3::steps::test::Execute({step1, step2, step3}),
                    std::exception);
}

// Test making a superstation out of nonexisting stations (should do nothing).
BOOST_AUTO_TEST_CASE(test_superstation_of_nonexisting_stations) {
  const int kNTime = 10;
  const int kNBl = 16;
  const int kNChan = 32;
  const int kNCorr = 4;

  auto step1 = std::make_shared<TestInput>(kNTime, kNBl, kNChan, kNCorr);
  dp3::common::ParameterSet parset;
  parset.add("stations", "{ns1:nonexistingstationpattern}");
  auto step2 = std::make_shared<StationAdder>(parset, "");
  auto step3 = std::make_shared<TestOutput4>(kNTime, kNBl, kNChan, kNCorr);
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_SUITE_END()
