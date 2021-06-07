// tPhaseShift.cc: Test program for class PhaseShift
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../PhaseShift.h"

#include "tStepCommon.h"

#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"

#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::PhaseShift;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(phaseshift)

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public InputStep {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {
    info().init(ncorr, 0, nchan, ntime, 0., 10., string(), string());
    casacore::MDirection phaseCenter(casacore::Quantity(45, "deg"),
                                     casacore::Quantity(30, "deg"),
                                     casacore::MDirection::J2000);
    info().set(casacore::MPosition(), phaseCenter, phaseCenter, phaseCenter);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 100000.);
    std::vector<double> chanFreqs;
    for (int i = 0; i < nchan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chanFreqs), std::move(chanWidth));
    // Fill the baseline stations.
    // Determine nr of stations using:  na*(na+1)/2 = nbl
    // If many baselines, divide into groups of 6 to test if
    // PhaseShift disentangles it correctly.
    int nant = int(-0.5 + sqrt(0.25 + 2 * nbl));
    if (nant * (nant + 1) / 2 < nbl) ++nant;
    int grpszant = 3;
    int grpszbl = grpszant * (grpszant + 1) / 2;
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
    for (int i = 0; i < nbl; ++i) {
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
    vector<double> antDiam(nant, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    itsStatUVW.resize(3, nant);
    for (int i = 0; i < nant; ++i) {
      itsStatUVW(0, i) = 0.01 + i * 0.02;
      itsStatUVW(1, i) = 0.05 + i * 0.03;
      itsStatUVW(2, i) = 0.015 + i * 0.025;
    }
  }

  void fillUVW(casacore::Matrix<double>& uvw, int count) {
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = (itsStatUVW(0, getInfo().getAnt2()[i]) + count * 0.002 -
                   (itsStatUVW(0, getInfo().getAnt1()[i]) + count * 0.002));
      uvw(1, i) = (itsStatUVW(1, getInfo().getAnt2()[i]) + count * 0.004 -
                   (itsStatUVW(1, getInfo().getAnt1()[i]) + count * 0.004));
      uvw(2, i) = (itsStatUVW(2, getInfo().getAnt2()[i]) + count * 0.006 -
                   (itsStatUVW(2, getInfo().getAnt1()[i]) + count * 0.006));
    }
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
          casacore::Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
    }
    DPBuffer buf;
    buf.setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags(flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    casacore::Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags(fullResFlags);
    casacore::Matrix<double> uvw(3, itsNBl);
    fillUVW(uvw, itsCount);
    buf.setUVW(uvw);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void updateInfo(const DPInfo&) {}
  virtual void finish() { getNextStep()->finish(); }
  virtual void show(std::ostream&) const {}

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
  casacore::Matrix<double> itsStatUVW;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput : public Step {
 public:
  TestOutput(TestInput* input, int ntime, int nbl, int nchan, int ncorr,
             bool flag)
      : itsInput(input),
        itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    casacore::Cube<casacore::Complex> result(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(result.size()); ++i) {
      result.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
    }
    casacore::Matrix<double> uvw(3, itsNBl);
    itsInput->fillUVW(uvw, itsCount);
    // Check the result.
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-7));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-7));
    BOOST_CHECK(allEQ(buf.getFlags(), itsFlag));
    BOOST_CHECK_CLOSE_FRACTION(buf.getTime(), 2. + 5 * itsCount, 1e-6);
    BOOST_CHECK(allNear(buf.getUVW(), uvw, 1e-7));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    casacore::MVDirection dir = infoIn.phaseCenter().getValue();
    BOOST_CHECK_CLOSE_FRACTION(dir.getLong("deg").getValue(), 45., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(dir.getLat("deg").getValue(), 30., 1e-6);
  }

  TestInput* itsInput;
  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result of null phase-shifted TestInput.
class TestOutput1 : public Step {
 public:
  TestOutput1(TestInput* input, int ntime, int nbl, int nchan, int ncorr,
              bool flag)
      : itsInput(input),
        itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    casacore::Cube<casacore::Complex> result(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(result.size()); ++i) {
      result.data()[i] =
          casacore::Complex(i + itsCount * 10, i - 1000 + itsCount * 6);
    }
    casacore::Matrix<double> uvw(3, itsNBl);
    itsInput->fillUVW(uvw, itsCount);
    // Check the result.
    BOOST_CHECK(!allNear(real(buf.getData()), real(result), 1e-5));
    BOOST_CHECK(!allEQ(real(buf.getData()), real(result)));
    BOOST_CHECK(!allNear(imag(buf.getData()), imag(result), 1e-5));
    BOOST_CHECK(!allEQ(imag(buf.getData()), imag(result)));
    BOOST_CHECK(allEQ(buf.getFlags(), itsFlag));
    BOOST_CHECK_CLOSE_FRACTION(buf.getTime(), 2. + 5 * itsCount, 1e-5);
    BOOST_CHECK(!allNear(buf.getUVW(), uvw, 1e-5));
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    casacore::MVDirection dir = infoIn.phaseCenter().getValue();
    BOOST_CHECK_CLOSE_FRACTION(dir.getLong("deg").getValue(), 50., 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(dir.getLat("deg").getValue(), 35., 1e-5);
  }

  TestInput* itsInput;
  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Test with a shift to the original center.
void test1(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  // Keep phase center the same to be able to check if data are correct.
  parset.add("phasecenter", "[45deg, 30deg]");
  Step::ShPtr step2(new PhaseShift(in, parset, ""));
  Step::ShPtr step3(new TestOutput(in, ntime, nbl, nchan, ncorr, flag));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test with a shift to another and then to the original phase center.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nbl, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  // First shift to another center, then back to original.
  ParameterSet parset;
  parset.add("phasecenter", "[50deg, 35deg]");
  ParameterSet parset1;
  parset1.add("phasecenter", "[]");
  Step::ShPtr step2(new PhaseShift(in, parset, ""));
  Step::ShPtr step3(new TestOutput1(in, ntime, nbl, nchan, ncorr, flag));
  Step::ShPtr step4(new PhaseShift(in, parset1, ""));
  Step::ShPtr step5(new TestOutput(in, ntime, nbl, nchan, ncorr, flag));
  dp3::steps::test::Execute({step1, step2, step3, step4, step5});
}

BOOST_AUTO_TEST_CASE(test1a) { test1(10, 3, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test1b) { test1(10, 10, 30, 1, true); }

BOOST_AUTO_TEST_CASE(test2a) { test2(10, 6, 32, 4, false); }

BOOST_AUTO_TEST_CASE(test2b) { test2(10, 6, 30, 1, true); }

BOOST_AUTO_TEST_SUITE_END()
