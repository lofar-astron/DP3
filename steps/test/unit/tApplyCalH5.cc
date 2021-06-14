// tApplyCalH5.cc: Test program for class ApplyCal
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/soltab.h>

#include "tStepCommon.h"
#include "../../Step.h"
#include "../../ApplyCal.h"
#include "../../InputStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using casacore::max;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::ApplyCal;
using dp3::steps::InputStep;
using dp3::steps::Step;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;
using std::complex;
using std::vector;

BOOST_AUTO_TEST_SUITE(applycalh5)

// Simple class to generate input arrays.
// 9 baselines, 3 antennas, 4 correlations
class TestInput : public InputStep {
 public:
  TestInput(unsigned int ntime, unsigned int nchan)
      : itsCount(0),
        itsNTime(ntime),
        itsNChan(nchan),
        itsNBl(9),
        itsNCorr(4),
        itsTimeInterval(5.),
        itsFirstTime(4472025740.0) {
    info().init(itsNCorr, 0, nchan, ntime, itsFirstTime, itsTimeInterval,
                string(), string());
    // Fill the baseline stations; use 3 stations.
    // So they are called 00 01 02 10 11 12 20 21 22, etc.

    vector<int> ant1(itsNBl);
    vector<int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (int i = 0; i < itsNBl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 3) {
        st2 = 0;
        if (++st1 == 3) {
          st1 = 0;
        }
      }
    }
    vector<string> antNames{"ant1", "ant2", "ant3"};
    // Define their positions (more or less WSRT RT0-3).
    vector<casacore::MPosition> antPos(3);
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
    vector<double> antDiam(3, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    vector<double> chanWidth(nchan, 100.e6);
    vector<double> chanFreqs(nchan);
    for (unsigned int ch = 0; ch < nchan; ++ch) {
      double freq = 100.e6 + ch * 10.e6;
      if (ch > 2) {
        // Make frequencies unevenly spaced
        freq += 4.e6;
      }
      if (ch > 4) {
        freq += 1.e6;
      }
      if (ch > 5) {
        freq += 15.e6;
      }
      chanFreqs[ch] = freq;
    }
    if (nchan == 2) {
      chanFreqs[0] = 100.e6;
      chanFreqs[1] = 101.e6;
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
      data.data()[i] = casacore::Complex(1, 0);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    weights = 1.;

    casacore::Matrix<double> uvw(3, itsNBl);
    for (int i = 0; i < itsNBl; ++i) {
      uvw(0, i) = 1 + itsCount + i;
      uvw(1, i) = 2 + itsCount + i;
      uvw(2, i) = 3 + itsCount + i;
    }
    DPBuffer buf;
    buf.setTime(itsCount * itsTimeInterval + itsFirstTime);
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

  int itsCount, itsNTime, itsNChan, itsNBl, itsNCorr, itsTimeInterval;
  double itsFirstTime;
};

// Class to check result of TestInput run by tests.
class TestOutput : public Step {
 public:
  enum tests {
    WeightsNotChanged = 1,
    DataNotChanged = 2,
    DataChanged = 4,
    WeightEquals = 8
  };
  TestOutput(int ntime, int nchan, int doTest, bool solsHadFreqAxis = true,
             bool solsHadTimeAxis = true,
             JonesParameters::MissingAntennaBehavior missingAntennaBehavior =
                 JonesParameters::MissingAntennaBehavior::kError)
      : itsCount(0),
        itsTimeStep(0),
        itsNTime(ntime),
        itsNBl(9),
        itsNChan(nchan),
        itsNCorr(4),
        itsTimeInterval(5.),
        itsDoTest(doTest),
        itsSolsHadFreqAxis(solsHadFreqAxis),
        itsSolsHadTimeAxis(solsHadTimeAxis),
        itsMissingAntennaBehavior(missingAntennaBehavior) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    casacore::Cube<casacore::Complex> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = casacore::Complex(1, 0);
    }
    casacore::Cube<float> weights(itsNCorr, itsNChan, itsNBl);
    indgen(weights, 1.0f, 0.0f);

    vector<double> rightTimes(std::max(itsNTime, 5));
    rightTimes[0] = 0;
    rightTimes[1] = 2;
    rightTimes[2] = 3;
    for (int t = 3; t < itsNTime; ++t) {
      rightTimes[t] = 4;
    }
    if (!itsSolsHadTimeAxis) {
      rightTimes.assign(itsNTime, 0);
    }

    vector<double> rightFreqs(std::max(itsNChan, 5));
    rightFreqs[0] = 1;
    rightFreqs[1] = 1;
    rightFreqs[2] = 2;
    rightFreqs[3] = 2;
    rightFreqs[4] = 2;
    for (int f = 5; f < itsNChan; ++f) {
      rightFreqs[f] = 3;
    }
    if (!itsSolsHadFreqAxis) {
      rightFreqs.assign(itsNChan, 1);
    }

    if (itsDoTest) {
      for (unsigned int bl = 0; bl < info().nbaselines(); ++bl) {
        unsigned int ant1 = info().getAnt1()[bl];
        unsigned int ant2 = info().getAnt2()[bl];

        for (int chan = 0; chan < itsNChan; ++chan) {
          // Square root of autocorrelation for first antenna
          complex<float> val = sqrt(
              buf.getData().data()[bl * itsNCorr * itsNChan + chan * itsNCorr]);

          bool flag =
              buf.getFlags().data()[bl * itsNCorr * itsNChan + chan * itsNCorr];
          if ((ant1 == 1 || ant2 == 1) && rightTimes[itsTimeStep] == 2 &&
              rightFreqs[chan] == 2) {
            BOOST_CHECK(flag);
          } else if (itsMissingAntennaBehavior ==
                         JonesParameters::MissingAntennaBehavior::kFlag &&
                     (ant1 == 2 || ant2 == 2)) {
            BOOST_CHECK(flag);
          } else if (itsMissingAntennaBehavior ==
                     JonesParameters::MissingAntennaBehavior::kUnit) {
            BOOST_CHECK(!flag);
            if (ant1 == 2 && ant2 == 2) {
              BOOST_CHECK_CLOSE(1., val, 1.e-3);
            } else if (ant1 != 2 && ant2 != 2) {
              BOOST_CHECK_CLOSE(
                  rightTimes[itsTimeStep] * 100 + rightFreqs[chan], val, 1.e-3);
            } else {
              BOOST_CHECK_CLOSE(
                  rightTimes[itsTimeStep] * 100 + rightFreqs[chan], val * val,
                  1.e-3);
            }
          } else {
            BOOST_CHECK(!flag);
            BOOST_CHECK_CLOSE(rightTimes[itsTimeStep] * 100 + rightFreqs[chan],
                              val, 1.e-3);
          }
        }
      }
    }

    if (itsDoTest & DataNotChanged) {
      BOOST_CHECK(allNear(buf.getData(), data, 1.e-7));
    }
    if (itsDoTest & DataChanged) {
      BOOST_CHECK(!(allNear(buf.getData(), data, 1.e-7)));
    }
    if (itsDoTest & WeightsNotChanged) {
      BOOST_CHECK(allNear(buf.getWeights(), weights, 1.e-7));
    }
    itsCount++;
    itsTimeStep++;
    return true;
  }

  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo& infoIn) {
    info() = infoIn;
    BOOST_CHECK_EQUAL(itsNChan, int(infoIn.origNChan()));
    BOOST_CHECK_EQUAL(itsNChan, int(infoIn.nchan()));
    BOOST_CHECK_EQUAL(itsNTime, int(infoIn.ntime()));
    BOOST_CHECK_EQUAL(itsTimeInterval, infoIn.timeInterval());
    BOOST_CHECK_EQUAL(itsNBl, int(infoIn.nbaselines()));
  }

  int itsCount;
  int itsTimeStep;
  int itsNTime, itsNBl, itsNChan, itsNCorr, itsTimeInterval, itsDoTest;
  bool itsSolsHadFreqAxis, itsSolsHadTimeAxis;
  JonesParameters::MissingAntennaBehavior itsMissingAntennaBehavior;
};

// Test amplitude correction
void testampl(int ntime, int nchan, bool freqaxis, bool timeaxis) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nchan);
  Step::ShPtr step1(in);

  ParameterSet parset1;
  parset1.add("correction", "myampl");
  parset1.add("parmdb", "tApplyCalH5_tmp.h5");
  Step::ShPtr step2(new ApplyCal(in, parset1, ""));

  Step::ShPtr step3(new TestOutput(ntime, nchan, TestOutput::WeightsNotChanged,
                                   freqaxis, timeaxis));

  dp3::steps::test::Execute({step1, step2, step3});
}

// Test with missing antenna option
void testmissingant(int ntime, int nchan, string missingant,
                    bool solshadfreqaxis = false,
                    bool solshadtimeaxis = false) {
  TestInput* in = new TestInput(ntime, nchan);
  Step::ShPtr step1(in);

  ParameterSet parset;
  parset.add("myapplycal.correction", "myampl");
  parset.add("myapplycal.parmdb", "tApplyCalH5_tmp.h5");
  parset.add("myapplycal.missingantennabehavior", missingant);

  Step::ShPtr step2(new ApplyCal(in, parset, "myapplycal."));
  Step::ShPtr step3(new TestOutput(
      ntime, nchan, TestOutput::WeightsNotChanged, solshadfreqaxis,
      solshadtimeaxis,
      JonesParameters::StringToMissingAntennaBehavior(missingant)));

  std::vector<std::string> unusedKeys = parset.unusedKeys();
  BOOST_CHECK(unusedKeys.empty());

  dp3::steps::test::Execute({step1, step2, step3});
}

// Write a temporary H5Parm
void createH5Parm(vector<double> times, vector<double> freqs,
                  bool missingAnt = false) {
  H5Parm h5parm("tApplyCalH5_tmp.h5", true);

  // Add antenna metadata
  std::vector<std::string> antNames;
  size_t nAntennas = 3;
  if (missingAnt) {
    nAntennas = 2;
  }
  std::vector<std::array<double, 3>> antPositions;
  const std::array<double, 3> oneAntPos{42.0};
  for (unsigned int i = 0; i < nAntennas; ++i) {
    std::stringstream antNameStr;
    antNameStr << "ant" << (i + 1);
    antNames.push_back(antNameStr.str());
    antPositions.push_back(oneAntPos);
  }
  h5parm.AddAntennas(antNames, antPositions);

  vector<schaapcommon::h5parm::AxisInfo> axes;
  axes.push_back(schaapcommon::h5parm::AxisInfo("ant", nAntennas));
  if (!times.empty()) {
    axes.push_back(schaapcommon::h5parm::AxisInfo("time", times.size()));
  }
  if (!freqs.empty()) {
    axes.push_back(schaapcommon::h5parm::AxisInfo("freq", freqs.size()));
  }

  SolTab soltab = h5parm.CreateSolTab("myampl", "amplitude", axes);
  BOOST_CHECK_EQUAL(size_t{1}, h5parm.NumSolTabs());
  BOOST_CHECK(h5parm.HasSolTab("myampl"));
  soltab.SetTimes(times);
  soltab.SetFreqs(freqs);
  soltab.SetAntennas(antNames);

  unsigned int ntimes = max(times.size(), 1);
  unsigned int nfreqs = max(freqs.size(), 1);
  vector<double> values(ntimes * nfreqs * nAntennas);
  vector<double> weights(ntimes * nfreqs * nAntennas);
  for (unsigned int ant = 0; ant < nAntennas; ++ant) {
    for (unsigned int t = 0; t < ntimes; ++t) {
      for (unsigned int f = 0; f < nfreqs; ++f) {
        values[ant * ntimes * nfreqs + t * nfreqs + f] =
            1. / (100. * (t % 100) + (1 + f));
        weights[ant * ntimes * nfreqs + t * nfreqs + f] = 1.;
        if (ant == 1 && t == 2 && f == 1) {
          weights[ant * ntimes * nfreqs + t * nfreqs + f] = 0.;
        }
      }
    }
  }
  soltab.SetValues(values, weights, "CREATE with DPPP tApplyCalH5");
}

BOOST_AUTO_TEST_CASE(testampl1) {
  const vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                             4472025748.0, 4472025762.0};
  const vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs);
  testampl(5, 7, true, true);
}

BOOST_AUTO_TEST_CASE(testampl2) {
  const vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                             4472025748.0, 4472025762.0};
  const vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs);
  testampl(5, 2, true, true);
}

BOOST_AUTO_TEST_CASE(testampl3) {
  const vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                             4472025748.0, 4472025762.0};
  createH5Parm(times, vector<double>());
  testampl(8, 9, false, true);
}

BOOST_AUTO_TEST_CASE(testampl4) {
  const vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(vector<double>(), freqs);
  testampl(13, 3, true, false);
}

BOOST_AUTO_TEST_CASE(testampl5) {
  createH5Parm(vector<double>(), vector<double>());
  testampl(9, 2, false, false);
}

// Check an exception message starts with a given string
bool checkMissingAntError(const std::exception& ex) {
  BOOST_CHECK_EQUAL(ex.what(),
                    std::string("SolTab has no element ant3 in ant"));
  return true;
}
BOOST_AUTO_TEST_CASE(testmissingant_error) {
  createH5Parm(vector<double>(), vector<double>(), true);
  BOOST_CHECK_EXCEPTION(testmissingant(9, 2, "error"), std::exception,
                        checkMissingAntError);
}

bool checkWrongArgError(const std::exception& ex) {
  BOOST_CHECK_EQUAL(
      ex.what(), std::string("missingantennabehavior should be one of 'error', "
                             "'flag' or 'unit', not 'wrongargument'"));
  return true;
}
BOOST_AUTO_TEST_CASE(testmissingant_wrongarg) {
  createH5Parm(vector<double>(), vector<double>(), true);
  BOOST_CHECK_EXCEPTION(testmissingant(9, 2, "wrongargument"), std::exception,
                        checkWrongArgError);
}

BOOST_AUTO_TEST_CASE(testmissingant_flag) {
  createH5Parm(vector<double>(), vector<double>(), true);
  testmissingant(9, 2, "flag");
}

BOOST_AUTO_TEST_CASE(testmissingant_unit) {
  const vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                             4472025748.0, 4472025762.0};
  const vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs, true);
  testmissingant(9, 2, "unit", true, true);
}

BOOST_AUTO_TEST_SUITE_END()
