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
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include "../../ApplyCal.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::ApplyCal;
using dp3::steps::Step;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

BOOST_AUTO_TEST_SUITE(applycalh5)

// Simple class to generate input arrays.
// 9 baselines, 3 antennas, 4 correlations
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(unsigned int ntime, unsigned int nchan)
      : itsCount(0),
        itsNTime(ntime),
        itsNChan(nchan),
        itsNBl(9),
        itsNCorr(4),
        itsTimeInterval(5.),
        itsFirstTime(4472025742.5) {
    GetWritableInfoOut() = DPInfo(itsNCorr, itsNChan);
    GetWritableInfoOut().setTimes(
        itsFirstTime, itsFirstTime + (itsNTime - 1) * itsTimeInterval,
        itsTimeInterval);
    // Fill the baseline stations; use 3 stations.
    // So they are called 00 01 02 10 11 12 20 21 22, etc.

    std::vector<int> ant1(itsNBl);
    std::vector<int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (std::size_t i = 0; i < itsNBl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 3) {
        st2 = 0;
        if (++st1 == 3) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antNames{"ant1", "ant2", "ant3"};
    // Define their positions (more or less WSRT RT0-3).
    std::vector<casacore::MPosition> antPos(3);
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
    std::vector<double> antDiam(3, 70.);
    GetWritableInfoOut().setAntennas(antNames, antDiam, antPos, ant1, ant2);
    // Define the frequencies.
    std::vector<double> chanWidth(nchan, 100.e6);
    std::vector<double> chanFreqs(nchan);
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
    GetWritableInfoOut().setChannels(std::move(chanFreqs),
                                     std::move(chanWidth));
  }

 private:
  bool process(std::unique_ptr<DPBuffer>) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    auto buffer = std::make_unique<DPBuffer>();
    buffer->SetTime(itsCount * itsTimeInterval + itsFirstTime);

    const std::array<std::size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    buffer->GetData().resize(shape);
    buffer->GetData().fill(std::complex<float>{1, 0});

    buffer->GetWeights().resize(shape);
    buffer->GetWeights().fill(1.0f);

    buffer->GetFlags().resize(shape);
    buffer->GetFlags().fill(false);

    buffer->GetUvw().resize({itsNBl, 3});
    for (std::size_t i = 0; i < itsNBl; ++i) {
      buffer->GetUvw()(i, 0) = 1 + itsCount + i;
      buffer->GetUvw()(i, 1) = 2 + itsCount + i;
      buffer->GetUvw()(i, 2) = 3 + itsCount + i;
    }

    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  std::size_t itsCount, itsNTime, itsNChan, itsNBl, itsNCorr, itsTimeInterval;
  double itsFirstTime;
};

// Class to check result of TestInput run by tests.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  enum tests {
    WeightsNotChanged = 1,
    DataNotChanged = 2,
    DataChanged = 4,
    WeightEquals = 8
  };
  TestOutput(std::size_t ntime, std::size_t nchan, bool doTest,
             bool solsHadFreqAxis = true, bool solsHadTimeAxis = true,
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
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    const std::array<std::size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> data(shape, std::complex<float>{1, 0});
    xt::xtensor<float, 3> weights(shape, 1.0f);

    std::vector<double> rightTimes(std::max<std::size_t>(itsNTime, 5));
    rightTimes[0] = 0;
    rightTimes[1] = 2;
    rightTimes[2] = 3;
    for (std::size_t t = 3; t < itsNTime; ++t) {
      rightTimes[t] = 4;
    }
    if (!itsSolsHadTimeAxis) {
      rightTimes.assign(itsNTime, 0);
    }

    std::vector<double> rightFreqs(std::max<std::size_t>(itsNChan, 5));
    rightFreqs[0] = 1;
    rightFreqs[1] = 1;
    rightFreqs[2] = 2;
    rightFreqs[3] = 2;
    rightFreqs[4] = 2;
    for (std::size_t f = 5; f < itsNChan; ++f) {
      rightFreqs[f] = 3;
    }
    if (!itsSolsHadFreqAxis) {
      rightFreqs.assign(itsNChan, 1);
    }

    if (itsDoTest) {
      for (unsigned int bl = 0; bl < getInfoOut().nbaselines(); ++bl) {
        unsigned int ant1 = getInfoOut().getAnt1()[bl];
        unsigned int ant2 = getInfoOut().getAnt2()[bl];

        for (std::size_t chan = 0; chan < itsNChan; ++chan) {
          // Square root of autocorrelation for first antenna
          std::complex<float> val = sqrt(buffer->GetData()(bl, chan, 0));
          bool flag = buffer->GetFlags()(bl, chan, 0);
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
      BOOST_CHECK(xt::allclose(buffer->GetData(), data, 1.0e-7));
    }
    if (itsDoTest & DataChanged) {
      BOOST_CHECK(!(xt::allclose(buffer->GetData(), data, 1.0e-7)));
    }
    if (itsDoTest & WeightsNotChanged) {
      BOOST_CHECK(xt::allclose(buffer->GetWeights(), weights, 1.0e-7));
    }
    itsCount++;
    itsTimeStep++;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    Step::updateInfo(infoIn);
    BOOST_CHECK_EQUAL(itsNChan, infoIn.origNChan());
    BOOST_CHECK_EQUAL(itsNChan, infoIn.nchan());
    BOOST_CHECK_EQUAL(itsNTime, infoIn.ntime());
    BOOST_CHECK_EQUAL(itsTimeInterval, infoIn.timeInterval());
    BOOST_CHECK_EQUAL(itsNBl, infoIn.nbaselines());
  }

  int itsCount;
  int itsTimeStep;
  std::size_t itsNTime;
  std::size_t itsNBl;
  std::size_t itsNChan;
  std::size_t itsNCorr;
  std::size_t itsTimeInterval;
  bool itsDoTest;
  bool itsSolsHadFreqAxis;
  bool itsSolsHadTimeAxis;
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
  auto step2 = std::make_shared<ApplyCal>(parset1, "");

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

  auto step2 = std::make_shared<ApplyCal>(parset, "myapplycal.");
  Step::ShPtr step3(new TestOutput(
      ntime, nchan, TestOutput::WeightsNotChanged, solshadfreqaxis,
      solshadtimeaxis,
      JonesParameters::StringToMissingAntennaBehavior(missingant)));

  std::vector<std::string> unusedKeys = parset.unusedKeys();
  BOOST_CHECK(unusedKeys.empty());

  dp3::steps::test::Execute({step1, step2, step3});
}

// Write a temporary H5Parm
void createH5Parm(std::vector<double> times, std::vector<double> freqs,
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

  std::vector<schaapcommon::h5parm::AxisInfo> axes;
  axes.push_back({"ant", static_cast<unsigned>(nAntennas)});
  if (!times.empty()) {
    axes.push_back({"time", static_cast<unsigned>(times.size())});
  }
  if (!freqs.empty()) {
    axes.push_back({"freq", static_cast<unsigned>(freqs.size())});
  }

  SolTab soltab = h5parm.CreateSolTab("myampl", "amplitude", axes);
  BOOST_CHECK_EQUAL(size_t{1}, h5parm.NumSolTabs());
  BOOST_CHECK(h5parm.HasSolTab("myampl"));
  soltab.SetTimes(times);
  soltab.SetFreqs(freqs);
  soltab.SetAntennas(antNames);

  const size_t ntimes = std::max<size_t>(times.size(), 1u);
  const size_t nfreqs = std::max<size_t>(freqs.size(), 1u);
  std::vector<double> values(ntimes * nfreqs * nAntennas);
  std::vector<double> weights(ntimes * nfreqs * nAntennas);
  for (size_t ant = 0; ant < nAntennas; ++ant) {
    for (size_t t = 0; t < ntimes; ++t) {
      for (size_t f = 0; f < nfreqs; ++f) {
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
  const std::vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                                  4472025748.0, 4472025762.0};
  const std::vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs);
  testampl(5, 7, true, true);
}

BOOST_AUTO_TEST_CASE(testampl2) {
  const std::vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                                  4472025748.0, 4472025762.0};
  const std::vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs);
  testampl(5, 2, true, true);
}

BOOST_AUTO_TEST_CASE(testampl3) {
  const std::vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                                  4472025748.0, 4472025762.0};
  createH5Parm(times, std::vector<double>());
  testampl(8, 9, false, true);
}

BOOST_AUTO_TEST_CASE(testampl4) {
  const std::vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(std::vector<double>(), freqs);
  testampl(13, 3, true, false);
}

BOOST_AUTO_TEST_CASE(testampl5) {
  createH5Parm(std::vector<double>(), std::vector<double>());
  testampl(9, 2, false, false);
}

// Check an exception message starts with a given string
bool checkMissingAntError(const std::exception& ex) {
  BOOST_CHECK_EQUAL(ex.what(),
                    std::string("SolTab has no element ant3 in ant"));
  return true;
}
BOOST_AUTO_TEST_CASE(testmissingant_error) {
  createH5Parm(std::vector<double>(), std::vector<double>(), true);
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
  createH5Parm(std::vector<double>(), std::vector<double>(), true);
  BOOST_CHECK_EXCEPTION(testmissingant(9, 2, "wrongargument"), std::exception,
                        checkWrongArgError);
}

BOOST_AUTO_TEST_CASE(testmissingant_flag) {
  createH5Parm(std::vector<double>(), std::vector<double>(), true);
  testmissingant(9, 2, "flag");
}

BOOST_AUTO_TEST_CASE(testmissingant_unit) {
  const std::vector<double> times{4472025742.0, 4472025745.0, 4472025747.5,
                                  4472025748.0, 4472025762.0};
  const std::vector<double> freqs{90.e6, 139.e6, 170.e6};
  createH5Parm(times, freqs, true);
  testmissingant(9, 2, "unit", true, true);
}

BOOST_AUTO_TEST_SUITE_END()
