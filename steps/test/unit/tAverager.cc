// tAverager.cc: Test program for class Averager
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <array>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include "steps/Averager.h"
#include "base/DPBuffer.h"
#include "base/DPInfo.h"
#include "common/ParameterSet.h"
#include "common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Averager;
using dp3::steps::Step;
using std::vector;

BOOST_AUTO_TEST_SUITE(averager)

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(int ntime, int nbl, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    buffer->SetTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo

    const std::array<std::size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    buffer->GetData().resize(shape);
    for (int i = 0; i < int(buffer->GetData().size()); ++i) {
      buffer->GetData().data()[i] =
          std::complex<float>(i + itsCount * 10, i - 1000 + itsCount * 6);
    }

    buffer->GetWeights().resize(shape);
    buffer->GetWeights().fill(1.0f);

    buffer->GetFlags().resize(shape);
    buffer->GetFlags().fill(itsFlag);

    buffer->GetUvw().resize({itsNBl, 3});
    xt::flatten(buffer->GetUvw()) =
        (itsCount * 100) + xt::arange(buffer->GetUvw().size());

    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Use timeInterval=5
    GetWritableInfoOut() = DPInfo(itsNCorr, itsNChan);
    GetWritableInfoOut().setTimes(100.0, 100.0 + (itsNTime - 1) * 5.0, 5.0);
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNChan, 100000.);
    for (std::size_t i = 0; i < itsNChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    GetWritableInfoOut().setChannels(std::move(chanFreqs),
                                     std::move(chanWidth));
  }

 private:
  int itsCount;
  int itsNTime;
  std::size_t itsNBl;
  std::size_t itsNChan;
  std::size_t itsNCorr;
  bool itsFlag;
};

// Class to check result of averaging TestInput.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr, int navgtime,
             int navgchan, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nbl),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsNAvgTime(navgtime),
        itsNAvgChan(navgchan),
        itsFlag(flag) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    const size_t nchan = 1 + (itsNChan - 1) / itsNAvgChan;
    const size_t navgtime =
        std::min(itsNAvgTime, itsNTime - itsCount * itsNAvgTime);
    // Fill expected result in similar way as TestInput.
    std::array<size_t, 3> shape{itsNBl, itsNChan, itsNCorr};
    xt::xtensor<std::complex<float>, 3> data(shape, 0.0f);
    xt::xtensor<float, 3> weights(shape, 0.0f);
    for (size_t j = itsCount * itsNAvgTime;
         j < itsCount * itsNAvgTime + navgtime; ++j) {
      for (size_t i = 0; i < data.size(); ++i) {
        data.data()[i] +=
            std::complex<float>(i + j * 10, int(i) - 1000 + int(j) * 6);
        weights.data()[i] += 1.0f;
      }
    }
    shape[1] = nchan;
    xt::xtensor<std::complex<float>, 3> result_data(shape, 0.0f);
    xt::xtensor<float, 3> result_weights(shape, 0.0f);
    xt::xtensor<bool, 3> result_flags(shape, itsFlag);
    // Average to get the true expected result.
    for (size_t bl = 0; bl < itsNBl; ++bl) {
      for (size_t corr = 0; corr < itsNCorr; ++corr) {
        for (size_t ch = 0; ch < nchan; ++ch) {
          size_t avg_chan;
          for (avg_chan = ch * itsNAvgChan;
               avg_chan < std::min((ch + 1) * itsNAvgChan, itsNChan);
               ++avg_chan) {
            result_data(bl, ch, corr) += data(bl, avg_chan, corr);
            result_weights(bl, ch, corr) += weights(bl, avg_chan, corr);
          }
          result_data(bl, ch, corr) /=
              float(navgtime * (avg_chan - ch * itsNAvgChan));
        }
      }
    }
    // Check the averaged result.
    BOOST_CHECK(xt::allclose(buffer->GetData(), result_data));
    BOOST_CHECK(buffer->GetFlags() == result_flags);
    BOOST_CHECK_CLOSE(
        buffer->GetTime(),
        2 + 5 * (itsCount * itsNAvgTime + (itsNAvgTime - 1) / 2.0), 1.0e-3);
    BOOST_CHECK(xt::allclose(buffer->GetWeights(), result_weights));
    if (navgtime == itsNAvgTime) {
      const xt::xtensor<double, 2> uvw =
          (100 * (itsCount * itsNAvgTime + 0.5 * (itsNAvgTime - 1))) +
          xt::arange(itsNBl * 3).reshape({itsNBl, 3});
      BOOST_CHECK(xt::allclose(buffer->GetUvw(), uvw));
    }

    ++itsCount;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(itsNChan, info.origNChan());
    BOOST_CHECK_EQUAL(1 + (itsNChan - 1) / itsNAvgChan, info.nchan());
    BOOST_CHECK_EQUAL(1 + (itsNTime - 1) / itsNAvgTime, info.ntime());
    BOOST_CHECK_EQUAL(5 * itsNAvgTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNAvgChan, info.nchanAvg());
    BOOST_CHECK_EQUAL(itsNAvgTime, info.ntimeAvg());
  }

 private:
  size_t itsCount;
  size_t itsNTime;
  size_t itsNBl;
  size_t itsNChan;
  size_t itsNCorr;
  size_t itsNAvgTime;
  size_t itsNAvgChan;
  bool itsFlag;
};

// More elaborate class which can set different flags and weights.
class TestInput3 : public dp3::steps::MockInput {
 public:
  TestInput3(int nrtime, int nrbl, int nrchan, int nrcorr)
      : itsCount(0),
        itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (itsCount == itsNrTime) {
      return false;
    }
    const std::array<std::size_t, 3> shape{itsNrBl, itsNrChan, itsNrCorr};
    buffer->GetData().resize(shape);
    buffer->GetWeights().resize(shape);
    buffer->GetFlags().resize(shape);
    int i = 0;
    for (std::size_t bl = 0; bl < itsNrBl; ++bl) {
      for (std::size_t ch = 0; ch < itsNrChan; ++ch) {
        for (std::size_t corr = 0; corr < itsNrCorr; ++corr) {
          buffer->GetData()(bl, ch, corr) =
              std::complex<float>(i + itsCount * 10, i - 1000 + itsCount * 6);
          buffer->GetWeights()(bl, ch, corr) =
              (1 + (itsCount + bl + ch) % 5) / 5.;
          buffer->GetFlags()(bl, ch, corr) =
              ((itsCount + 2 * bl + 3 * ch) % 7 == 0);
          i++;
        }
      }
    }
    buffer->SetTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    casacore::Vector<dp3::common::rownr_t> rownrs(1, itsCount);
    buffer->SetRowNumbers(rownrs);
    buffer->GetUvw().resize({itsNrBl, 3});
    xt::flatten(buffer->GetUvw()) = xt::arange(buffer->GetUvw().size());
    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Use timeInterval=5
    GetWritableInfoOut() = DPInfo(itsNrCorr, itsNrChan);
    GetWritableInfoOut().setTimes(100.0, 100.0 + (itsNrTime - 1) * 5.0, 5.0);
    // Define the frequencies.
    std::vector<double> chanFreqs;
    std::vector<double> chanWidth(itsNrChan, 100000.);
    for (std::size_t i = 0; i < itsNrChan; i++) {
      chanFreqs.push_back(1050000. + i * 100000.);
    }
    GetWritableInfoOut().setChannels(std::move(chanFreqs),
                                     std::move(chanWidth));
  }

 private:
  int itsCount;
  int itsNrTime;
  std::size_t itsNrBl;
  std::size_t itsNrChan;
  std::size_t itsNrCorr;
};

// Class to check result of averaging TestInput3.
// All input must be averaged (in one or more steps) to a single value
// per corr/baseline.
class TestOutput3 : public dp3::steps::test::ThrowStep {
 public:
  TestOutput3(int nrtime, int nrbl, int nrchan, int nrcorr)
      : itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    xt::xtensor<std::complex<float>, 3> result({itsNrBl, 1, itsNrCorr}, 0.0f);
    xt::xtensor<float, 3> weights({itsNrBl, 1, itsNrCorr}, 0.0f);
    xt::xtensor<bool, 3> flags({itsNrBl, 1, itsNrCorr}, true);
    // Create data in the same way as in TestInput3.
    for (size_t time = 0; time < itsNrTime; ++time) {
      int i = 0;
      for (size_t bl = 0; bl < itsNrBl; ++bl) {
        for (size_t ch = 0; ch < itsNrChan; ++ch) {
          for (size_t corr = 0; corr < itsNrCorr; ++corr) {
            if ((time + 2 * bl + 3 * ch) % 7 != 0) {
              float weight = (1 + (time + bl + ch) % 5) / 5.;
              result(bl, 0, corr) +=
                  weight *
                  std::complex<float>(i + time * 10, i - 1000 + int(time) * 6);
              weights(bl, 0, corr) += weight;
              flags(bl, 0, corr) = false;
            }
            i++;
          }
        }
      }
    }
    BOOST_CHECK(xt::all(xt::not_equal(weights, 0.0f)));
    for (unsigned int i = 0; i < result.size(); ++i) {
      result.data()[i] /= weights.data()[i];
    }
    // Check the averaged result.
    BOOST_CHECK(xt::allclose(buffer->GetData(), result));
    BOOST_CHECK(buffer->GetFlags() == flags);
    BOOST_CHECK_CLOSE(buffer->GetTime(), 2.0 + 5 * (itsNrTime - 1) / 2.0,
                      1.0e-3);
    BOOST_CHECK(xt::allclose(buffer->GetWeights(), weights));
    const xt::xtensor<double, 2> uvw =
        xt::arange(itsNrBl * 3).reshape({itsNrBl, 3});
    BOOST_CHECK(xt::allclose(buffer->GetUvw(), uvw));
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(itsNrChan, info.origNChan());
    BOOST_CHECK_EQUAL(size_t{1}, info.nchan());
    BOOST_CHECK_EQUAL(size_t{1}, info.ntime());
    BOOST_CHECK_EQUAL(5 * itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNrChan, info.nchanAvg());
    BOOST_CHECK_EQUAL(itsNrTime, info.ntimeAvg());
  }

 private:
  size_t itsNrTime;
  size_t itsNrBl;
  size_t itsNrChan;
  size_t itsNrCorr;
};

// Simple class to flag every step-th XX point.
class TestFlagger : public dp3::steps::test::ThrowStep {
 public:
  TestFlagger(int step) : itsCount(0), itsStep(step) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    int ncorr = buffer->GetFlags().shape(2);
    int np = buffer->GetFlags().size() / ncorr;
    bool* flagPtr = buffer->GetFlags().data();
    for (int i = 0; i < np; ++i) {
      if ((i + itsCount) % itsStep == 0) {
        for (int j = 0; j < ncorr; ++j) {
          flagPtr[i * ncorr + j] = true;
        }
      }
    }
    getNextStep()->process(std::move(buffer));
    ++itsCount;
    return true;
  }

  void updateInfo(const DPInfo& info) override { Step::updateInfo(info); }
  void finish() override { getNextStep()->finish(); }

 private:
  int itsCount;
  int itsStep;
};

// Class to check result of averaging and flagging TestInput3.
// First the data are averaged from 8,4 to 4,2, then every step-th point
// is flagged, and finally it is averaged to 1,1.
class TestOutput4 : public dp3::steps::test::ThrowStep {
 public:
  TestOutput4(int nrtime, int nrbl, int nrchan, int nrcorr, int step)
      : itsNrTime(nrtime),
        itsNrBl(nrbl),
        itsNrChan(nrchan),
        itsNrCorr(nrcorr),
        itsStep(step) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    xt::xtensor<std::complex<float>, 3> result({itsNrBl, 1, itsNrCorr}, 0.0f);
    xt::xtensor<float, 3> weights({itsNrBl, 1, itsNrCorr}, 0.0f);
    xt::xtensor<bool, 3> flags({itsNrBl, 1, itsNrCorr}, true);
    // Create data in the same way as in TestInput3.
    for (size_t time = 0; time < itsNrTime; ++time) {
      int i = 0;
      for (size_t bl = 0; bl < itsNrBl; ++bl) {
        for (size_t ch = 0; ch < itsNrChan; ++ch) {
          // TestFlagger flags every step-th point of 2x2 averaged data.
          size_t tf = time / 2;  // same as itsCount in testFlagger
          if (((bl * itsNrChan + ch) / 2 + tf) % itsStep == 0) {
            i += itsNrCorr;
          } else {
            for (size_t corr = 0; corr < itsNrCorr; ++corr) {
              if ((time + 2 * bl + 3 * ch) % 7 != 0) {
                float weight = (1 + (time + bl + ch) % 5) / 5.;
                result(bl, 0, corr) +=
                    weight * std::complex<float>(i + time * 10,
                                                 i - 1000 + int(time) * 6);
                weights(bl, 0, corr) += weight;
                flags(bl, 0, corr) = false;
              }
              i++;
            }
          }
        }
      }
    }
    for (unsigned int i = 0; i < result.size(); ++i) {
      if (!flags.data()[i]) {
        result.data()[i] /= weights.data()[i];
      }
    }
    // Check the averaged result.
    BOOST_CHECK(xt::allclose(buffer->GetData(), result));
    BOOST_CHECK(buffer->GetFlags() == flags);
    BOOST_CHECK_CLOSE(buffer->GetTime(), 2.0 + 5 * (itsNrTime - 1) / 2.0,
                      1.0e-3);
    BOOST_CHECK(xt::allclose(buffer->GetWeights(), weights));
    const xt::xtensor<double, 2> uvw =
        xt::arange(itsNrBl * 3).reshape({itsNrBl, 3});
    BOOST_CHECK(xt::allclose(buffer->GetUvw(), uvw));
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(itsNrChan, info.origNChan());
    BOOST_CHECK_EQUAL(size_t{1}, info.nchan());
    BOOST_CHECK_EQUAL(size_t{1}, info.ntime());
    BOOST_CHECK_EQUAL(5 * itsNrTime, info.timeInterval());
    BOOST_CHECK_EQUAL(itsNrChan, info.nchanAvg());
    BOOST_CHECK_EQUAL(itsNrTime, info.ntimeAvg());
  }

  const size_t itsNrTime;
  const size_t itsNrBl;
  const size_t itsNrChan;
  const size_t itsNrCorr;
  const size_t itsStep;
};

// Test simple averaging without flagged points.
void test1(int ntime, int nbl, int nchan, int ncorr, int navgtime, int navgchan,
           bool flag) {
  // Create the steps.
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqstep", std::to_string(navgchan));
  parset.add("timestep", std::to_string(navgtime));
  auto step2 = std::make_shared<Averager>(parset, "");
  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, navgtime,
                                            navgchan, flag);
  dp3::steps::test::Execute({step1, step2, step3});
}

// Like test 1, but specify target resolution
void test1resolution(int ntime, int nbl, int nchan, int ncorr,
                     double timeresolution, double freqresolution,
                     string frequnit, bool flag) {
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset;
  parset.add("freqresolution", std::to_string(freqresolution) + frequnit);
  parset.add("timeresolution", std::to_string(timeresolution));
  auto step2 = std::make_shared<Averager>(parset, "");

  if (!frequnit.empty()) {
    casacore::Quantity q(freqresolution, frequnit);
    freqresolution = q.getValue("Hz", true);
  }

  int navgchan = std::max(1, int(freqresolution / 100000 + 0.5));
  int navgtime = std::max(1, int(timeresolution / 5. + 0.5));

  auto step3 = std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, navgtime,
                                            navgchan, flag);
  dp3::steps::test::Execute({step1, step2, step3});
}

// Like test1, but the averaging is done in two steps.
void test2(int ntime, int nbl, int nchan, int ncorr, bool flag) {
  // Create the steps.
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, ncorr, flag);
  ParameterSet parset1, parset2;
  parset1.add("freqstep", "4");
  parset2.add("timestep", "2");
  auto step2a = std::make_shared<Averager>(parset1, "");
  auto step2b = std::make_shared<Averager>(parset2, "");
  auto step3 =
      std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr, 2, 4, flag);
  dp3::steps::test::Execute({step1, step2a, step2b, step3});
}

// Do tests with weighting and some flagged points.
void test3(int nrbl, int nrcorr) {
  {
    // Create the steps.
    auto step1 = std::make_shared<TestInput3>(2, nrbl, 2, nrcorr);
    ParameterSet parset1;
    parset1.add("freqstep", "2");
    parset1.add("timestep", "2");
    auto step2 = std::make_shared<Averager>(parset1, "");
    auto step3 = std::make_shared<TestOutput3>(2, nrbl, 2, nrcorr);
    dp3::steps::test::Execute({step1, step2, step3});
  }
  {
    // Create the steps.
    auto step1 = std::make_shared<TestInput3>(4, nrbl, 8, nrcorr);
    ParameterSet parset1, parset2;
    parset1.add("freqstep", "4");
    parset1.add("timestep", "2");
    parset2.add("freqstep", "2");
    parset2.add("timestep", "2");
    auto step2a = std::make_shared<Averager>(parset1, "");
    auto step2b = std::make_shared<Averager>(parset2, "");
    auto step3 = std::make_shared<TestOutput3>(4, nrbl, 8, nrcorr);
    dp3::steps::test::Execute({step1, step2a, step2b, step3});
  }
}

void TestAveragingAndFlagging(int nrbl, int nrcorr, int flagstep) {
  auto step1 = std::make_shared<TestInput3>(4, nrbl, 8, nrcorr);
  ParameterSet parset1, parset2;
  parset1.add("freqstep", "2");
  parset1.add("timestep", "2");
  parset2.add("freqstep", "4");
  parset2.add("timestep", "2");
  auto step2a = std::make_shared<Averager>(parset1, "");
  auto step2b = std::make_shared<TestFlagger>(flagstep);
  auto step2c = std::make_shared<Averager>(parset2, "");
  auto step3 = std::make_shared<TestOutput4>(4, nrbl, 8, nrcorr, flagstep);
  dp3::steps::test::Execute({step1, step2a, step2b, step2c, step3});
}

BOOST_AUTO_TEST_CASE(testaverager1) { test1(10, 3, 32, 4, 2, 4, false); }

BOOST_AUTO_TEST_CASE(testaverager2) { test1(10, 3, 30, 1, 3, 3, true); }

BOOST_AUTO_TEST_CASE(testaverager3) { test1(10, 3, 30, 1, 3, 3, false); }

BOOST_AUTO_TEST_CASE(testaverager4) { test1(11, 3, 30, 2, 3, 3, false); }

BOOST_AUTO_TEST_CASE(testaverager5) { test1(10, 3, 32, 4, 1, 32, false); }

BOOST_AUTO_TEST_CASE(testaverager6) { test1(10, 3, 32, 1, 1, 1, false); }

BOOST_AUTO_TEST_CASE(testaverager7) { test2(10, 3, 32, 2, true); }

BOOST_AUTO_TEST_CASE(testaverager8) { test2(10, 3, 32, 2, false); }

BOOST_AUTO_TEST_CASE(testaverager9) { test3(1, 1); }

BOOST_AUTO_TEST_CASE(testaverager10) { test3(10, 4); }

BOOST_AUTO_TEST_CASE(testaverager11) { TestAveragingAndFlagging(1, 4, 3); }

BOOST_AUTO_TEST_CASE(testaverager12) { TestAveragingAndFlagging(20, 4, 5); }

BOOST_AUTO_TEST_CASE(testresolution1) {
  test1resolution(10, 3, 32, 4, 10., 100000, "Hz", false);
}

BOOST_AUTO_TEST_CASE(testresolution2) {
  test1resolution(11, 3, 32, 4, 1., 800, "kHz", false);
}

BOOST_AUTO_TEST_CASE(testresolution3) {
  test1resolution(11, 3, 32, 4, 15., 0.4, "MHz", false);
}

BOOST_AUTO_TEST_SUITE_END()
