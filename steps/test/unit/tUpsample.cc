// tAverager.cc: Test program for class Averager
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Upsample.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../base/UVWCalculator.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Upsample;

namespace {
const std::size_t kNCorr = 4;
const std::size_t kNChannels = 5;
const std::size_t kNBaselines = 3;

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::InputStep {
 public:
  TestInput(const std::vector<double>& times, const std::vector<bool>& flags,
            const std::vector<double>& uvws, double time_interval)
      : time_step_(0),
        times_(times),
        flags_(flags),
        uvws_(uvws),
        time_interval_(time_interval) {}

  bool process(const DPBuffer&) override {
    // Stop when all times are done.
    if (time_step_ == times_.size()) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(kNCorr, kNChannels, kNBaselines);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] =
          casacore::Complex(i + time_step_ * 10, i - 1000 + time_step_ * 6);
    }
    DPBuffer buf;
    buf.setTime(times_[time_step_]);
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Cube<bool> flags(data.shape());
    flags = flags_[time_step_];
    buf.setFlags(flags);
    buf.setExposure(time_interval_);

    if (!uvws_.empty()) {
      casacore::Matrix<double> uvw(3, kNBaselines);
      indgen(uvw, uvws_[time_step_]);
      buf.setUVW(uvw);
    }
    getNextStep()->process(buf);
    ++time_step_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }

  void show(std::ostream&) const override {}

  void updateInfo(const DPInfo&) override {
    info().init(kNCorr, 0, kNChannels, times_.size(), times_[0], time_interval_,
                "", "");
    // Define the frequencies.
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(kNChannels, 100000.);
    for (unsigned int i = 0; i < kNChannels; i++) {
      chan_freqs.push_back(1050000. + i * 100000.);
    }
    info().set(std::move(chan_freqs), std::move(chan_width));

    // Define antennas and baselines.
    const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
    const std::vector<casacore::MPosition> kAntPos{
        casacore::MVPosition{100, 0, 0}, casacore::MVPosition{200, 0, 0},
        casacore::MVPosition{300, 0, 0}};
    const std::vector<double> kAntDiam(kNBaselines, 1.0);
    // Baseline 0 has auto-correlations.
    const std::vector<int> kAnt1{0, 0, 1};
    const std::vector<int> kAnt2{0, 2, 2};
    info().set(kAntNames, kAntDiam, kAntPos, kAnt1, kAnt2);
  }

 private:
  unsigned int time_step_;
  const std::vector<double> times_;
  const std::vector<bool> flags_;
  const std::vector<double> uvws_;
  const double time_interval_;
};

// Class to check result of upsampling TestInput
class TestOutput : public dp3::steps::Step {
 public:
  TestOutput(const std::vector<double>& times, const std::vector<bool>& flags,
             const std::vector<double>& uvws, double time_interval,
             bool update_uvw)
      : time_step_(0),
        times_(times),
        flags_(flags),
        uvws_(uvws),
        time_interval_(time_interval),
        update_uvw_(update_uvw) {}

  bool process(const DPBuffer& buf) override {
    BOOST_CHECK_SMALL(buf.getTime() - times_[time_step_],
                      time_interval_ * 0.01);
    BOOST_CHECK(allTrue(buf.getFlags()) == flags_[time_step_]);

    const casacore::Matrix<double>& buf_uvw = buf.getUVW();
    BOOST_REQUIRE(buf_uvw.shape() == casacore::IPosition(2, 3, kNBaselines));
    if (update_uvw_) {
      // The first baseline has auto-correlations. UVW should be zero.
      BOOST_CHECK_CLOSE(buf_uvw(0, 0), 0.0, 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(1, 0), 0.0, 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(2, 0), 0.0, 1.0e-6);
      // For the other baselines, check that Upsample set the UVW correctly.
      dp3::base::UVWCalculator calc(info().phaseCenter(), info().arrayPos(),
                                    info().antennaPos());
      const std::array<double, 3> uvw_0_2 = calc.getUVW(0, 2, buf.getTime());
      BOOST_CHECK_CLOSE(buf_uvw(0, 1), uvw_0_2[0], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(1, 1), uvw_0_2[1], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(2, 1), uvw_0_2[2], 1.0e-6);
      const std::array<double, 3> uvw_1_2 = calc.getUVW(1, 2, buf.getTime());
      BOOST_CHECK_CLOSE(buf_uvw(0, 2), uvw_1_2[0], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(1, 2), uvw_1_2[1], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(2, 2), uvw_1_2[2], 1.0e-6);
    } else {
      // Upsample should have copied the original UVW values.
      casacore::Matrix<double> expected_uvw(3, kNBaselines);
      indgen(expected_uvw, uvws_[time_step_]);
      BOOST_CHECK_EQUAL_COLLECTIONS(buf_uvw.begin(), buf_uvw.end(),
                                    expected_uvw.begin(), expected_uvw.end());
    }

    ++time_step_;
    return true;
  }

  void finish() override {}
  void show(std::ostream&) const override {}
  void updateInfo(const DPInfo& info) override {
    Step::updateInfo(info);
    BOOST_CHECK(casacore::near(info.timeInterval(), time_interval_));
  }

 private:
  unsigned int time_step_;
  const std::vector<double> times_;
  const std::vector<bool> flags_;
  const std::vector<double> uvws_;
  const double time_interval_;
  const bool update_uvw_;
};

void TestUpsample(bool update_uvw) {
  const double kTimeInterval = 2.01327;

  const std::vector<double> times{5020763030.74, 5020763032.75, 5020763034.76,
                                  5020763035.77, 5020763037.78, 5020763039.8};
  const std::vector<bool> flags{false, false, true, false, false, false};
  std::vector<double> uvws{4100, 4200, 4300, 4400, 4500, 4600};
  if (update_uvw) {
    // Test that Upsample sets UVW values even if the input has no UVW values.
    uvws.clear();
  }
  auto in_step = std::make_shared<TestInput>(times, flags, uvws, kTimeInterval);

  ParameterSet parset;
  parset.add("timestep", "2");
  if (update_uvw) parset.add("updateuvw", "true");
  auto upsample = std::make_shared<Upsample>(parset, "");

  const std::vector<double> new_times{
      5020763030.23, 5020763031.24, 5020763032.25, 5020763033.25,
      5020763034.26, 5020763035.27, 5020763036.27, 5020763037.28,
      5020763038.29, 5020763039.29, 5020763040.3};
  const std::vector<bool> new_flags{false, false, false, false, true, false,
                                    false, false, false, false, false};
  // new_uvws are only used when update_uvw is false.
  const std::vector<double> new_uvws{4100, 4100, 4200, 4200, 4300, 4400,
                                     4400, 4500, 4500, 4600, 4600};
  auto out_step = std::make_shared<TestOutput>(new_times, new_flags, new_uvws,
                                               0.5 * kTimeInterval, update_uvw);

  in_step->setNextStep(upsample);
  upsample->setNextStep(out_step);

  in_step->setInfo(DPInfo());

  DPBuffer buf;
  while (in_step->process(buf)) {
  }
  in_step->finish();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(upsample)

BOOST_AUTO_TEST_CASE(copy_uvw) { TestUpsample(false); }

BOOST_AUTO_TEST_CASE(update_uvw) { TestUpsample(true); }

BOOST_AUTO_TEST_SUITE_END()
