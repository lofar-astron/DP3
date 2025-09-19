// tAverager.cc: Test program for class Averager
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Upsample.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Quanta/Quantum.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "tStepCommon.h"
#include "mock/MockInput.h"
#include "mock/ThrowStep.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../base/UVWCalculator.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

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
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(const std::vector<double>& times, const std::vector<bool>& flags,
            const std::vector<double>& uvws, double time_interval)
      : time_step_(0),
        times_(times),
        flags_(flags),
        uvws_(uvws),
        time_interval_(time_interval) {}

  bool process(std::unique_ptr<DPBuffer>) override {
    // Stop when all times are done.
    if (time_step_ == times_.size()) {
      return false;
    }
    std::array<size_t, 3> data_shape{kNBaselines, kNChannels, kNCorr};
    auto buffer =
        std::make_unique<DPBuffer>(times_[time_step_], time_interval_);
    buffer->GetData().resize(data_shape);
    for (std::size_t i = 0; i < buffer->GetData().size(); ++i) {
      buffer->GetData().data()[i] = std::complex<float>(
          i + time_step_ * 10.0,
          static_cast<int>(i) - 1000 + static_cast<int>(time_step_) * 6.0);
    }
    buffer->GetWeights().resize(data_shape);
    buffer->GetWeights().fill(1.0);
    buffer->GetFlags().resize(data_shape);
    buffer->GetFlags().fill(flags_[time_step_]);

    if (!uvws_.empty()) {
      buffer->GetUvw().resize({kNBaselines, 3});
      buffer->GetUvw() = xt::arange(uvws_[time_step_],
                                    uvws_[time_step_] + (kNBaselines * 3), 1)
                             .reshape(buffer->GetUvw().shape());
    }
    getNextStep()->process(std::move(buffer));
    ++time_step_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }

  void show(std::ostream&) const override {}

  void updateInfo(const DPInfo&) override {
    GetWritableInfoOut() = DPInfo(kNCorr, kNChannels);
    GetWritableInfoOut().setTimes(times_.front(), times_.back(),
                                  time_interval_);
    // Define the frequencies.
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(kNChannels, 100000.0);
    for (std::size_t i = 0; i < kNChannels; i++) {
      chan_freqs.push_back(1050000.0 + i * 100000.0);
    }
    GetWritableInfoOut().setChannels(std::move(chan_freqs),
                                     std::move(chan_width));

    // Define antennas and baselines.
    const std::vector<std::string> kAntNames{"ant0", "ant1", "ant2"};
    const std::vector<casacore::MPosition> kAntPos{
        casacore::MVPosition{100, 0, 0}, casacore::MVPosition{200, 0, 0},
        casacore::MVPosition{300, 0, 0}};
    const std::vector<double> kAntDiam(kNBaselines, 1.0);
    // Baseline 0 has auto-correlations.
    const std::vector<int> kAnt1{0, 0, 1};
    const std::vector<int> kAnt2{0, 2, 2};
    GetWritableInfoOut().setAntennas(kAntNames, kAntDiam, kAntPos, kAnt1,
                                     kAnt2);
  }

 private:
  unsigned int time_step_;
  const std::vector<double> times_;
  const std::vector<bool> flags_;
  const std::vector<double> uvws_;
  const double time_interval_;
};

// Class to check result of upsampling TestInput
class TestOutput : public dp3::steps::test::ThrowStep {
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

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    BOOST_CHECK_SMALL(buffer->GetTime() - times_[time_step_],
                      time_interval_ * 0.01);
    BOOST_CHECK(xt::all(xt::equal(buffer->GetFlags(), flags_[time_step_])));

    DPBuffer::UvwType buf_uvw = buffer->GetUvw();
    BOOST_TEST(buf_uvw.shape() == (std::array<std::size_t, 2>{kNBaselines, 3}),
               boost::test_tools::per_element());
    if (update_uvw_) {
      // The first baseline has auto-correlations. UVW should be zero.
      BOOST_CHECK_CLOSE(buf_uvw(0, 0), 0.0, 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(0, 1), 0.0, 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(0, 2), 0.0, 1.0e-6);
      // For the other baselines, check that Upsample set the UVW correctly.
      dp3::base::UVWCalculator calc(getInfoOut().phaseCenter(),
                                    getInfoOut().arrayPos(),
                                    getInfoOut().antennaPos());
      const std::array<double, 3> uvw_0_2 =
          calc.getUVW(0, 2, buffer->GetTime());
      BOOST_CHECK_CLOSE(buf_uvw(1, 0), uvw_0_2[0], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(1, 1), uvw_0_2[1], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(1, 2), uvw_0_2[2], 1.0e-6);
      const std::array<double, 3> uvw_1_2 =
          calc.getUVW(1, 2, buffer->GetTime());
      BOOST_CHECK_CLOSE(buf_uvw(2, 0), uvw_1_2[0], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(2, 1), uvw_1_2[1], 1.0e-6);
      BOOST_CHECK_CLOSE(buf_uvw(2, 2), uvw_1_2[2], 1.0e-6);
    } else {
      // Upsample should have copied the original UVW values.
      xt::xtensor<double, 2> expected_uvw =
          xt::arange(uvws_[time_step_], uvws_[time_step_] + (kNBaselines * 3),
                     1)
              .reshape(buf_uvw.shape());
      BOOST_CHECK_EQUAL_COLLECTIONS(buf_uvw.begin(), buf_uvw.end(),
                                    expected_uvw.begin(), expected_uvw.end());
    }

    ++time_step_;
    return true;
  }

  void finish() override {}
  void updateInfo(const DPInfo& info) override {
    Step::updateInfo(info);
    BOOST_CHECK_CLOSE(info.timeInterval(), time_interval_, 1.0e-3);
  }

 private:
  unsigned int time_step_;
  const std::vector<double> times_;
  const std::vector<bool> flags_;
  const std::vector<double> uvws_;
  const double time_interval_;
  const bool update_uvw_;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(upsample)

BOOST_AUTO_TEST_CASE(fields) {
  using dp3::steps::Step;

  const Upsample no_update_uvw("no_update_uvw", 2, false);
  BOOST_TEST(no_update_uvw.getRequiredFields() == Step::kFlagsField);
  BOOST_TEST(no_update_uvw.getProvidedFields() == dp3::common::Fields());

  const Upsample updates_uvw("updates_uvw", 2, true);
  BOOST_TEST(updates_uvw.getRequiredFields() ==
             (Step::kFlagsField | Step::kUvwField));
  BOOST_TEST(updates_uvw.getProvidedFields() == Step::kUvwField);
}

BOOST_DATA_TEST_CASE(execute, boost::unit_test::data::make({true, false}),
                     update_uvw) {
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

  dp3::steps::test::Execute({in_step, upsample, out_step});
}

BOOST_AUTO_TEST_SUITE_END()
