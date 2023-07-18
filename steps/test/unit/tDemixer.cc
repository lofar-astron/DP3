// tDemixer.cc: Test program for class Demixer
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Demixer.h"

#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <string>

#include <boost/test/unit_test.hpp>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>

#include "tPredict.h"
#include "tStepCommon.h"
#include "mock/ThrowStep.h"

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

#include "../../Averager.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Demixer;
using dp3::steps::Step;

namespace {

// Demixer only works with 4 correlations
const size_t kNCorr = 4;

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, baselines.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(size_t ntime, size_t nbl, size_t nchan, bool flag)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        flag_data_(flag) {}

  bool process(std::unique_ptr<DPBuffer> buffer) override {
    // Stop when all times are done.
    if (count_ == n_times_) {
      return false;
    }

    buffer->setTime(count_ * 5 + 2);  // same interval as in updateAverageInfo

    const std::array<size_t, 3> shape{n_baselines_, n_channels_, kNCorr};
    buffer->GetData().resize(shape);
    for (int i = 0; i < static_cast<int>(buffer->GetData().size()); ++i) {
      buffer->GetData().data()[i] =
          std::complex<float>(i + count_ * 10.0f, i - 1000.0f + count_ * 6.0f);
    }

    buffer->GetWeights().resize(shape);
    buffer->GetWeights().fill(1.0);

    buffer->ResizeFlags(shape);
    buffer->GetFlags().fill(flag_data_);

    buffer->ResizeUvw(n_baselines_);
    buffer->GetUvw().fill(double(count_ * 100));
    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      buffer->GetUvw()(bl, 0) += bl * 3;
      buffer->GetUvw()(bl, 1) += bl * 3 + 1;
      buffer->GetUvw()(bl, 2) += bl * 3 + 2;
    }

    getNextStep()->process(std::move(buffer));
    ++count_;
    return true;
  }

  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override
  // Use startchan=8 and timeInterval=5
  {
    info() = DPInfo(kNCorr, n_channels_);
    info().setTimes(2.5, 2.5 + (n_times_ - 1) * 5.0, 5.0);

    std::vector<std::string> ant_names = {"antenna_1", "antenna_2"};
    std::vector<double> ant_diameter = {0.0, 1.0};
    const std::vector<casacore::MPosition> ant_position{
        casacore::MVPosition{0, 0, 0}, casacore::MVPosition{1, 1, 1}};
    std::vector<int> antenna1 = {0, 0, 1};
    std::vector<int> antenna2 = {0, 1, 1};
    std::vector<unsigned int> baselines = {0, 1, 2};
    info().setAntennas(ant_names, ant_diameter, ant_position, antenna1,
                       antenna2);

    double start_channel = 8.0;
    std::vector<double> chan_freqs;
    std::vector<double> chan_width(n_channels_, 100000.);
    for (size_t i = 0; i < n_channels_; i++) {
      chan_freqs.push_back(start_channel + i * 100000.);
    }
    info().setChannels(std::move(chan_freqs), std::move(chan_width));
    info().update(start_channel, n_channels_, baselines, false);
  }

 private:
  size_t count_;
  size_t n_times_;
  size_t n_baselines_;
  size_t n_channels_;
  bool flag_data_;
};

// Class to check result of averaging TestInput.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, size_t nbl, size_t nchan, size_t navgtime,
             size_t navgchan, bool flag)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        n_average_time_(navgtime),
        n_average_channel_(navgchan),
        flag_data_(flag) {}

 private:
  bool process(std::unique_ptr<DPBuffer> buffer) override {
    size_t nchan = 1 + (n_channels_ - 1) / n_average_channel_;
    assert(n_times_ >= count_ * n_average_time_);
    size_t navgtime =
        std::min(n_average_time_, n_times_ - count_ * n_average_time_);

    // Fill expected result in similar way as TestInput.
    std::array<size_t, 3> shape{n_baselines_, n_channels_, kNCorr};
    xt::xtensor<std::complex<float>, 3> data(shape, 0.0f);
    xt::xtensor<float, 3> weights(shape, 0.0f);
    if (!flag_data_) {
      for (size_t j = count_ * n_average_time_;
           j < count_ * n_average_time_ + navgtime; ++j) {
        for (size_t i = 0; i < data.size(); ++i) {
          data.data()[i] += std::complex<float>(
              i + j * 10, static_cast<int>(i) - 1000 + static_cast<int>(j) * 6);
          weights.data()[i] += 1.0f;
        }
      }
    }

    shape = {n_baselines_, nchan, kNCorr};
    xt::xtensor<std::complex<float>, 3> result_data(shape, 0.0f);
    xt::xtensor<float, 3> result_weights(shape, 0.0f);
    // Average to get the true expected result_data.
    for (size_t bl = 0; bl < n_baselines_; ++bl) {
      for (size_t corr = 0; corr < kNCorr; ++corr) {
        for (size_t ch = 0; ch < nchan; ++ch) {
          size_t avg_ch = ch * n_average_channel_;
          for (; avg_ch < std::min((ch + 1) * n_average_channel_, n_channels_);
               ++avg_ch) {
            result_data(bl, ch, corr) += data(bl, avg_ch, corr);
            result_weights(bl, ch, corr) += weights(bl, avg_ch, corr);
          }
          result_data(bl, ch, corr) /=
              float(navgtime * (avg_ch - ch * n_average_channel_));
        }
      }
    }

    // Check the averaged result_data. When all the flags are set, do not check
    // the values of the data and weights, since DP3 should not use those values
    // anyway.
    if (!flag_data_) {
      BOOST_CHECK(xt::allclose(buffer->GetData(), result_data));
      BOOST_CHECK(xt::allclose(buffer->GetWeights(), result_weights));
    }
    if (navgtime == n_average_time_) {
      xt::xtensor<double, 2> result_uvw(
          {n_baselines_, 3},
          100 * (count_ * n_average_time_ + 0.5 * (n_average_time_ - 1)));
      for (size_t bl = 0; bl < n_baselines_; ++bl) {
        result_uvw(bl, 0) += bl * 3;
        result_uvw(bl, 1) += bl * 3 + 1;
        result_uvw(bl, 2) += bl * 3 + 2;
      }
      BOOST_CHECK(xt::allclose(buffer->GetUvw(), result_uvw));
    }
    xt::xtensor<bool, 3> result_flags(shape, flag_data_);
    BOOST_CHECK_EQUAL(buffer->GetFlags(), result_flags);
    BOOST_CHECK_CLOSE(
        buffer->getTime(),
        2.0 + 5.0 * (count_ * n_average_time_ + (n_average_time_ - 1) / 2.0),
        1.0e-3);

    ++count_;
    return true;
  }

  void finish() override {}
  virtual void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(size_t{8}, info.startchan());
    BOOST_CHECK_EQUAL(n_channels_, info.origNChan());
    BOOST_CHECK_EQUAL(1 + (n_channels_ - 1) / n_average_channel_, info.nchan());
    BOOST_CHECK_EQUAL(1 + (n_times_ - 1) / n_average_time_, info.ntime());
    BOOST_CHECK_EQUAL(5 * n_average_time_, info.timeInterval());
    BOOST_CHECK_EQUAL(n_average_channel_, info.nchanAvg());
    BOOST_CHECK_EQUAL(n_average_time_, info.ntimeAvg());
  }

  int count_;
  int n_times_;
  size_t n_baselines_;
  size_t n_channels_;
  size_t n_average_time_;
  size_t n_average_channel_;
  bool flag_data_;
};

// This test only tests the averager functionality of the Demixer.
void TestDemixer(int ntime, int nbl, int nchan, int navgtime, int navgchan,
                 bool flag) {
  auto step1 = std::make_shared<TestInput>(ntime, nbl, nchan, flag);
  ParameterSet parset;
  parset.add("freqstep", std::to_string(navgchan));
  parset.add("timestep", std::to_string(navgtime));
  parset.add("skymodel", dp3::steps::test::kPredictSourceDB);
  auto step2 = std::make_shared<Demixer>(parset, "");
  auto step3 =
      std::make_shared<TestOutput>(ntime, nbl, nchan, navgtime, navgchan, flag);
  dp3::steps::test::Execute({step1, step2, step3});
}

}  // namespace

BOOST_AUTO_TEST_SUITE(demixer)

BOOST_AUTO_TEST_CASE(execute_1) { TestDemixer(10, 3, 32, 2, 4, false); }

BOOST_AUTO_TEST_CASE(execute_2) { TestDemixer(10, 3, 30, 3, 3, true); }

BOOST_AUTO_TEST_CASE(execute_3) { TestDemixer(10, 3, 30, 3, 3, false); }

BOOST_AUTO_TEST_CASE(execute_4) { TestDemixer(11, 3, 30, 3, 3, false); }

BOOST_AUTO_TEST_CASE(execute_5) { TestDemixer(10, 3, 32, 1, 32, false); }

BOOST_AUTO_TEST_CASE(fields) {
  using dp3::steps::Averager;

  ParameterSet parset;
  parset.add("skymodel", dp3::steps::test::kPredictSourceDB);
  Demixer demixer(parset, "");

  BOOST_TEST(demixer.getRequiredFields() == Averager::kRequiredFields);
  BOOST_TEST(demixer.getProvidedFields() == Averager::kProvidedFields);
}

BOOST_AUTO_TEST_SUITE_END()
