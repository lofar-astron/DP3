// tDemixer.cc: Test program for class Demixer
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../Demixer.h"
#include "../../DemixerNew.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

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
using dp3::steps::DemixerNew;
using dp3::steps::Step;

namespace {

// Demixer only works with 4 correlations
const int kNCorr = 4;

// Simple class to generate input arrays.
// It can only set all flags to true or all false.
// Weights are always 1.
// It can be used with different nr of times, channels, baselines.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(int ntime, int nbl, int nchan, bool flag)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        flag_data_(flag) {}

 private:
  bool process(const DPBuffer&) override {
    // Stop when all times are done.
    if (count_ == n_times_) {
      return false;
    }
    casacore::Cube<casacore::Complex> data(kNCorr, n_channels_, n_baselines_);
    for (int i = 0; i < int{data.size()}; ++i) {
      data.data()[i] =
          casacore::Complex(i + count_ * 10, i - 1000 + count_ * 6);
    }
    DPBuffer buf;
    buf.setTime(count_ * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    casacore::Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    casacore::Cube<bool> flags(data.shape(), flag_data_);
    buf.setFlags(flags);
    casacore::Matrix<double> uvw(3, n_baselines_);
    indgen(uvw, double(count_ * 100));
    buf.setUVW(uvw);
    getNextStep()->process(buf);
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
    for (int i = 0; i < n_channels_; i++) {
      chan_freqs.push_back(start_channel + i * 100000.);
    }
    info().setChannels(std::move(chan_freqs), std::move(chan_width));
    info().update(start_channel, n_channels_, baselines, false);
  }

  int count_;
  int n_times_;
  int n_baselines_;
  int n_channels_;
  bool flag_data_;
};

// Class to check result of averaging TestInput.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, int nbl, int nchan, int navgtime, int navgchan,
             bool flag)
      : count_(0),
        n_times_(ntime),
        n_baselines_(nbl),
        n_channels_(nchan),
        n_average_time_(navgtime),
        n_average_channel_(navgchan),
        flag_data_(flag) {}

 private:
  bool process(const DPBuffer& buf) override {
    int nchan = 1 + (n_channels_ - 1) / n_average_channel_;
    int navgtime =
        std::min(n_average_time_, n_times_ - count_ * n_average_time_);
    // Fill expected result in similar way as TestInput.
    casacore::Cube<casacore::Complex> data(kNCorr, n_channels_, n_baselines_);
    casacore::Cube<float> weights(kNCorr, n_channels_, n_baselines_);
    weights = 0;
    if (!flag_data_) {
      for (int j = count_ * n_average_time_;
           j < count_ * n_average_time_ + navgtime; ++j) {
        for (int i = 0; i < int{data.size()}; ++i) {
          data.data()[i] += casacore::Complex(i + j * 10, i - 1000 + j * 6);
          weights.data()[i] += float(1);
        }
      }
    }
    casacore::Cube<casacore::Complex> result(kNCorr, nchan, n_baselines_);
    casacore::Cube<float> resultw(kNCorr, nchan, n_baselines_);
    resultw = 0;
    // Average to get the true expected result.
    for (int k = 0; k < n_baselines_; ++k) {
      for (int i = 0; i < kNCorr; ++i) {
        for (int j = 0; j < nchan; ++j) {
          int jc;
          for (jc = j * n_average_channel_;
               jc < std::min((j + 1) * n_average_channel_, n_channels_); ++jc) {
            result(i, j, k) += data(i, jc, k);
            resultw(i, j, k) += weights(i, jc, k);
          }
          result(i, j, k) /= float(navgtime * (jc - j * n_average_channel_));
        }
      }
    }
    // Check the averaged result. When all the flags are set, do not check the
    // values of the data and weights, since DP3 should not use those values
    // anyway.
    if (!flag_data_) {
      BOOST_CHECK(allNear(real(buf.GetCasacoreData()), real(result), 1.0e-5));
      BOOST_CHECK(allNear(imag(buf.GetCasacoreData()), imag(result), 1.0e-5));
      BOOST_CHECK(allNear(buf.GetCasacoreWeights(), resultw, 1.0e-5));
    }
    if (navgtime == n_average_time_) {
      casacore::Matrix<double> uvw(3, n_baselines_);
      indgen(uvw,
             100 * (count_ * n_average_time_ + 0.5 * (n_average_time_ - 1)));
      BOOST_CHECK(allNear(buf.GetCasacoreUvw(), uvw, 1e-5));
    }
    BOOST_CHECK(allEQ(buf.GetCasacoreFlags(), flag_data_));
    BOOST_CHECK_CLOSE(
        buf.getTime(),
        2.0 + 5.0 * (count_ * n_average_time_ + (n_average_time_ - 1) / 2.0),
        1.0e-3);

    ++count_;
    return true;
  }

  void finish() override {}
  virtual void updateInfo(const DPInfo& info) override {
    BOOST_CHECK_EQUAL(size_t{8}, info.startchan());
    BOOST_CHECK_EQUAL(n_channels_, int{info.origNChan()});
    BOOST_CHECK_EQUAL(1 + (n_channels_ - 1) / n_average_channel_,
                      int{info.nchan()});
    BOOST_CHECK_EQUAL(1 + (n_times_ - 1) / n_average_time_, int{info.ntime()});
    BOOST_CHECK_EQUAL(5 * n_average_time_, info.timeInterval());
    BOOST_CHECK_EQUAL(n_average_channel_, int{info.nchanAvg()});
    BOOST_CHECK_EQUAL(n_average_time_, int{info.ntimeAvg()});
  }

  int count_;
  int n_times_;
  int n_baselines_;
  int n_channels_;
  int n_average_time_;
  int n_average_channel_;
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

  ParameterSet parset_new;
  parset_new.add("ateam.skymodel", dp3::steps::test::kPredictSourceDB);
  parset_new.add("target.skymodel", dp3::steps::test::kPredictSourceDB);
  parset_new.add("sources", dp3::steps::test::kPredictDirection);
  DemixerNew demixer_new(parset_new, "");

  BOOST_TEST(demixer.getRequiredFields() == Averager::kRequiredFields);
  BOOST_TEST(demixer.getProvidedFields() == Averager::kProvidedFields);
  BOOST_TEST(demixer_new.getRequiredFields() == Averager::kRequiredFields);
  BOOST_TEST(demixer_new.getProvidedFields() == Averager::kProvidedFields);
}

BOOST_AUTO_TEST_SUITE_END()
