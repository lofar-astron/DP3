// tUVWFlagger.cc: Test program for class UVWFlagger
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "../../UVWFlagger.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "mock/ThrowStep.h"
#include <dp3/base/BDABuffer.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Step;
using dp3::steps::UVWFlagger;
using std::vector;

namespace {

// Simple class to generate input arrays (Regular or BDA).
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
template <class T>
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(size_t n_times, size_t n_baselines, size_t n_channels,
            size_t n_correlations)
      : count_(0),
        n_times_(n_times),
        n_baselines_(n_baselines),
        n_channels_(n_channels),
        n_correlations_(n_correlations) {
    info() = DPInfo(n_correlations, n_channels);
    info().setTimes(0.0, (n_times - 1) * time_interval_, time_interval_);
    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    vector<int> ant1(n_baselines);
    vector<int> ant2(n_baselines);
    int st1 = 0;
    int st2 = 0;
    for (size_t i = 0; i < n_baselines; ++i) {
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
    vector<double> antDiam(4, 70.);
    info().setAntennas(antNames, antDiam, antPos, ant1, ant2);
  }

 private:
  bool process(const DPBuffer&) override {
    // Stop when all times are done.
    if (count_ == n_times_) {
      return false;
    }
    T buffer = CreateInputBuffer();
    getNextStep()->process(std::move(buffer));
    ++count_;
    return true;
  };

  T CreateInputBuffer() {
    T t;
    return t;
  };
  void finish() override { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  size_t count_, n_times_, n_baselines_, n_channels_, n_correlations_;
  double start_frequency_ = 10000000.0;
  double min_channel_width_ = 1000000.0;
  double time_interval_ = 5.0;
  double first_time_ = 4472025740.0;
};

template <>
DPBuffer TestInput<DPBuffer>::CreateInputBuffer() {
  casacore::Cube<casacore::Complex> data(n_correlations_, n_channels_,
                                         n_baselines_);
  for (size_t i = 0; i < data.size(); ++i) {
    data.data()[i] = casacore::Complex(i + count_ * 10, i - 10 + count_ * 6);
  }
  casacore::Matrix<double> uvw(3, n_baselines_);
  for (size_t i = 0; i < n_baselines_; ++i) {
    uvw(0, i) = 1 + count_ + i;
    uvw(1, i) = 2 + count_ + i;
    uvw(2, i) = 3 + count_ + i;
  }
  DPBuffer buf;
  buf.setTime(count_ * 30 + first_time_);
  buf.setData(data);
  buf.setUVW(uvw);
  casacore::Cube<bool> flags(data.shape());
  flags = false;
  buf.setFlags(flags);

  return buf;
}

template <>
std::unique_ptr<BDABuffer>
TestInput<std::unique_ptr<BDABuffer>>::CreateInputBuffer() {
  std::unique_ptr<BDABuffer> buffer = std::make_unique<BDABuffer>(
      n_correlations_ * n_channels_ * n_baselines_ * n_times_);

  std::vector<std::vector<double>> channel_frequencies = info().BdaChanFreqs();

  const double bda_first_time = count_ * 30 + first_time_;

  for (size_t i = 0; i < n_baselines_; ++i) {
    size_t n_channels = channel_frequencies[i].size();
    std::vector<std::complex<float>> data(n_correlations_ * n_channels);
    for (size_t i = 0; i < data.size(); ++i) {
      data.data()[i] = casacore::Complex(i + count_ * 10, i - 10 + count_ * 6);
    }

    const double uvw[3]{1.0 + count_ + i, 2.0 + count_ + i, 3.0 + count_ + i};
    const std::vector<float> weights(n_correlations_ * n_channels, 1.0);
    // AddRow requires a contiguous set of bool elements. (std::vector<bool>
    // does not satisfy this requirement.)
    auto flags = std::make_unique<bool[]>(n_correlations_ * n_channels);
    buffer->AddRow(bda_first_time, time_interval_, time_interval_, i,
                   n_channels, n_correlations_, data.data(), flags.get(),
                   weights.data(), nullptr, uvw);
  }
  return buffer;
}

template <>
void TestInput<DPBuffer>::updateInfo(const DPInfo&) {
  // Define the frequencies for the regular buffer.
  std::vector<double> channel_widths(n_channels_, min_channel_width_);
  std::vector<double> channel_frequencies;
  for (size_t i = 0; i < n_channels_; i++) {
    channel_frequencies.push_back(start_frequency_ +
                                  (i + 0.5) * min_channel_width_);
  }
  info().setChannels(std::move(channel_frequencies), std::move(channel_widths));
}

template <>
void TestInput<std::unique_ptr<BDABuffer>>::updateInfo(const DPInfo&) {
  // Define the frequencies for the BDA buffer.
  // Even baselines have a channel averaging factor of 5
  // Odd baselines have a channel averaging factor of 3
  // For a precise mapping between regular and BDA results, the original number
  // of channel should be divisible by 3 and 5.
  std::vector<std::vector<double>> channel_frequencies(n_baselines_);
  std::vector<std::vector<double>> channel_widths(n_baselines_);

  for (size_t k = 0; k < n_baselines_; k++) {
    if (k % 2 == 0) {
      channel_widths[k] = std::vector<double>(static_cast<int>(n_channels_ / 5),
                                              5 * min_channel_width_);
    } else {
      channel_widths[k] = std::vector<double>(static_cast<int>(n_channels_ / 3),
                                              3 * min_channel_width_);
    }
    for (size_t i = 0; i < channel_widths[k].size(); ++i) {
      channel_frequencies[k].push_back(start_frequency_ +
                                       (i + 0.5) * channel_widths[k][i]);
    }
  }

  info().setChannels(std::move(channel_frequencies), std::move(channel_widths));
}

// Class to check result of flagged, unaveraged TestInput run by test1.
template <class T>
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(size_t n_times, size_t n_baselines, size_t n_channels,
             size_t n_correlations, size_t test_id)
      : count_(0),
        n_times_(n_times),
        n_baselines_(n_baselines),
        n_channels_(n_channels),
        n_correlations_(n_correlations),
        test_id_(test_id) {}

 private:
  bool process(const DPBuffer&) override { return false; }
  bool process(std::unique_ptr<BDABuffer>) override { return false; }
  const casacore::Cube<bool> GetResult() const {
    casacore::Cube<bool> result;
    switch (test_id_) {
      case 1:
        result = GetResultTest1();
        break;
      case 2:
        result = GetResultTest2();
        break;
      case 3:
        result = GetResultTest3();
        break;
    }
    return result;
  }

  const casacore::Cube<bool> GetResultTest1() const {
    casacore::Cube<bool> result(n_correlations_, n_channels_, n_baselines_);
    result = false;
    for (size_t i = 0; i < n_baselines_; ++i) {
      double u = 1 + i + count_;
      double v = 2 + i + count_;
      double w = 3 + i + count_;
      double uv = sqrt(u * u + v * v);
      if ((uv > 5.5 && uv < 8.5) || (u > 20.5 && u < 23.5) ||
          (u > 31.5 && u < 40.5) || (v > 11.5 && v < 14.5) || w < 3.5 ||
          w > 44.5) {
        for (size_t j = 0; j < n_channels_; ++j) {
          for (size_t k = 0; k < n_correlations_; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    return result;
  }

  const casacore::Cube<bool> GetResultTest2() const {
    casacore::Cube<bool> result(n_correlations_, n_channels_, n_baselines_);
    result = false;
    for (size_t i = 0; i < n_baselines_; ++i) {
      for (size_t j = 0; j < n_channels_; ++j) {
        double wavel = 2.99792458e+08 / (10.5e6 + j * 1e6);
        double u = (1 + i + count_) / wavel;
        double v = (2 + i + count_) / wavel;
        double w = (3 + i + count_) / wavel;
        double uv = sqrt(u * u + v * v);
        if ((uv > 0.2 && uv < 0.31) || (u > 1.55 && u < 1.485) ||
            (u > 0.752 && u < 0.862) || (v > 0.42 && v < 0.53) || w < 0.12 ||
            w > 1.63) {
          for (size_t k = 0; k < n_correlations_; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    return result;
  }

  const casacore::Cube<bool> GetResultTest3() const {
    // These are the UVW coordinates as calculated by UVWFlagger for the
    // station positions and times defined in TestInput and phase center
    // defined in test3.
    double uvwvals[] = {
        0,         0,        0,        0.423756,  -127.372, 67.1947,
        0.847513,  -254.744, 134.389,  0.277918,  -382.015, 201.531,
        -0.423756, 127.372,  -67.1947, 0,         0,        0,
        0.423756,  -127.372, 67.1947,  -0.145838, -254.642, 134.336,
        -0.847513, 254.744,  -134.389, -0.423756, 127.372,  -67.1947,
        0,         0,        0,        -0.569594, -127.27,  67.1417,
        -0.277918, 382.015,  -201.531, 0.145838,  254.642,  -134.336,
        0.569594,  127.27,   -67.1417, 0,         0,        0,
        0,         0,        0,        0.738788,  -127.371, 67.1942,
        1.47758,   -254.742, 134.388,  1.22276,   -382.013, 201.53,
        -0.738788, 127.371,  -67.1942, 0,         0,        0,
        0.738788,  -127.371, 67.1942,  0.483976,  -254.642, 134.336,
        -1.47758,  254.742,  -134.388, -0.738788, 127.371,  -67.1942,
        0,         0,        0,        -0.254812, -127.271, 67.1421,
        -1.22276,  382.013,  -201.53,  -0.483976, 254.642,  -134.336,
        0.254812,  127.271,  -67.1421, 0,         0,        0};
    casacore::Cube<double> uvws(casacore::IPosition(3, 3, 16, 2), uvwvals,
                                casacore::SHARE);
    // Flag where u,v,w matches intervals given in test3.
    casacore::Cube<bool> result(n_correlations_, n_channels_, n_baselines_);
    result = false;
    for (size_t i = 0; i < n_baselines_; ++i) {
      double u = uvws(0, i, count_);
      double v = uvws(1, i, count_);
      double w = uvws(2, i, count_);
      double uv = sqrt(u * u + v * v);
      if ((uv > 5.5 && uv < 8.5) || (u > 20.5 && u < 23.5) ||
          (u > 31.5 && u < 40.5) || (v > 11.5 && v < 14.5) || w < 3.5 ||
          w > 44.5) {
        for (size_t j = 0; j < n_channels_; ++j) {
          for (size_t k = 0; k < n_correlations_; ++k) {
            result(k, j, i) = true;
          }
        }
      }
    }
    return result;
  }

  void finish() override {}
  void updateInfo(const DPInfo& infoIn) override {
    info() = infoIn;
    BOOST_CHECK_EQUAL(infoIn.origNChan(), n_channels_);
    BOOST_CHECK_EQUAL(infoIn.ntime(), n_times_);
    BOOST_CHECK_EQUAL(infoIn.timeInterval(), 5.0);
    BOOST_CHECK_EQUAL(infoIn.ntimeAvg(), 1u);
  }

  size_t count_;
  size_t n_times_, n_baselines_, n_channels_, n_correlations_, test_id_;
};

template <>
bool TestOutput<DPBuffer>::process(const DPBuffer& buf) {
  // Flag where u,v,w matches intervals given in the requested test.
  const casacore::Cube<bool> expected_result = GetResult();

  BOOST_CHECK(allEQ(buf.getFlags(), expected_result));
  count_++;
  return true;
}

template <>
bool TestOutput<std::unique_ptr<BDABuffer>>::process(
    const std::unique_ptr<BDABuffer> buf) {
  // Flag where u,v,w matches intervals given in the requested test.
  std::vector<std::vector<double>> channel_frequencies = info().BdaChanFreqs();
  const casacore::Cube<bool> expected_result = GetResult();

  for (auto row : buf->GetRows()) {
    auto flag_ptr = row.flags;
    // Use the scaling factor to get the index corresponding to the current
    // frequency channel in the non-averaged result cube
    double scaling_factor =
        double(n_channels_) / channel_frequencies[row.baseline_nr].size();
    for (size_t j = 0; j < channel_frequencies[row.baseline_nr].size(); ++j) {
      double frequency_index = (j + 0.5) * scaling_factor;
      for (size_t k = 0; k < n_correlations_; ++k) {
        BOOST_CHECK(*flag_ptr ==
                    expected_result(k, static_cast<int>(frequency_index),
                                    row.baseline_nr));
        ++flag_ptr;
      }
    }
  }
  ++count_;
  return true;
}

// Test flagging a few baselines on UV in m.
template <class T>
void test1(size_t n_times, size_t n_baselines, size_t n_channels,
           size_t n_correlations, Step::MsType input_type) {
  // Create the steps.
  auto in = std::make_shared<TestInput<T>>(n_times, n_baselines, n_channels,
                                           n_correlations);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add("vmrange", "[11.5..14.5]");
  parset.add("wmmax", "44.5");
  parset.add("wmmin", "3.5");
  auto flagger = std::make_shared<UVWFlagger>(parset, "", input_type);
  auto out = std::make_shared<TestOutput<T>>(n_times, n_baselines, n_channels,
                                             n_correlations, 1);
  dp3::steps::test::Execute({in, flagger, out});
}

// Test flagging a few baselines on UV in wavelengths.
template <class T>
void test2(size_t n_times, size_t n_baselines, size_t n_channels,
           size_t n_correlations, Step::MsType input_type) {
  // Create the steps.
  auto in = std::make_shared<TestInput<T>>(n_times, n_baselines, n_channels,
                                           n_correlations);
  dp3::common::ParameterSet parset;
  parset.add("uvlambdarange", "[0.2..0.31]");
  parset.add("ulambdarange", "[1.55..1.485, 0.807+-0.055]");
  parset.add("vlambdarange", "[0.42..0.53]");
  parset.add("wlambdamax", "1.63");
  parset.add("wlambdamin", "0.12");
  auto flagger = std::make_shared<UVWFlagger>(parset, "", input_type);
  auto out = std::make_shared<TestOutput<T>>(n_times, n_baselines, n_channels,
                                             n_correlations, 2);
  dp3::steps::test::Execute({in, flagger, out});
}

// Test flagging a few baselines on UV in m with a different phase center.
template <class T>
void test3(size_t n_times, size_t n_baselines, size_t n_channels,
           size_t n_correlations, Step::MsType input_type) {
  // Create the steps.
  auto in = std::make_shared<TestInput<T>>(n_times, n_baselines, n_channels,
                                           n_correlations);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("umrange", "[31.5..40.5, 22+-1.5]");
  parset.add("vmrange", "[11.5..14.5]");
  parset.add("wmmax", "44.5");
  parset.add("wmmin", "3.5");
  parset.add("phasecenter", "[-1.92653768rad, 1.09220917rad, j2000]");
  auto flagger = std::make_shared<UVWFlagger>(parset, "", input_type);
  BOOST_REQUIRE_EQUAL(flagger->isDegenerate(), false);

  auto out = std::make_shared<TestOutput<T>>(n_times, n_baselines, n_channels,
                                             n_correlations, 3);
  dp3::steps::test::Execute({in, flagger, out});
}

// Test flagging a few baselines on UV in m with a different phase center.
template <class T>
void test_constructor(size_t n_times, size_t n_baselines, size_t n_channels,
                      size_t n_correlations, Step::MsType input_type) {
  // Create the steps.
  TestInput<T> in(n_times, n_baselines, n_channels, n_correlations);
  dp3::common::ParameterSet parset;
  UVWFlagger uvw_flagger_step(parset, "", input_type);
  BOOST_REQUIRE_EQUAL(uvw_flagger_step.isDegenerate(), true);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(uvwflagger)

BOOST_AUTO_TEST_CASE(constructor) {
  test_constructor<DPBuffer>(10, 16, 32, 4, Step::MsType::kRegular);
  test_constructor<std::unique_ptr<BDABuffer>>(10, 16, 32, 4,
                                               Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(fields) {
  using dp3::steps::Step;

  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  const UVWFlagger degenerate(parset, "", Step::MsType::kRegular);
  BOOST_TEST(degenerate.isDegenerate());
  BOOST_TEST(degenerate.getRequiredFields() == dp3::common::Fields());
  BOOST_TEST(degenerate.getProvidedFields() == dp3::common::Fields());

  parset.add("ummax", "42");
  const UVWFlagger without_center(parset, "", Step::MsType::kBda);
  BOOST_TEST(without_center.getRequiredFields() ==
             (Step::kFlagsField | Step::kUvwField));
  BOOST_TEST(without_center.getProvidedFields() == Step::kFlagsField);

  parset.add("phasecenter", "Jupiter");
  const UVWFlagger with_center(parset, "", Step::MsType::kRegular);
  BOOST_TEST(with_center.getRequiredFields() == Step::kFlagsField);
  BOOST_TEST(with_center.getProvidedFields() == Step::kFlagsField);
}

BOOST_AUTO_TEST_CASE(test1_regular) {
  test1<DPBuffer>(10, 16, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test2_regular) {
  test1<DPBuffer>(100, 105, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test3_regular) {
  test2<DPBuffer>(2, 16, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test4_regular) {
  test2<DPBuffer>(2, 36, 16, 2, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test5_regular) {
  test2<DPBuffer>(10, 16, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test6_regular) {
  test2<DPBuffer>(100, 105, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test7_regular) {
  test3<DPBuffer>(2, 16, 32, 4, Step::MsType::kRegular);
}

BOOST_AUTO_TEST_CASE(test1_bda) {
  test1<std::unique_ptr<BDABuffer>>(10, 16, 32, 4, Step::MsType::kBda);
}
BOOST_AUTO_TEST_CASE(test2_bda) {
  test1<std::unique_ptr<BDABuffer>>(100, 105, 32, 4, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(test3_bda) {
  test2<std::unique_ptr<BDABuffer>>(2, 16, 15, 4, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(test4_bda) {
  test2<std::unique_ptr<BDABuffer>>(2, 36, 15, 2, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(test5_bda) {
  test2<std::unique_ptr<BDABuffer>>(10, 16, 30, 4, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(test6_bda) {
  test2<std::unique_ptr<BDABuffer>>(100, 105, 30, 4, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(test7_bda) {
  test3<std::unique_ptr<BDABuffer>>(2, 16, 32, 4, Step::MsType::kBda);
}

BOOST_AUTO_TEST_CASE(sun_as_phase_center) {
  auto in = std::make_shared<TestInput<DPBuffer>>(1, 1, 1, 1);
  dp3::common::ParameterSet parset;
  parset.add("uvmrange", "[5.5..8.5]");
  parset.add("phasecenter", "Sun");
  BOOST_CHECK_NO_THROW(
      std::make_shared<UVWFlagger>(parset, "", Step::MsType::kRegular));
}

BOOST_AUTO_TEST_SUITE_END()
