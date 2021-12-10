// tScaleData.cc: Test program for class ScaleData
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Lars Krombeen

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../../ScaleData.h"
#include "../../InputStep.h"
#include "../../../base/DPBuffer.h"
#include "../../../base/BDABuffer.h"
#include "../../../base/DPInfo.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"
#include "../../../common/StreamUtil.h"

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::ScaleData;
using std::vector;

namespace {
const float kFreq = 10.5;  // MHz
}

BOOST_AUTO_TEST_SUITE(scaledata_bda)

// Class to check result of TestInput run by test1.
class TestOutput : public dp3::steps::Step {
 public:
  TestOutput(int ntime, int nbl, int nchan, int ncorr)
      : count_(0), ntime_(ntime), nbl_(nbl), nchan_(nchan), ncorr_(ncorr) {}

 public:
  std::unique_ptr<BDABuffer> results_;

 private:
  virtual bool process(const DPBuffer&) { return true; }
  virtual bool process(std::unique_ptr<BDABuffer> results) {
    results_ = std::move(results);
    return true;
  }
  virtual void finish() {}
  virtual void show(std::ostream&) const {}
  virtual void updateInfo(const DPInfo&) {}

 private:
  int count_;
  int ntime_, nbl_, nchan_, ncorr_;
};

// Generate DP Info for 2 antennas
DPInfo GenerateDPInfo(int ntime, int nbl, int nchan, int ncorr) {
  DPInfo info = DPInfo();
  info.init(ncorr, 0, nchan, ntime, 0., 5., string(), string());
  // Fill the baseline stations; use 2 stations.
  // So they are called 00 01 10 11 00 01 ...
  vector<int> ant1(nbl);
  vector<int> ant2(nbl);
  for (int i = 0; i < nbl; ++i) {
    ant1[i] = i / 2;
    ant2[i] = i % 2;
  }
  vector<string> ant_names(2);
  ant_names[0] = "rs01.s01";
  ant_names[1] = "rs02.s01";
  // Define their positions (more or less WSRT RT0-3).
  vector<casacore::MPosition> ant_pos(2);
  vector<double> vals(3);
  vals[0] = 3828763;
  vals[1] = 442449;
  vals[2] = 5064923;
  ant_pos[0] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vals[0] = 3828746;
  vals[1] = 442592;
  vals[2] = 5064924;
  ant_pos[1] = casacore::MPosition(
      casacore::Quantum<casacore::Vector<double>>(vals, "m"),
      casacore::MPosition::ITRF);
  vector<double> antDiam(2, 70.);
  info.set(ant_names, antDiam, ant_pos, ant1, ant2);
  // Define the frequencies.
  std::vector<double> chan_width(nchan, 1e6);
  std::vector<double> chan_freqs(nchan, kFreq * 1e6);
  info.set(std::move(chan_freqs), std::move(chan_width));

  return info;
}

BOOST_AUTO_TEST_CASE(test_processing_for_bda_buffer) {
  int ntime{10};
  int nbl{4};
  int nchan{1};
  int ncorr{2};
  const int datasize{nbl * nchan * ncorr};
  const std::vector<float> coeffs{2 + 1 * kFreq,
                                  3 + 2 * kFreq + 1 * kFreq * kFreq};

  // Preparation
  ParameterSet parset;
  parset.add("stations", "[rs01.s01, *]");
  parset.add("coeffs", "[[2,1],[3,2,1]]");
  parset.add("scalesize", "false");

  DPInfo info = GenerateDPInfo(ntime, nbl, nchan, ncorr);

  auto step_scale_data =
      std::make_shared<ScaleData>(nullptr, parset, "", ScaleData::MsType::kBda);
  auto step_test_output =
      std::make_shared<TestOutput>(ntime, nbl, nchan, ncorr);
  step_scale_data->setNextStep(step_test_output);
  step_scale_data->setInfo(info);

  // Initialize buffer
  std::unique_ptr<BDABuffer> bda_buffer{new BDABuffer(datasize)};
  std::vector<std::complex<float>> data;
  std::vector<std::complex<float>> expected_data;
  for (int bl = 0; bl < nbl; ++bl) {
    const float scale =
        std::sqrt(coeffs[info.getAnt1()[bl]] * coeffs[info.getAnt2()[bl]]);
    for (int i = 0; i < nchan * ncorr; ++i) {
      data.emplace_back(10.0 * bl + i + 1, -10.0 * bl - i - 1);
      expected_data.push_back(data.back() * scale);
    }
    bda_buffer->AddRow(ntime, 5., 5., bl, nchan, ncorr,
                       data.data() + bl * nchan * ncorr);
  }

  // Execution
  step_scale_data->process(std::move(bda_buffer));

  // Assertion
  BOOST_REQUIRE_EQUAL(expected_data.size(),
                      step_test_output->results_->GetNumberOfElements());
  const auto* results = step_test_output->results_->GetData();
  for (std::size_t i = 0; i < expected_data.size(); ++i) {
    BOOST_TEST(expected_data[i] == results[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
