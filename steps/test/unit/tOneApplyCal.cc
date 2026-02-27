// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/OneApplyCal.h"

#include <array>
#include <filesystem>

#include <boost/test/unit_test.hpp>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/soltab.h>

#include "common/ParameterSet.h"

#include "mock/MockStep.h"

using dp3::steps::OneApplyCal;
using dp3::steps::Step;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace {
const std::string kTempH5ParmName = "tOneApplyCal_tmp.h5";
const std::string kAntennaName = "dummy antenna";
constexpr double kFirstTime = 42.0;
constexpr double kLastTime = 43.0;
constexpr double kInterval = 1.0;
constexpr double kFrequency = 100.0;
constexpr double kChannelWidth = 1.0;

struct TestOneApplyCalFixture {
 public:
  TestOneApplyCalFixture() {
    CreateTempH5Parm();

    parset.add("parmdb", kTempH5ParmName);
    parset.add("correction", "amplitude000");
  }

  ~TestOneApplyCalFixture() {
    std::filesystem::remove(kTempH5ParmName.c_str());
  }

  void CreateTempH5Parm() {
    H5Parm h5parm(kTempH5ParmName, true);

    const std::vector<std::string> names{kAntennaName};
    const std::vector<std::array<double, 3>> positions{{0.0, 0.0, 0.0}};
    h5parm.AddAntennas(names, positions);

    const std::vector<double> times{kFirstTime, kLastTime};
    const std::vector<double> frequencies{kFrequency};
    const std::vector<double> values(
        times.size() * frequencies.size() * names.size(), 1.0);
    const std::vector<double> weights(values.size(), 1.0);

    std::vector<schaapcommon::h5parm::AxisInfo> axes;
    axes.push_back({"ant", static_cast<unsigned>(names.size())});
    axes.push_back({"time", static_cast<unsigned>(times.size())});
    axes.push_back({"freq", static_cast<unsigned>(frequencies.size())});

    SolTab soltab = h5parm.CreateSolTab("amplitude000", "amplitude", axes);
    soltab.SetTimes(times);
    soltab.SetFreqs(frequencies);
    soltab.SetAntennas(names);
    soltab.SetValues(values, weights, "temporary OneApplyCal test table");
  }

  dp3::common::ParameterSet parset;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(one_apply_cal)

BOOST_FIXTURE_TEST_CASE(fields, TestOneApplyCalFixture) {
  const OneApplyCal one_apply_cal(parset, "", "");
  BOOST_TEST(one_apply_cal.getRequiredFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
  BOOST_TEST(one_apply_cal.getProvidedFields() ==
             (Step::kDataField | Step::kFlagsField));

  parset.add("updateweights", "true");
  const OneApplyCal updates_weights(parset, "", "");
  BOOST_TEST(updates_weights.getRequiredFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
  BOOST_TEST(updates_weights.getProvidedFields() ==
             (Step::kDataField | Step::kWeightsField | Step::kFlagsField));
}

// Test a corner case where the time of the buffer is beyond the last time.
BOOST_FIXTURE_TEST_CASE(process_beyond_end, TestOneApplyCalFixture) {
  // Make ApplyCal update parameters every time slot, since that triggered
  // the original bug (GEC-279).
  parset.add("timeslotsperparmupdate", "1");

  OneApplyCal one_apply_cal(parset, "", "");
  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  one_apply_cal.setNextStep(mock_step);

  // The bug from GEC-279 occured when the last time is a multiple of the
  // time interval + a value between between 0.0 and 0.5 * kInterval.
  constexpr double kBufferTime = kLastTime + 0.1 * kInterval;
  constexpr size_t kNCorrelations = 4;

  dp3::base::DPInfo info(kNCorrelations);
  // Also use the updated last time in DPInfo.
  info.setTimes(kFirstTime, kBufferTime, kInterval);
  info.setChannels(std::vector<double>(1, kFrequency),
                   std::vector<double>(1, kChannelWidth));
  info.setAntennas(std::vector<std::string>(1, kAntennaName),
                   std::vector<double>(1, 1.0),
                   std::vector<casacore::MPosition>(1), std::vector<int>(1, 0),
                   std::vector<int>(1, 0));
  one_apply_cal.updateInfo(info);

  // Pass an empty buffer the adjusted last time.
  auto buffer = std::make_unique<dp3::base::DPBuffer>();
  buffer->SetTime(kBufferTime);
  BOOST_CHECK_NO_THROW(one_apply_cal.process(std::move(buffer)));

  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers().size(), 1);
  BOOST_CHECK_EQUAL(mock_step->GetRegularBuffers()[0]->GetTime(), kBufferTime);
}

BOOST_AUTO_TEST_SUITE_END()
