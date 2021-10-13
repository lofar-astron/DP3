// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../BdaDdeCal.h"
#include "../../InputStep.h"

#include <aocommon/threadpool.h>

BOOST_AUTO_TEST_SUITE(bdaddecal)

/// Helper forwarder to improve test failure output.
static void TestShow(
    const std::string& expected,
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  const dp3::common::ParameterSet parset =
      dp3::steps::test::CreateParameterSet(parameters);
  BOOST_CHECK_EQUAL(expected,
                    dp3::steps::test::Show<dp3::steps::BdaDdeCal>(parset));

  BOOST_CHECK_MESSAGE(parset.unusedKeys().empty(),
                      "Not all keys are used, is there a typo?");
}

BOOST_AUTO_TEST_CASE(show_default) {
  TestShow(
      R"(BdaDdeCal prefix.
  mode (constraints):  diagonal
  directions:          [[center]]
  solver algorithm:    directionsolve
  H5Parm:              tDDECal.MS/instrument.h5
  subtract model:      false
  solution interval:   0 s
  #channels/block:     1
  #channel blocks:     18446744073709551615
  tolerance:           0.0001
  max iter:            50
  flag unconverged:    false
     diverged only:    false
  propagate solutions: false
       converged only: false
  detect stalling:     true
  step size:           0.2
Model steps for direction [center]
Predict
BDAExpander prefix.
OnePredict prefix.
  sourcedb:           tDDECal.MS/sky
   number of patches: 1
   number of sources: 1
   all unpolarized:   true
   correct freq smearing: false
  apply beam:         false
  operation:          replace
  threads:            )" +
          std::to_string(aocommon::ThreadPool::NCPUs()) + R"(
BDAAverager prefix.
  timebase:        0s
  max interval:    0s
  frequencybase:   0
  min channels:    1
  max time factor: 1
  max freq factor: 1

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky"}});
}

BOOST_AUTO_TEST_CASE(show_modified) {
  TestShow(
      R"(BdaDdeCal prefix.
  mode (constraints):  diagonal
  directions:          [[center]]
  solver algorithm:    hybrid
  H5Parm:              tDDECal.MS/instrument.h5
  subtract model:      true
  solution interval:   0 s
  #channels/block:     44
  #channel blocks:     18446744073709551615
  tolerance:           1e-05
  max iter:            49
  flag unconverged:    true
     diverged only:    true
  propagate solutions: true
       converged only: true
  detect stalling:     true
  step size:           0.2
  coreconstraint:      45.123
  smoothnessconstraint:46.123
  smoothnessreffrequency:47.123
  smoothnessrefdistance:48.123
  tecscreen.coreconstraint:49.123
Model steps for direction [center]
Predict
BDAExpander prefix.
OnePredict prefix.
  sourcedb:           tDDECal.MS/sky
   number of patches: 1
   number of sources: 1
   all unpolarized:   true
   correct freq smearing: false
  apply beam:         false
  operation:          replace
  threads:            )" +
          std::to_string(aocommon::ThreadPool::NCPUs()) + R"(
BDAAverager prefix.
  timebase:        0s
  max interval:    0s
  frequencybase:   0
  min channels:    1
  max time factor: 1
  max freq factor: 1

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.propagatesolutions", "true"},
       {"prefix.propagateconvergedonly", "true"},
       {"prefix.flagunconverged", "true"},
       {"prefix.flagdivergedonly", "true"},
       {"prefix.subtract", "true"},
       {"prefix.solveralgorithm", "hybrid"},
       {"prefix.sourcedb", "tDDECal.MS/sky"},
       {"prefix.nchan", "44"},
       {"prefix.coreconstraint", "45.123"},
       {"prefix.smoothnessconstraint", "46.123"},
       {"prefix.smoothnessreffrequency", "47.123"},
       {"prefix.smoothnessrefdistance", "48.123"},
       {"prefix.tecscreen.coreconstraint", "49.123"},
       {"prefix.maxiter", "49"}});
}

BOOST_AUTO_TEST_SUITE_END()
