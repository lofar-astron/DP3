// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../DDECal.h"

using dp3::steps::DDECal;

BOOST_AUTO_TEST_SUITE(ddecal)

BOOST_AUTO_TEST_CASE(provided_fields) {
  using dp3::steps::Step;

  dp3::steps::MockInput input;
  dp3::common::ParameterSet parset;
  parset.add("msin", "tDDECal.MS");
  parset.add("prefix.directions", "[[center]]");
  parset.add("prefix.sourcedb", "tDDECal.MS/sky.txt");
  // Errors occur when creating multiple DDECal's that write the same H5 file.
  // Putting each DDECal object in a different scope also works.
  parset.add("prefix.h5parm", "tddecal_fields.h5");
  const DDECal ddecal(&input, parset, "prefix.");
  BOOST_TEST(ddecal.getProvidedFields() == dp3::common::Fields());

  dp3::common::ParameterSet parset_only_predict = parset;
  parset_only_predict.add("prefix.onlypredict", "true");
  parset_only_predict.replace("prefix.h5parm",
                              "tddecal_fields_only_predict.h5");
  const DDECal ddecal_only_predict(&input, parset_only_predict, "prefix.");
  BOOST_TEST(ddecal_only_predict.getProvidedFields() == Step::kDataField);

  dp3::common::ParameterSet parset_subtract = parset;
  parset_subtract.add("prefix.subtract", "true");
  parset_subtract.replace("prefix.h5parm", "tddecal_fields_subtract.h5");
  const DDECal ddecal_subtract(&input, parset_subtract, "prefix.");
  BOOST_TEST(ddecal_subtract.getProvidedFields() == Step::kDataField);
}

/// Helper forwarder to improve test failure output.
static void TestShow(
    const std::string& expected,
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  dp3::steps::MockInput input;
  const dp3::common::ParameterSet parset =
      dp3::steps::test::CreateParameterSet(parameters);
  const DDECal ddecal(&input, parset, "prefix.");
  BOOST_CHECK_EQUAL(expected, dp3::steps::test::Show(ddecal));

  BOOST_CHECK_MESSAGE(parset.unusedKeys().empty(),
                      "Not all keys are used, is there a typo?");
}

BOOST_AUTO_TEST_CASE(show_default) {
  TestShow(
      R"(DDECal prefix.
  mode (constraints):  diagonal
  algorithm:           directionsolve
  H5Parm:              tDDECal.MS/instrument.h5
  write sol to buffer: false
  solint:              1
  nchan:               1
  directions:          [[center]]
  sols per direction:  [1]
  tolerance:           0.0001
  max iter:            50
  flag unconverged:    false
     diverged only:    false
  propagate solutions: false
       converged only: false
  detect stalling:     true
  step size:           0.2
  approximate fitter:  false
  only predict:        false
  subtract model:      false
Model steps for direction center
Predict
OnePredict prefix.
  sourcedb:                tDDECal.MS/sky.txt
   number of patches:      1
   patches clustered:      false
   number of components:   1
   absolute orientation:   false
   all unpolarized:        true
   correct freq smearing:  false
  apply beam:              false
  operation:               replace
  threads:                 )" +
          std::to_string(aocommon::system::ProcessorCount()) + R"(

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky.txt"}});
}

BOOST_AUTO_TEST_CASE(show_modified) {
  TestShow(
      R"(DDECal prefix.
  mode (constraints):  diagonal
  algorithm:           hybrid
  H5Parm:              tDDECal.MS/instrument.h5
  write sol to buffer: false
  solint:              42
  nchan:               44
  directions:          [[center]]
  sols per direction:  [1]
  min visib. ratio:    43.123
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
  approximate fitter:  false
  only predict:        true
  subtract model:      true
Model steps for direction center
Predict
OnePredict prefix.
  sourcedb:                tDDECal.MS/sky.txt
   number of patches:      1
   patches clustered:      false
   number of components:   1
   absolute orientation:   false
   all unpolarized:        true
   correct freq smearing:  false
  apply beam:              false
  operation:               replace
  threads:                 )" +
          std::to_string(aocommon::system::ProcessorCount()) + R"(

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky.txt"},
       {"prefix.propagatesolutions", "true"},
       {"prefix.propagateconvergedonly", "true"},
       {"prefix.flagunconverged", "true"},
       {"prefix.flagdivergedonly", "true"},
       {"prefix.onlypredict", "true"},
       {"prefix.subtract", "true"},
       {"prefix.solveralgorithm", "hybrid"},
       {"prefix.solint", "42"},
       {"prefix.minvisratio", "43.123"},
       {"prefix.nchan", "44"},
       {"prefix.coreconstraint", "45.123"},
       {"prefix.smoothnessconstraint", "46.123"},
       {"prefix.smoothnessreffrequency", "47.123"},
       {"prefix.smoothnessrefdistance", "48.123"},
       {"prefix.tecscreen.coreconstraint", "49.123"},
       {"prefix.maxiter", "49"}});
}

BOOST_AUTO_TEST_SUITE_END()
