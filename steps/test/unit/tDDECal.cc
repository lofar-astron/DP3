// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../DDECal.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "tStepCommon.h"
#include "../../MSReader.h"
#include "../../MultiResultStep.h"

using dp3::steps::DDECal;
using dp3::steps::test::CreateParameterSet;

BOOST_AUTO_TEST_SUITE(ddecal)

BOOST_AUTO_TEST_CASE(provided_fields) {
  using dp3::steps::Step;

  const dp3::common::ParameterSet parset = CreateParameterSet(
      {{"msin", "TODO(AST-1271): Remove msin"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky.txt"},
       // Errors occur when creating multiple DDECal's that write the same H5
       // file. Putting each DDECal object in a different scope also works.
       {"prefix.h5parm", "tddecal_fields.h5"}});
  const DDECal ddecal(parset, "prefix.");
  BOOST_TEST(ddecal.getProvidedFields() == dp3::common::Fields());

  dp3::common::ParameterSet parset_only_predict = parset;
  parset_only_predict.add("prefix.onlypredict", "true");
  parset_only_predict.replace("prefix.h5parm",
                              "tddecal_fields_only_predict.h5");
  const DDECal ddecal_only_predict(parset_only_predict, "prefix.");
  BOOST_TEST(ddecal_only_predict.getProvidedFields() == Step::kDataField);

  dp3::common::ParameterSet parset_subtract = parset;
  parset_subtract.add("prefix.subtract", "true");
  parset_subtract.replace("prefix.h5parm", "tddecal_fields_subtract.h5");
  const DDECal ddecal_subtract(parset_subtract, "prefix.");
  BOOST_TEST(ddecal_subtract.getProvidedFields() == Step::kDataField);
}

/// Helper forwarder to improve test failure output.
static void TestShow(
    const std::string& expected,
    const std::vector<std::pair<std::string, std::string>>& parameters) {
  const dp3::common::ParameterSet parset = CreateParameterSet(parameters);
  const DDECal ddecal(parset, "prefix.");
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

BOOST_DATA_TEST_CASE(store_solutions_in_buffer,
                     boost::unit_test::data::xrange(1, 5), solution_interval) {
  const size_t kNTimes = 6;  // tDDECal.MS has 6 time slots.
  const size_t kNChannelBlocks = 8;
  const size_t kNSolutionsPerChannelBlock = 16;

  auto in = std::make_shared<dp3::steps::MSReader>(
      casacore::MeasurementSet("tDDECal.MS"), dp3::common::ParameterSet(), "");

  auto ddecal = std::make_shared<DDECal>(
      CreateParameterSet({{"msin", "TODO(AST-1271): Remove msin"},
                          {"directions", "[[center]]"},
                          {"sourcedb", "tDDECal.MS/sky.txt"},
                          {"h5parm", "tddecal_storebuffer.h5"},
                          {"storebuffer", "true"},
                          {"solint", std::to_string(solution_interval)}}),
      "");

  auto out = std::make_shared<dp3::steps::MultiResultStep>(kNTimes);

  in->setFieldsToRead(ddecal->getRequiredFields());

  dp3::steps::test::Execute({in, ddecal, out});

  BOOST_REQUIRE(out->get().size() == kNTimes);
  for (size_t i = 0; i < kNTimes; ++i) {
    const std::vector<std::vector<std::complex<double>>>& solution =
        out->get()[i].GetSolution();

    if (i % solution_interval == 0) {
      // Only the first buffer in each solution interval has solutions.
      // Check that the solutions have the expected shape.
      // test_oneapplycal_from_buffer in tDDECal.py checks if the solution
      // values in DPBuffer are correct.
      BOOST_REQUIRE_EQUAL(solution.size(), kNChannelBlocks);
      for (const std::vector<std::complex<double>>& channel_block_solution :
           solution) {
        BOOST_CHECK_EQUAL(channel_block_solution.size(),
                          kNSolutionsPerChannelBlock);
      }
    } else {
      // Other buffers should not have solutions.
      BOOST_CHECK(solution.empty());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
