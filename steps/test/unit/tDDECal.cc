// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../DDECal.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

#include <dp3/base/DP3.h>

#include "../../MSReader.h"
#include "../../MultiResultStep.h"
#include "../../ResultStep.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"
#include "tStepCommon.h"

using dp3::common::test::FixtureDirectory;
using dp3::steps::DDECal;
using dp3::steps::test::CreateParameterSet;

BOOST_AUTO_TEST_SUITE(ddecal)

BOOST_AUTO_TEST_CASE(provided_fields) {
  using dp3::steps::Step;

  const dp3::common::ParameterSet parset = CreateParameterSet(
      {{"prefix.directions", "[[center]]"},
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
  keep model:          false
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
  directions:          [[REUSE],[center]]
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
  keep model:          true
Direction REUSE reuses data from otherddecal.REUSE
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
       {"prefix.keepmodel", "true"},
       {"prefix.reusemodel", "[otherddecal.REUSE]"},
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
      CreateParameterSet({{"directions", "[[center]]"},
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

namespace {

// Fixture with common code for the keep_model_data tests.
class KeepModelDataFixture : public FixtureDirectory {
 public:
  KeepModelDataFixture()
      : reader(std::make_shared<dp3::steps::MSReader>(
            casacore::MeasurementSet(kMsName), dp3::common::ParameterSet(),
            "")),
        result_step(std::make_shared<dp3::steps::MultiResultStep>(kNTimes)) {}

  void Execute(std::shared_ptr<DDECal>& ddecal) {
    reader->setFieldsToRead(ddecal->getRequiredFields());
    dp3::steps::test::Execute({reader, ddecal, result_step});
  }

  void CheckOutput(const std::vector<std::string>& expected_names) {
    BOOST_CHECK(result_step->get().size() == kNTimes);
    for (const dp3::base::DPBuffer& buffer : result_step->get()) {
      for (const std::string& name : expected_names) {
        BOOST_REQUIRE(buffer.HasData(name));
        const std::array<std::size_t, 3> shape = buffer.GetData(name).shape();
        BOOST_CHECK_EQUAL(shape[0], reader->getInfo().nbaselines());
        BOOST_CHECK_EQUAL(shape[1], reader->getInfo().nchan());
        BOOST_CHECK_EQUAL(shape[2], reader->getInfo().ncorr());
      }
    }
  }

 public:
  const std::string kPrefix = "FooPrefix.";
  const std::string kMsName = "../tDDECal.MS";

 private:
  const std::size_t kNTimes = 6;  // tDDECal.MS has 6 time slots.
  std::shared_ptr<dp3::steps::MSReader> reader;
  std::shared_ptr<dp3::steps::MultiResultStep> result_step;
};
}  // namespace

// Because of AST-1281, use separate tests for "keepmodel": One with
// "modeldatacolumns", one with "directions", and one with "idg.regions".
BOOST_FIXTURE_TEST_CASE(keep_model_data_columnreader, KeepModelDataFixture) {
  // tDDECal.MS has an extra foursources_DATA column (see tIDGPredict.py).
  const std::string kModelDataColumn = "foursources_DATA";

  auto ddecal = std::make_shared<DDECal>(
      CreateParameterSet({{kPrefix + "modeldatacolumns", kModelDataColumn},
                          {kPrefix + "keepmodel", "true"},
                          {kPrefix + "h5parm", "keepmodel_columnreader.h5"}}),
      kPrefix);
  Execute(ddecal);
  CheckOutput({kPrefix + kModelDataColumn});
}

BOOST_FIXTURE_TEST_CASE(keep_model_data_directions, KeepModelDataFixture) {
  auto ddecal = std::make_shared<DDECal>(
      CreateParameterSet(
          {{kPrefix + "sourcedb", kMsName + "/sky.txt"},
           {kPrefix + "directions", "[[center,dec_off],[ra_off],[radec_off]]"},
           {kPrefix + "keepmodel", "true"},
           {kPrefix + "h5parm", "keepmodel_directions.h5"}}),
      kPrefix);
  Execute(ddecal);

  // The expected names in the output buffer are the prefix followed by the
  // names of the first directions.
  CheckOutput({kPrefix + "center", kPrefix + "ra_off", kPrefix + "radec_off"});
}

BOOST_FIXTURE_TEST_CASE(keep_model_data_idg, KeepModelDataFixture) {
  auto ddecal = std::make_shared<DDECal>(
      CreateParameterSet({{kPrefix + "idg.regions", "../sources.reg"},
                          {kPrefix + "idg.images", "../sources-model.fits"},
                          {kPrefix + "keepmodel", "true"},
                          {kPrefix + "h5parm", "keepmodel_idg.h5"}}),
      kPrefix);
  Execute(ddecal);

  // The expected names in the output buffer are the prefix followed by "dir"
  // followed by an increasing number.
  CheckOutput(
      {kPrefix + "dir0", kPrefix + "dir1", kPrefix + "dir2", kPrefix + "dir3"});
}

BOOST_FIXTURE_TEST_CASE(model_data_is_corrected, FixtureDirectory) {
  const std::string kModelName{"model name"};
  const std::size_t kNCorrelations{4};
  const std::size_t kNChannels{1};
  const std::size_t kNStations{5};
  const std::size_t kNBaselines{(kNStations * (kNStations - 1)) / 2};
  const std::array<std::size_t, 3> kShape{kNBaselines, kNChannels,
                                          kNCorrelations};
  const std::vector<std::string> kAntennaNames(kNStations, "");
  const std::vector<double> kAntennaDiameters(kNStations, 1);
  const std::vector<casacore::MPosition> kAntennaPositions(
      kNStations, casacore::MPosition());
  const std::complex<float> kDataValue{2.0f, 2.0f};
  const std::complex<float> kStation0DataValue{4.0f, -4.0f};
  const std::complex<float> kModelDataValue{8.0f, 8.0f};
  const std::vector<std::complex<double>> kExpectedSolution{
      {0.32391142, -0.94158080},
      {0.47395513, 0.16304438},
      {0.47395513, 0.16304438},
      {0.47395513, 0.16304438},
      {0.47395513, 0.16304438}};
  // Proof that kExpectedSolution is correct:
  BOOST_CHECK_CLOSE(kExpectedSolution[0] *
                        std::complex<double>(kModelDataValue) *
                        std::conj(kExpectedSolution[1]),
                    std::complex<double>(kStation0DataValue), 1.0);
  BOOST_CHECK_CLOSE(kExpectedSolution[1] *
                        std::complex<double>(kModelDataValue) *
                        std::conj(kExpectedSolution[1]),
                    std::complex<double>(kDataValue), 1.0);

  std::vector<int> antenna1;
  std::vector<int> antenna2;
  for (std::size_t station1 = 0; station1 < kNStations; ++station1) {
    for (std::size_t station2 = station1 + 1; station2 < kNStations;
         ++station2) {
      antenna1.push_back(static_cast<int>(station1));
      antenna2.push_back(static_cast<int>(station2));
    }
  }

  dp3::base::DPInfo info(kNCorrelations, kNChannels);
  info.setAntennas(kAntennaNames, kAntennaDiameters, kAntennaPositions,
                   antenna1, antenna2);
  info.setChannels(std::vector<double>(kNChannels, 42.0e6),
                   std::vector<double>(kNChannels, 1.0e6));
  info.setNThreads(1);

  auto ddecal = std::make_shared<DDECal>(
      CreateParameterSet({{"keepmodel", "true"},
                          {"reusemodel", "[" + kModelName + "]"},
                          {"subtract", "false"},
                          {"storebuffer", "true"},
                          {"mode", "scalar"},
                          {"h5parm", "h5parm_output_is_mandatory.h5"}}),
      "");
  auto result_step = std::make_shared<dp3::steps::ResultStep>();
  ddecal->setNextStep(result_step);
  ddecal->setInfo(info);

  auto buffer = std::make_unique<dp3::base::DPBuffer>();

  xt::xtensor<std::complex<float>, 3> input_data(kShape, kDataValue);
  xt::view(input_data, xt::range(0, kNStations - 1), xt::all(), xt::all())
      .fill(kStation0DataValue);

  buffer->ResizeData(kShape);
  buffer->GetData().assign(input_data);
  buffer->AddData(kModelName);
  buffer->GetData(kModelName).fill(kModelDataValue);

  // Use default flags and weights.
  buffer->ResizeFlags(kShape);
  buffer->GetFlags().fill(false);
  buffer->ResizeWeights(kShape);
  buffer->GetWeights().fill(1.0f);

  ddecal->process(std::move(buffer));
  ddecal->finish();

  buffer = result_step->take();
  BOOST_REQUIRE(buffer);

  // Check model data buffer. It should be close to input data buffer.
  BOOST_REQUIRE(buffer->HasData(kModelName));
  BOOST_CHECK(xt::allclose(buffer->GetData(kModelName), input_data, 1.0e-2));

  // Check that (main) data, weights and flags remained equal.
  BOOST_CHECK_EQUAL(buffer->GetData(), input_data);
  BOOST_CHECK_EQUAL(buffer->GetFlags(), (xt::xtensor<bool, 3>(kShape, false)));
  BOOST_CHECK_EQUAL(buffer->GetWeights(),
                    (xt::xtensor<float, 3>(kShape, 1.0f)));

  // Check the gains in the solutions.
  const std::vector<std::vector<std::complex<double>>>& solution =
      buffer->GetSolution();
  BOOST_REQUIRE_EQUAL(solution.size(), 1);
  BOOST_REQUIRE_EQUAL(solution.front().size(), kNStations);
  for (size_t i = 0; i < kNStations; ++i) {
    BOOST_CHECK_CLOSE(solution.front()[i], kExpectedSolution[i], 1.0e-3);
  }
}

BOOST_AUTO_TEST_SUITE_END()
