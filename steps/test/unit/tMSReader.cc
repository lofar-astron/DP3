// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSReader.h"

#include <memory>

#include <boost/test/unit_test.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <aocommon/xt/utensor.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/TableCopy.h>

#include "mock/MockStep.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"

using casacore::Table;

using dp3::base::DPBuffer;
using dp3::common::ParameterSet;
using dp3::common::test::FixtureDirectory;
using dp3::steps::MSReader;

const std::string kInputMs = "../tNDPPP_tmp.MS";
const std::string kCopyMs = "tNDPPP_tmp.copy.MS";
const std::string kCopyMsPol = "tNDPPP_tmp.copy.MS/POLARIZATION";

BOOST_AUTO_TEST_SUITE(msreader)

class FixtureCopyAndUpdatePol : FixtureDirectory {
 public:
  FixtureCopyAndUpdatePol() : FixtureDirectory() {
    Table(kInputMs).deepCopy(kCopyMs, Table::New);
    casacore::tableCommand("update " + kCopyMsPol +
                           " set CORR_TYPE=[9, 12], "
                           "CORR_PRODUCT=[[0,0],[1,1]], NUM_CORR=2");
  }
};

class FixtureCopyAndAddExtraData : FixtureDirectory {
 public:
  FixtureCopyAndAddExtraData() : FixtureDirectory() {
    Table(kInputMs).deepCopy(kCopyMs, Table::New);
    casacore::tableCommand(
        "ALTER TABLE " + kCopyMs +
        " ADD COLUMN MODEL_DATA C4 shape=[4,16], MODEL_DATA_2 C4 shape=[4,16]");
    casacore::tableCommand("UPDATE " + kCopyMs + " SET MODEL_DATA=2*DATA+5");
    casacore::tableCommand("UPDATE " + kCopyMs + " SET MODEL_DATA_2=3*DATA+8");
  }
};

BOOST_FIXTURE_TEST_CASE(missing_data, FixtureDirectory) {
  const casacore::MeasurementSet ms(kInputMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kInputMs);
  parset.add("msin.datacolumn", "THISDOESNOTEXIST");
  BOOST_CHECK_THROW(MSReader reader(ms, parset, "msin."), std::runtime_error);

  // With missingData set, it should just print a warning.
  MSReader reader(ms, parset, "msin.", true);
}

BOOST_FIXTURE_TEST_CASE(missing_extra_data, FixtureCopyAndAddExtraData) {
  // Single non-existing extra-data column.
  const casacore::MeasurementSet ms(kCopyMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset1;
  parset1.add("msin", kCopyMs);
  parset1.add("msin.extradatacolumns", "[MODEL_DATA, NONEXISTING]");
  BOOST_CHECK_THROW(MSReader reader1(ms, parset1, "msin."), std::runtime_error);

  // Two non-existing extra-data columns.
  ParameterSet parset2;
  parset2.add("msin", kCopyMs);
  parset2.add("msin.extradatacolumns", "[THISDOESNOTEXIST, NONEXISTINGTOO]");
  BOOST_CHECK_THROW(MSReader reader2(ms, parset2, "msin."), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(process, FixtureDirectory,
                        *boost::unit_test::tolerance(0.0001) *
                            boost::unit_test::tolerance(0.0001f)) {
  // Make step chain.
  const casacore::MeasurementSet ms(kInputMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kInputMs);
  parset.add("msin.weightcolumn", "WEIGHT");

  MSReader reader(ms, parset, "msin.");
  reader.setFieldsToRead(dp3::steps::Step::kDataField |
                         dp3::steps::Step::kWeightsField |
                         dp3::steps::Step::kUvwField);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  // Read the MS.
  reader.updateInfo(dp3::base::DPInfo());
  while (reader.process(std::make_unique<dp3::base::DPBuffer>()))
    ;
  reader.finish();

  const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers =
      mock_step->GetRegularBuffers();
  BOOST_REQUIRE(!result_buffers.empty());

  // The following constants were determined with taql with:
  // select gmean(DATA), gmin(real(DATA)), gmin(imag(DATA)) from tNDPPP_tmp.MS
  //     ...gmax(real(DATA)), gmax(imag(DATA)), gvariance(DATA)...
  std::complex<double> kExpectedDataMean(1.45104, 0.0262766);
  std::complex<float> kExpectedDataMin(-0.325074, -1.3423);
  std::complex<float> kExpectedDataMax(4.00101, 1.46311);
  // Except for the variance, seems taql computes it differently than xtensor
  std::complex<double> kExpectedDataVariance(2.10023, 0.0383692);
  xt::xtensor<double, 1> kExpectedUVW = {-0.104612, -639.098, 331.494};

  // Collect all data. Concatenate along the n_baselines dimension. Do this per
  // DPBuffer at a time for simplicity.
  DPBuffer::DataType all_results = result_buffers[0]->GetData();
  for (std::size_t i = 1; i < result_buffers.size(); i++) {
    DPBuffer::DataType curdata = result_buffers[i]->GetData();
    // Exclude buffers from missing time slots (all zeros), otherwise this will
    // cause mismatches with the taql statistics.
    if (xt::sum(curdata)() != std::complex<float>())
      all_results = xt::concatenate(xt::xtuple(all_results, curdata), 0);
  }

  // Tests.
  BOOST_CHECK_CLOSE(xt::mean(all_results)(), kExpectedDataMean, 1e-4);
  BOOST_TEST(xt::amin(xt::real(all_results))() == kExpectedDataMin.real());
  BOOST_TEST(xt::amin(xt::imag(all_results))() == kExpectedDataMin.imag());
  BOOST_TEST(xt::amax(xt::real(all_results))() == kExpectedDataMax.real());
  BOOST_TEST(xt::amax(xt::imag(all_results))() == kExpectedDataMax.imag());
  BOOST_CHECK_CLOSE(xt::variance(all_results)(), kExpectedDataVariance, 1e-4);

  BOOST_CHECK(
      xt::allclose(xt::row(result_buffers[0]->GetUvw(), 5), kExpectedUVW));
  for (std::size_t i = 0; i < result_buffers.size(); i++) {
    DPBuffer::WeightsType curweights = result_buffers[i]->GetWeights();
    if (xt::sum(result_buffers[i]->GetData())() == std::complex<float>()) {
      // Buffer from a missing time slot
      BOOST_TEST(curweights ==
                 xt::zeros<std::complex<float>>(curweights.shape()));
    } else {
      BOOST_TEST(curweights ==
                 xt::ones<std::complex<float>>(curweights.shape()));
    }
  }
}

BOOST_FIXTURE_TEST_CASE(process_extra_data, FixtureCopyAndAddExtraData,
                        *boost::unit_test::tolerance(0.0001) *
                            boost::unit_test::tolerance(0.0001f)) {
  // Make step chain.
  const casacore::MeasurementSet ms(kCopyMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kCopyMs);
  parset.add("msin.extradatacolumns", "[MODEL_DATA, MODEL_DATA_2]");

  MSReader reader(ms, parset, "msin.");
  reader.setFieldsToRead(dp3::steps::Step::kDataField);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  // Read the MS.
  reader.updateInfo(dp3::base::DPInfo());
  while (reader.process(std::make_unique<dp3::base::DPBuffer>()))
    ;
  reader.finish();

  const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers =
      mock_step->GetRegularBuffers();
  BOOST_REQUIRE(!result_buffers.empty());

  // Collect all data. Concatenate along the n_baselines dimension. Do this per
  // DPBuffer at a time for simplicity.
  using BufferDataType = DPBuffer::DataType;
  BufferDataType all_results = result_buffers[0]->GetData();
  BufferDataType all_results_m1 = result_buffers[0]->GetData("MODEL_DATA");
  BufferDataType all_results_m2 = result_buffers[0]->GetData("MODEL_DATA_2");
  for (std::size_t i = 1; i < result_buffers.size(); i++) {
    BufferDataType curdata = result_buffers[i]->GetData();
    BufferDataType curdata_m1 = result_buffers[i]->GetData("MODEL_DATA");
    BufferDataType curdata_m2 = result_buffers[i]->GetData("MODEL_DATA_2");
    // Exclude buffers from missing time slots (all zeros), this will cause
    // mismatches with the embedded relation between these data columns.
    if (xt::sum(curdata)() != std::complex<float>()) {
      all_results = xt::concatenate(xt::xtuple(all_results, curdata), 0);
      all_results_m1 =
          xt::concatenate(xt::xtuple(all_results_m1, curdata_m1), 0);
      all_results_m2 =
          xt::concatenate(xt::xtuple(all_results_m2, curdata_m2), 0);
    }
  }

  BOOST_CHECK(
      xt::allclose(2 * all_results + std::complex<float>(5), all_results_m1));
  BOOST_CHECK(
      xt::allclose(3 * all_results + std::complex<float>(8), all_results_m2));
}

BOOST_FIXTURE_TEST_CASE(polarization_initialization, FixtureDirectory) {
  const casacore::MeasurementSet ms(kInputMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kInputMs);
  const MSReader reader(ms, parset, "msin.");

  const std::set<aocommon::PolarizationEnum> expected_polarizations{
      aocommon::PolarizationEnum::XX, aocommon::PolarizationEnum::XY,
      aocommon::PolarizationEnum::YX, aocommon::PolarizationEnum::YY};

  const std::set<aocommon::PolarizationEnum> actual_polarizations =
      reader.getInfo().polarizations();

  BOOST_CHECK_EQUAL_COLLECTIONS(
      actual_polarizations.begin(), actual_polarizations.end(),
      expected_polarizations.begin(), expected_polarizations.end());
}

BOOST_FIXTURE_TEST_CASE(expect_four_polarizations, FixtureCopyAndUpdatePol) {
  casacore::MeasurementSet ms(kCopyMs, casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kCopyMs);

  BOOST_CHECK_THROW(std::make_unique<MSReader>(ms, parset, "msin."),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
