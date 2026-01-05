// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "steps/MsReader.h"
#include "steps/MultiMsReader.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <aocommon/xt/utensor.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/TableCopy.h>
#include "base/DP3.h"

#include "mock/MockStep.h"
#include "common/ParameterSet.h"
#include "common/test/unit/fixtures/fDirectory.h"

using casacore::Table;

using dp3::base::DPBuffer;
using dp3::common::ParameterSet;
using dp3::common::test::FixtureDirectory;
using dp3::steps::MsReader;
using dp3::steps::MultiMsReader;

const std::string kInputMs = "../tNDPPP_tmp.MS";
const std::string kCopyMs = "tNDPPP_tmp.copy.MS";
const std::string kCopyMsSplit = "tNDPPP_tmp.copy_chunk";  // prefix part only
const std::string kCopyMsPol = "tNDPPP_tmp.copy.MS/POLARIZATION";

// Expected properties of the main data buffer. Were determined with taql with:
// select gmean(DATA), gmin(real(DATA)), gmin(imag(DATA)) from tNDPPP_tmp.MS
//     ...gmax(real(DATA)), gmax(imag(DATA)), gvariance(DATA)...
const std::complex<double> kExpectedDataMean(1.45104, 0.0262766);
const std::complex<float> kExpectedDataMin(-0.325074, -1.3423);
const std::complex<float> kExpectedDataMax(4.00101, 1.46311);
// Except for the variance, as taql computes it differently than xtensor.
const std::complex<double> kExpectedDataVariance(2.10023, 0.0383692);
const xt::xtensor<double, 1> kExpectedUVW = {-0.104612, -639.098, 331.494};
const size_t kNTimeslots = 20;

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

class FixtureSplitChannelCopyExtraData : FixtureCopyAndAddExtraData {
 public:
  FixtureSplitChannelCopyExtraData() : FixtureCopyAndAddExtraData() {
    const std::string kParsetFile = "tDP3.parset";

    // Split input MS into four frequency chunks (to test the MultiMsReader).
    for (std::size_t startchan = 0; startchan < 16; startchan += 4) {
      std::string msout = kCopyMsSplit + std::to_string(startchan / 4) + ".MS";
      {
        std::ofstream ostr(kParsetFile, std::ios::trunc);
        ostr << "msin=" << kInputMs << "\n";
        ostr << "msout=" << msout << "\n";
        ostr << "msin.startchan=" << startchan << "\n";
        ostr << "msin.nchan=4\n";
        ostr << "steps=[]\n";
        ostr << "verbosity=quiet\n";
      }
      dp3::base::Execute(kParsetFile);

      // Add extra data
      casacore::tableCommand(
          "ALTER TABLE " + msout +
          " ADD COLUMN MODEL_DATA C4 shape=[4,4], MODEL_DATA_2 C4 shape=[4,4]");
      casacore::tableCommand("UPDATE " + msout + " SET MODEL_DATA=2*DATA+5");
      casacore::tableCommand("UPDATE " + msout + " SET MODEL_DATA_2=3*DATA+8");
    }
  }
};

/**
 * Collect all data from all DPBuffers gathered by the trailing result step in
 * the chain of (2) steps from the tests. Collects the main data buffer as well
 * as two extra data buffers when @p has_extra_data is set. It concatenates each
 * into a single xtensor cube. Once concatenated, it is verified that main data
 * buffer is as expected. If @p has_extra_data is set, it also verifies that the
 * two extra data buffers satisfy the criteria the extra data has as compared to
 * the main data buffer. So, (2*data+5, and 3*data+8) as set by taql in the two
 * extra-data fixtures.
 *
 * @param has_extra_data Whether or not there are extra data buffers, or that
 *                       there is just the main data buffer.
 */
void ConcatenateAndCheckAllData(
    const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers,
    const bool has_extra_data) {
  BOOST_REQUIRE(!result_buffers.empty());

  // Collect all data. Concatenate along the n_baselines dimension. Do this per
  // DPBuffer at a time for simplicity.
  DPBuffer::DataType all_results = result_buffers[0]->GetData();
  DPBuffer::DataType all_results_m1;
  DPBuffer::DataType all_results_m2;
  if (has_extra_data) {
    all_results_m1 = result_buffers[0]->GetData("MODEL_DATA");
    all_results_m2 = result_buffers[0]->GetData("MODEL_DATA_2");
  }

  for (std::size_t i = 1; i < result_buffers.size(); i++) {
    DPBuffer::DataType curdata = result_buffers[i]->GetData();
    // Exclude buffers from missing time slots (all zeros), this will cause
    // mismatches with the embedded relation between these data columns.
    if (xt::sum(curdata)() != std::complex<float>()) {
      all_results = xt::concatenate(xt::xtuple(all_results, curdata), 0);
      if (has_extra_data) {
        all_results_m1 = xt::concatenate(
            xt::xtuple(all_results_m1,
                       result_buffers[i]->GetData("MODEL_DATA")),
            0);
        all_results_m2 = xt::concatenate(
            xt::xtuple(all_results_m2,
                       result_buffers[i]->GetData("MODEL_DATA_2")),
            0);
      }
    }
  }

  // Test main data buffer.
  BOOST_CHECK_CLOSE(xt::mean(all_results)(), kExpectedDataMean, 1e-4);
  BOOST_TEST(xt::amin(xt::real(all_results))() == kExpectedDataMin.real());
  BOOST_TEST(xt::amin(xt::imag(all_results))() == kExpectedDataMin.imag());
  BOOST_TEST(xt::amax(xt::real(all_results))() == kExpectedDataMax.real());
  BOOST_TEST(xt::amax(xt::imag(all_results))() == kExpectedDataMax.imag());
  BOOST_CHECK_CLOSE(xt::variance(all_results)(), kExpectedDataVariance, 1e-4);

  // Test extra data buffers
  if (has_extra_data) {
    BOOST_CHECK(
        xt::allclose(2 * all_results + std::complex<float>(5), all_results_m1));
    BOOST_CHECK(
        xt::allclose(3 * all_results + std::complex<float>(8), all_results_m2));
  }
}

BOOST_FIXTURE_TEST_CASE(missing_data, FixtureDirectory) {
  const casacore::MeasurementSet ms(kInputMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kInputMs);
  parset.add("msin.datacolumn", "THISDOESNOTEXIST");
  BOOST_CHECK_THROW(MsReader reader(ms, parset, "msin."), std::runtime_error);

  // With missingData set, it should just print a warning.
  MsReader reader(ms, parset, "msin.", true);
}

BOOST_FIXTURE_TEST_CASE(missing_extra_data, FixtureCopyAndAddExtraData) {
  // Single non-existing extra-data column.
  const casacore::MeasurementSet ms(kCopyMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset1;
  parset1.add("msin", kCopyMs);
  parset1.add("msin.extradatacolumns", "[MODEL_DATA, NONEXISTING]");
  BOOST_CHECK_THROW(MsReader reader1(ms, parset1, "msin."), std::runtime_error);

  // Two non-existing extra-data columns.
  ParameterSet parset2;
  parset2.add("msin", kCopyMs);
  parset2.add("msin.extradatacolumns", "[THISDOESNOTEXIST, NONEXISTINGTOO]");
  BOOST_CHECK_THROW(MsReader reader2(ms, parset2, "msin."), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(starttimeslot, FixtureDirectory) {
  // kInputMs is irregular at the beginning, which allows testing various cases.
  const std::vector<double> kSlotTimes = {
      4472025680.0,  // Extra time slot at the beginning.
      4472025710.0,  // Extra time slot at the beginning.
      4472025740.0,  // First actual time slot of kInputMs.
      4472025765.0,  // Second time slot, 25 seconds (!) later.
      4472025795.0,  // Third time slot, 30 seconds later.
      4472025885.0,  // Fourth time slot, 90 seconds (!) later.
      4472025915.0   // Fifth time slot, 30 seconds later.
  };
  for (int i = 0; i < int(kSlotTimes.size()); ++i) {
    const int starttimeslot = i - 2;
    const casacore::MeasurementSet ms(kInputMs,
                                      casacore::TableLock::AutoNoReadLocking);
    ParameterSet parset;
    parset.add("msin.starttimeslot", std::to_string(starttimeslot));
    const MsReader reader(ms, parset, "msin.");
    BOOST_CHECK_CLOSE(reader.getInfoOut().firstTime(), kSlotTimes[i], 1.0e-10);
  }
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

  MsReader reader(ms, parset, "msin.");
  reader.setFieldsToRead(dp3::steps::Step::kDataField |
                         dp3::steps::Step::kWeightsField |
                         dp3::steps::Step::kUvwField);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  // Read the MS.
  reader.updateInfo(dp3::base::DPInfo());
  size_t count = 0;
  while (reader.process(std::make_unique<dp3::base::DPBuffer>())) count++;
  reader.finish();
  BOOST_REQUIRE(count == kNTimeslots);

  const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers =
      mock_step->GetRegularBuffers();

  // Test data buffer.
  ConcatenateAndCheckAllData(result_buffers, false);

  // Test UVW and weights
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

  MsReader reader(ms, parset, "msin.");
  reader.setFieldsToRead(dp3::steps::Step::kDataField);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  // Read the MS.
  reader.updateInfo(dp3::base::DPInfo());
  size_t count = 0;
  while (reader.process(std::make_unique<dp3::base::DPBuffer>())) count++;
  reader.finish();
  BOOST_REQUIRE(count == kNTimeslots);

  const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers =
      mock_step->GetRegularBuffers();

  ConcatenateAndCheckAllData(result_buffers, true);
}

BOOST_FIXTURE_TEST_CASE(polarization_initialization, FixtureDirectory) {
  const casacore::MeasurementSet ms(kInputMs,
                                    casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kInputMs);
  const MsReader reader(ms, parset, "msin.");

  const std::set<aocommon::PolarizationEnum> expected_polarizations{
      aocommon::PolarizationEnum::XX, aocommon::PolarizationEnum::XY,
      aocommon::PolarizationEnum::YX, aocommon::PolarizationEnum::YY};

  const std::set<aocommon::PolarizationEnum> actual_polarizations =
      reader.getInfoOut().polarizations();

  BOOST_CHECK_EQUAL_COLLECTIONS(
      actual_polarizations.begin(), actual_polarizations.end(),
      expected_polarizations.begin(), expected_polarizations.end());
}

BOOST_FIXTURE_TEST_CASE(expect_four_polarizations, FixtureCopyAndUpdatePol) {
  casacore::MeasurementSet ms(kCopyMs, casacore::TableLock::AutoNoReadLocking);
  ParameterSet parset;
  parset.add("msin", kCopyMs);

  BOOST_CHECK_THROW(std::make_unique<MsReader>(ms, parset, "msin."),
                    std::runtime_error);
}

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(0.0001) *
                     boost::unit_test::tolerance(0.0001f))
BOOST_DATA_TEST_CASE_F(FixtureSplitChannelCopyExtraData, process_multiple_ms,
                       boost::unit_test::data::xrange(2),
                       test_with_extra_data) {
  // Bundle input MS names for parset and MultiMsReader constructor
  std::vector<std::string> msNames;
  std::string msNames_string;
  for (std::size_t startchan = 0; startchan < 16; startchan += 4) {
    std::string ms_in = kCopyMsSplit + std::to_string(startchan / 4) + ".MS";
    msNames.emplace_back(ms_in);
    msNames_string += ms_in + ",";
  }
  msNames_string.erase(msNames_string.size() - 1);

  // Make step chain.
  ParameterSet parset;
  parset.add("msin", '[' + msNames_string + ']');
  if (test_with_extra_data)
    parset.add("msin.extradatacolumns", "[MODEL_DATA, MODEL_DATA_2]");

  MultiMsReader reader(msNames, parset, "msin.");
  reader.setFieldsToRead(dp3::steps::Step::kDataField);

  auto mock_step = std::make_shared<dp3::steps::MockStep>();
  reader.setNextStep(mock_step);

  // Read the MS.
  reader.updateInfo(dp3::base::DPInfo());
  size_t count = 0;
  while (reader.process(std::make_unique<dp3::base::DPBuffer>())) count++;
  reader.finish();
  BOOST_REQUIRE(count == kNTimeslots);

  const std::vector<std::unique_ptr<dp3::base::DPBuffer>>& result_buffers =
      mock_step->GetRegularBuffers();

  ConcatenateAndCheckAllData(result_buffers, test_with_extra_data);
}

BOOST_AUTO_TEST_SUITE_END()
