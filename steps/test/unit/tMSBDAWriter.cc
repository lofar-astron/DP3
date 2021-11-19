// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

#include "../../../base/BDABuffer.h"
#include "../../MSBDAWriter.h"
#include "../../MSReader.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"

using dp3::base::BDABuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::MSBDAWriter;
using dp3::steps::MSReader;

using casacore::MeasurementSet;
using casacore::Table;
using casacore::TableDesc;
using casacore::TableExprNode;
using casacore::TableLock;

const std::string prefix = "msout.";

BOOST_AUTO_TEST_SUITE(msbdawriter)

// Test that the output measurementset is correct for simple data.
BOOST_FIXTURE_TEST_CASE(process_simple, FixtureDirectory,
                        *boost::unit_test::label("slow")) {
  const unsigned int ncorr(1);
  const unsigned int nchan(1);
  const double kTime(3.0);
  const double kInterval(1.5);
  const double kExposure(1);
  const std::complex<float> kData(42.0, 43.0);
  const float kWeight(44.0);
  const bool kFlag(false);
  const double kUVW[3]{45.0, 46.0, 47.0};
  const std::string kMsName = "bda_simple.MS";

  // TODO: Remove reader, when the writer can create subtables etc.
  const casacore::MeasurementSet ms_in("../tNDPPP_tmp.MS");
  MSReader reader(ms_in, ParameterSet(), "");
  MSBDAWriter writer(&reader, kMsName, ParameterSet(), "");

  DPInfo info;
  info.init(ncorr, 0, 1, nchan, kTime, kInterval, "", "");
  info.set(std::vector<std::string>{"ant"}, std::vector<double>{1.0},
           {casacore::MVPosition{0, 0, 0}}, std::vector<int>{0},
           std::vector<int>{0});
  info.set(std::vector<std::vector<double>>{{1.}},
           std::vector<std::vector<double>>{{10.}});
  info.set(std::vector<double>(nchan, 1.), std::vector<double>(nchan, 5000.));
  writer.updateInfo(info);

  auto buffer = boost::make_unique<BDABuffer>(1);
  buffer->AddRow(kTime, kInterval, kExposure, 0, 1, 1, &kData, &kFlag, &kWeight,
                 nullptr, kUVW);
  writer.process(std::move(buffer));
  writer.finish();

  // Created MS should be valid
  BOOST_CHECK_NO_THROW(MeasurementSet(kMsName, TableLock::AutoNoReadLocking));

  MeasurementSet ms(kMsName, TableLock::AutoNoReadLocking);

  // Assert the visibility data
  BOOST_TEST(ms.col("TIME").getDouble(0) == kTime);
  BOOST_TEST(ms.col("TIME_CENTROID").getDouble(0) == kTime);
  BOOST_TEST(ms.col("EXPOSURE").getDouble(0) == kExposure);
  BOOST_TEST(ms.col("ANTENNA1").getInt(0) == 0);
  BOOST_TEST(ms.col("ANTENNA2").getInt(0) == 0);
  casacore::IPosition dim(2, ncorr, nchan);
  BOOST_TEST(ms.col("DATA").getArrayDComplex(0).tovector()[0].imag() ==
             kData.imag());
  BOOST_TEST(ms.col("DATA").getArrayDComplex(0).tovector()[0].real() ==
             kData.real());
  BOOST_TEST(ms.col("WEIGHT_SPECTRUM").getArrayDouble(0).tovector()[0] ==
             kWeight);
  BOOST_TEST(ms.col("FLAG").getArrayBool(0).tovector()[0] == kFlag);
  BOOST_TEST(ms.col("UVW").getArrayDouble(0).tovector() ==
             std::vector<double>(std::begin(kUVW), std::end(kUVW)));

  // Assert default values
  BOOST_TEST(ms.col("FEED1").getInt(0) == 0);
  BOOST_TEST(ms.col("FEED2").getInt(0) == 0);
  BOOST_TEST(ms.col("DATA_DESC_ID").getInt(0) == 0);
  BOOST_TEST(ms.col("PROCESSOR_ID").getInt(0) == 0);
  BOOST_TEST(ms.col("INTERVAL").getDouble(0) == kInterval);
  BOOST_TEST(ms.col("SCAN_NUMBER").getInt(0) == 0);
  BOOST_TEST(ms.col("ARRAY_ID").getInt(0) == 0);
  BOOST_TEST(ms.col("OBSERVATION_ID").getInt(0) == 0);
  BOOST_TEST(ms.col("STATE_ID").getInt(0) == 0);

  // Assert sigma and weight
  std::vector<double> sigma_weight(ncorr, 1);
  BOOST_TEST(ms.col("WEIGHT").getArrayDouble(0).tovector() == sigma_weight);
  BOOST_TEST(ms.col("SIGMA").getArrayDouble(0).tovector() == sigma_weight);

  // Assert if the correct columns are created
  Table table(kMsName, TableLock::AutoNoReadLocking);
  BOOST_TEST(table.keywordSet().isDefined("BDA_TIME_AXIS"));
  const TableDesc td_bda =
      table.keywordSet().asTable("BDA_TIME_AXIS").tableDesc();
  BOOST_TEST(td_bda.isColumn("BDA_TIME_AXIS_ID"));
  BOOST_TEST(td_bda.isColumn("IS_BDA_APPLIED"));
  BOOST_TEST(td_bda.isColumn("SINGLE_FACTOR_PER_BASELINE"));
  BOOST_TEST(td_bda.isColumn("MAX_TIME_INTERVAL"));
  BOOST_TEST(td_bda.isColumn("MIN_TIME_INTERVAL"));
  BOOST_TEST(td_bda.isColumn("UNIT_TIME_INTERVAL"));
  BOOST_TEST(td_bda.isColumn("INTEGER_INTERVAL_FACTORS"));
  BOOST_TEST(td_bda.isColumn("HAS_BDA_ORDERING"));
  BOOST_TEST(td_bda.isColumn("BDA_FREQ_AXIS_ID"));
  BOOST_TEST(td_bda.isColumn("FIELD_ID"));

  // Check that the version is present and set
  BOOST_TEST(td_bda.keywordSet().asString("BDA_TIME_AXIS_VERSION") == "1.0");

  // Assert if the correct columns are created
  BOOST_TEST(table.keywordSet().isDefined("SPECTRAL_WINDOW"));
  const Table td = table.keywordSet().asTable("SPECTRAL_WINDOW");
  BOOST_TEST(!td.tableDesc().isColumn("BDA_FREQ_AXIS_ID"));
  BOOST_TEST(td.tableDesc().isColumn("BDA_SET_ID"));

  BOOST_TEST(td.col("BDA_SET_ID").getInt(0) == 0);
  // Assert that the metadata is filled correctly
  Table t = table.keywordSet().asTable("BDA_TIME_AXIS");
  BOOST_TEST(t.nrow() == 1U);
  BOOST_TEST(t.col("IS_BDA_APPLIED").getBool(0));
  BOOST_TEST(t.col("SINGLE_FACTOR_PER_BASELINE").getBool(0));
  // Assert that the interval is correct
  BOOST_TEST(t.col("MIN_TIME_INTERVAL").getDouble(0) == kInterval);
  BOOST_TEST(t.col("MAX_TIME_INTERVAL").getDouble(0) == kInterval);
  BOOST_TEST(t.col("UNIT_TIME_INTERVAL").getDouble(0) == kInterval);
  BOOST_TEST(t.col("INTEGER_INTERVAL_FACTORS").getBool(0));
  BOOST_TEST(t.col("HAS_BDA_ORDERING").getBool(0));
  BOOST_TEST(t.col("BDA_FREQ_AXIS_ID").getInt(0) == -1);
  BOOST_TEST(t.col("FIELD_ID").getInt(0) == -1);
}

// Test that an exception is thrown when the base lines are not correct
BOOST_FIXTURE_TEST_CASE(exception_when_mismatch, FixtureDirectory) {
  const double kTime(3.0);
  const double kInterval(1.5);
  const std::string kMsName = "bda_exception.MS";

  MSBDAWriter writer(nullptr, kMsName, ParameterSet(), "");

  DPInfo info;
  info.init(1, 0, 1, 1, kTime, kInterval, "", "");
  info.set(std::vector<std::string>{"ant", "ant2"},
           std::vector<double>{1.0, 2.0},
           {casacore::MVPosition{0, 0, 0}, casacore::MVPosition{0, 0, 0}},
           std::vector<int>{0, 1}, std::vector<int>{0, 1});
  info.set(std::vector<double>(1, 1.), std::vector<double>(1, 5000.));

  // ntimeAvgs is 1, nbaselines 2 so we expect an exception
  BOOST_TEST(info.ntimeAvgs().size() == size_t(1));
  BOOST_TEST(info.nbaselines() == static_cast<unsigned int>(2));
  BOOST_CHECK_THROW(writer.updateInfo(info), std::invalid_argument);
}

// Test that empty default subtables are created when there is no reader data
// available.
BOOST_FIXTURE_TEST_CASE(create_default_subtables, FixtureDirectory,
                        *boost::unit_test::label("slow")) {
  DPInfo info;
  info.init(1, 0, 1, 1, 3.0, 1.5, "", "");
  info.set(std::vector<std::string>{"ant"}, std::vector<double>{1.0},
           {casacore::MVPosition{0, 0, 0}}, std::vector<int>{0},
           std::vector<int>{0});
  info.set(std::vector<double>(1, 1.), std::vector<double>(1, 5000.));
  const std::string kMsName = "default_tables.MS";

  MSBDAWriter writer(nullptr, kMsName, ParameterSet(), "");
  writer.updateInfo(info);

  Table table(kMsName, TableLock::AutoNoReadLocking);
  BOOST_TEST(table.keywordSet().isDefined("ANTENNA"));
  BOOST_TEST(table.keywordSet().isDefined("DATA_DESCRIPTION"));
  BOOST_TEST(table.keywordSet().isDefined("FEED"));
  BOOST_TEST(table.keywordSet().isDefined("FIELD"));
  BOOST_TEST(table.keywordSet().isDefined("HISTORY"));
  BOOST_TEST(table.keywordSet().isDefined("OBSERVATION"));
  BOOST_TEST(table.keywordSet().isDefined("SPECTRAL_WINDOW"));
  BOOST_CHECK_NO_THROW(
      MeasurementSet ms(kMsName, TableLock::AutoNoReadLocking));
}

// Test BDA_TIME_AXIS for different max and min intervals.
BOOST_FIXTURE_TEST_CASE(different_bda_intervals, FixtureDirectory,
                        *boost::unit_test::label("slow")) {
  // Setup test
  const std::string msOutName = "bda_multiple_ms_out.MS";
  const double timeInterval = 1;
  const unsigned int kMinTimeInterval = 2 * timeInterval;
  const unsigned int kMaxTimeInterval = 10 * timeInterval;
  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  info.init(1, 0, 1, 1, 3.0, timeInterval, "", "");
  info.set(std::vector<std::string>{"ant", "ant2"},
           std::vector<double>{1.0, 1.0},
           {casacore::MVPosition{0, 0, 0}, casacore::MVPosition{10, 10, 0}},
           std::vector<int>{0, 1}, std::vector<int>{0, 1});
  info.set(std::vector<std::vector<double>>{{1.}, {1.}},
           std::vector<std::vector<double>>{{5000.}, {5000.}});
  info.update(std::vector<unsigned int>{kMinTimeInterval, kMaxTimeInterval});
  casacore::MeasurementSet ms_in("../tNDPPP_tmp.MS");
  MSReader reader(ms_in, parset, prefix);
  MSBDAWriter writer(&reader, msOutName, parset, prefix);

  // Execute
  writer.setInfo(info);

  // Assert if the correct columns are created
  Table table(msOutName, TableLock::AutoNoReadLocking);
  const TableDesc td = table.keywordSet().asTable("BDA_TIME_AXIS").tableDesc();

  // Assert that the data is filled correctly
  Table t = table.keywordSet().asTable("BDA_TIME_AXIS");
  // Assert that the interval is correct
  BOOST_TEST(t.col("MIN_TIME_INTERVAL").getDouble(0) == kMinTimeInterval);
  BOOST_TEST(t.col("MAX_TIME_INTERVAL").getDouble(0) == kMaxTimeInterval);
  BOOST_TEST(t.col("UNIT_TIME_INTERVAL").getDouble(0) == timeInterval);
}

BOOST_AUTO_TEST_SUITE_END()
