// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

#include "../../BDABuffer.h"
#include "../../MSBDAWriter.h"
#include "../../MSReader.h"
#include "../../../Common/ParameterSet.h"
#include "./fixtures/fDirectory.cc"

using DP3::ParameterSet;
using DP3::DPPP::BDABuffer;
using DP3::DPPP::DPInfo;
using DP3::DPPP::MSBDAWriter;
using DP3::DPPP::MSReader;

using casacore::Bool;
using casacore::MeasurementSet;
using casacore::TableDesc;
using casacore::TableExprNode;
using casacore::TableLock;

const std::string prefix = "msout.";

BOOST_AUTO_TEST_SUITE(msbdawriter)

BOOST_FIXTURE_TEST_CASE(CreateBDATimeAxis, FixtureDirectory) {
  // Setup test
  const std::string msOutName = "bda_ms_out.MS";
  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  info.init(1, 0, 1, 1, 3.0, 1.5, "", "");
  info.set(std::vector<std::string>{"ant"}, std::vector<double>{1.0},
           {casacore::MVPosition{0, 0, 0}}, std::vector<int>{0},
           std::vector<int>{0});
  MSReader reader("../tNDPPP_tmp.MS", parset, prefix);
  MSBDAWriter writer(&reader, msOutName, parset, prefix);

  // Execute
  writer.setInfo(info);

  // Assert if the correct columns are created
  Table table(msOutName, TableLock::AutoNoReadLocking);
  BOOST_TEST(table.keywordSet().isDefined("BDA_TIME_AXIS"));
  const TableDesc td = table.keywordSet().asTable("BDA_TIME_AXIS").tableDesc();
  BOOST_TEST(td.isColumn("BDA_TIME_AXIS_ID"));
  BOOST_TEST(td.isColumn("IS_BDA_APPLIED"));
  BOOST_TEST(td.isColumn("SINGLE_FACTOR_PER_BASELINE"));
  BOOST_TEST(td.isColumn("MAX_TIME_INTERVAL"));
  BOOST_TEST(td.isColumn("MIN_TIME_INTERVAL"));
  BOOST_TEST(td.isColumn("UNIT_TIME_INTERVAL"));
  BOOST_TEST(td.isColumn("INTEGER_INTERVAL_FACTORS"));
  BOOST_TEST(td.isColumn("HAS_BDA_ORDERING"));

  BOOST_TEST(!td.isColumn("BDA_FREQ_AXIS_ID"));
  BOOST_TEST(!td.isColumn("FIELD_ID"));
}

BOOST_FIXTURE_TEST_CASE(CreateMetaDataFrequencyColumns, FixtureDirectory) {
  // Setup test
  const std::string msOutName = "bda_ms_out_freq.MS";
  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  info.init(1, 0, 1, 1, 3.0, 1.5, "", "");
  info.set(std::vector<std::string>{"ant"}, std::vector<double>{1.0},
           {casacore::MVPosition{0, 0, 0}}, std::vector<int>{0},
           std::vector<int>{0});
  MSReader reader("../tNDPPP_tmp.MS", parset, prefix);
  MSBDAWriter writer(&reader, msOutName, parset, prefix);

  // Execute
  writer.setInfo(info);

  // Assert if the correct columns are created
  Table table(msOutName, TableLock::AutoNoReadLocking);
  BOOST_TEST(table.keywordSet().isDefined("SPECTRAL_WINDOW"));
  const TableDesc td =
      table.keywordSet().asTable("SPECTRAL_WINDOW").tableDesc();
  BOOST_TEST(td.isColumn("BDA_FREQ_AXIS_ID"));
}

BOOST_FIXTURE_TEST_CASE(process_simple, FixtureDirectory) {
  // BOOST_AUTO_TEST_CASE(process_simple) {
  const double kTime(3.0);
  const double kInterval(1.5);
  const std::complex<float> kData(42.0, 43.0);
  const float kWeight(44.0);
  const bool kFlag(false);
  const double kUVW[3]{45.0, 46.0, 47.0};
  const std::string kMsName = "bda_simple.MS";

  // TODO: Remove reader, when the writer can create subtables etc.
  MSReader reader("../tNDPPP_tmp.MS", ParameterSet(), "");
  MSBDAWriter writer(&reader, kMsName, ParameterSet(), "");

  DPInfo info;
  info.init(1, 0, 1, 1, kTime, kInterval, "", "");
  info.set(std::vector<std::string>{"ant"}, std::vector<double>{1.0},
           {casacore::MVPosition{0, 0, 0}}, std::vector<int>{0},
           std::vector<int>{0});
  writer.updateInfo(info);

  auto buffer = boost::make_unique<BDABuffer>(1);
  buffer->AddRow(kTime, kInterval, 0, 0, 1, 1, &kData, &kFlag, &kWeight,
                 nullptr, kUVW);
  writer.process(std::move(buffer));
  writer.finish();

  Table table(kMsName, TableLock::AutoNoReadLocking);
}

BOOST_AUTO_TEST_SUITE_END()
