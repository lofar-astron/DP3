#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <boost/test/unit_test.hpp>

#include "../../MSBDAWriter.h"
#include "../../MSReader.h"
#include "../../../Common/ParameterSet.h"

using DP3::ParameterSet;
using DP3::DPPP::DPInfo;
using DP3::DPPP::MSBDAWriter;
using DP3::DPPP::MSReader;

using casacore::Bool;
using casacore::MeasurementSet;
using casacore::TableDesc;
using casacore::TableExprNode;
using casacore::TableLock;

BOOST_AUTO_TEST_SUITE(msbdawriter)

BOOST_AUTO_TEST_CASE(CreateBDATimeAxis) {
  // Setup test
  std::string prefix = "msout.";
  std::string msOutName = "bda_ms_out.MS";

  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  MSReader reader("tNDPPP_tmp.MS", parset, prefix);
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

BOOST_AUTO_TEST_CASE(AddMetaDataFrequency) {
  // Setup test
  std::string prefix = "msout.";
  std::string msOutName = "bda_ms_out_frequency.MS";

  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  MSReader reader("tNDPPP_tmp.MS", parset, prefix);
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

BOOST_AUTO_TEST_SUITE_END()
