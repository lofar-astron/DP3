#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <boost/test/unit_test.hpp>

#include "../../MSBDAWriter.h"
#include "../../MSReader.h"
#include "../../../Common/ParameterSet.h"

using DP3::ParameterSet;
using DP3::DPPP::DPInfo;
using DP3::DPPP::MSBDAWriter;
using DP3::DPPP::MSReader;

using casacore::Bool;
using casacore::TableLock;
using casacore::MeasurementSet;
using casacore::TableDesc;

BOOST_AUTO_TEST_SUITE(msbdawriter)

BOOST_AUTO_TEST_CASE(createBDATimeAxis) {
  // Setup test
  std::string prefix =  "msout.";
  std::string msOutName = "bda_ms_out.MS";

  ParameterSet parset;
  parset.add(prefix + "overwrite", "true");
  DPInfo info;
  MSReader reader("tNDPPP_tmp.MS", parset,prefix);
  MSBDAWriter writer(&reader, msOutName, parset, prefix);

  // Execute
  writer.setInfo(info);

  // Assert if table exists and no error is thrown
  MeasurementSet ms(msOutName, TableLock::AutoNoReadLocking);
  BOOST_TEST(Bool(true) == ms.keywordSet().isDefined("BDA_TIME_AXIS"));
  BOOST_REQUIRE_NO_THROW(ms.keywordSet().asTable("BDA_TIME_AXIS"));
}

BOOST_AUTO_TEST_SUITE_END()
