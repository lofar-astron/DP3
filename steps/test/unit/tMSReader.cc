// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MSReader.h"

#include <memory>

#include <boost/test/unit_test.hpp>

#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/TableCopy.h>

#include "../../../common/ParameterSet.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"

using casacore::Table;

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
