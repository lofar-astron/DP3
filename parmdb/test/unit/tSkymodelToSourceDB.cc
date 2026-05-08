// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "parmdb/SkymodelToSourceDB.h"

#include "parmdb/SourceDB.h"
#include "base/test/LoggerFixture.h"
#include "common/test/unit/fixtures/fDirectory.h"
#include "common/test/unit/fixtures/fSkymodel.h"
#include "parmdb/SkymodelToSourceDB.h"
#include "model/SourceDBUtil.h"

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>

static const char* kSkymodelName = "unittest.skymodel";

BOOST_AUTO_TEST_SUITE(
    skymodeltosourcedb,
    *boost::unit_test::fixture<dp3::base::test::LoggerFixture>() *
        boost::unit_test::fixture<dp3::common::test::FixtureDirectory>() *
        boost::unit_test::fixture<FixtureSkymodel>(FixtureSkymodel::Arguments{
            kSkymodelName}))

BOOST_AUTO_TEST_CASE(read_format) {
  std::string skymodel_filename = "read_format.skymodel";
  std::ofstream skymodel_file(skymodel_filename);

  std::string default_format =
      "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
  skymodel_file
      << "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \n";

  skymodel_file.close();

  std::string format =
      dp3::parmdb::skymodel_to_source_db::ReadFormat("", skymodel_filename);

  BOOST_REQUIRE_EQUAL(format, default_format);
}

BOOST_AUTO_TEST_CASE(read_skymodel) {
  const std::string format =
      dp3::parmdb::skymodel_to_source_db::ReadFormat("", kSkymodelName);

  BOOST_REQUIRE(!format.empty());
  const dp3::parmdb::SourceDBSkymodel source_db =
      dp3::parmdb::skymodel_to_source_db::MakeSourceDBSkymodel(kSkymodelName,
                                                               format);

  {
    const std::vector<dp3::parmdb::PatchInfo>& patches = source_db.GetPatches();

    BOOST_REQUIRE_EQUAL(patches.size(), 3u);
    BOOST_CHECK_EQUAL(patches[0].getName(), "center");
    BOOST_CHECK_EQUAL(patches[1].getName(), "ra_off");
    BOOST_CHECK_EQUAL(patches[2].getName(), "radec_off");
  }

  {
    const std::vector<std::string> patch_list =
        dp3::model::MakePatchList(source_db, std::vector<std::string>{});
    BOOST_REQUIRE_EQUAL(patch_list.size(), 3u);
    BOOST_CHECK_EQUAL(patch_list[0], "center");
    BOOST_CHECK_EQUAL(patch_list[1], "ra_off");
    BOOST_CHECK_EQUAL(patch_list[2], "radec_off");

    const std::vector<std::shared_ptr<dp3::model::Patch>> patches =
        dp3::model::MakePatches(source_db, patch_list);
    CheckEqual(*patches[0], test_source_db::Expected[0]);
    CheckEqual(*patches[1], test_source_db::Expected[1]);
    CheckEqual(*patches[2], test_source_db::Expected[2]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
