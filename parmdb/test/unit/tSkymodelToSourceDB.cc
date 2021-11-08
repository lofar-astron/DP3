// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SkymodelToSourceDB.h"

#include "../../SourceDB.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"
#include "../../../common/test/unit/fixtures/fSkymodel.h"
#include "../../../base/SourceDBUtil.h"

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>

static const char* kSkymodelName = "unittest.skymodel";

BOOST_AUTO_TEST_SUITE(skymodeltosourcedb,
                      // clang-format off
                      *boost::unit_test::fixture<FixtureDirectory>()
                      *boost::unit_test::fixture<FixtureSkymodel>(
                              FixtureSkymodel::Arguments{kSkymodelName}))
// clang-format on

BOOST_AUTO_TEST_CASE(make_source_db) {
  std::string format =
      dp3::parmdb::skymodel_to_source_db::ReadFormat("", kSkymodelName);
  if (format.empty()) {
    format = "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
  }

  dp3::parmdb::SourceDB sourceDB =
      dp3::parmdb::skymodel_to_source_db::MakeSourceDb(
          kSkymodelName, "read_skymodel_sourcedb", std::string(), format, "",
          "", false, true, false,
          dp3::parmdb::skymodel_to_source_db::GetSearchInfo("", "", ""));

  std::vector<std::string> patches = sourceDB.getPatches();

  BOOST_REQUIRE_EQUAL(patches.size(), 3u);
  BOOST_REQUIRE_EQUAL(patches[0], "center");
  BOOST_REQUIRE_EQUAL(patches[1], "ra_off");
  BOOST_REQUIRE_EQUAL(patches[2], "radec_off");

  patches = dp3::base::makePatchList(sourceDB, std::vector<std::string>{});
  BOOST_REQUIRE_EQUAL(patches.size(), 3u);
  BOOST_REQUIRE_EQUAL(patches[0], "center");
  BOOST_REQUIRE_EQUAL(patches[1], "ra_off");
  BOOST_REQUIRE_EQUAL(patches[2], "radec_off");

  std::vector<dp3::base::Patch::ConstPtr> foo =
      dp3::base::makePatches(sourceDB, patches, patches.size());

  BOOST_REQUIRE_EQUAL(foo.size(), 3u);
  CheckEqual(*foo[0], test_source_db::Expected[0]);
  CheckEqual(*foo[1], test_source_db::Expected[1]);
  CheckEqual(*foo[2], test_source_db::Expected[2]);
}

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
        dp3::base::MakePatchList(source_db, std::vector<std::string>{});
    BOOST_REQUIRE_EQUAL(patch_list.size(), 3u);
    BOOST_CHECK_EQUAL(patch_list[0], "center");
    BOOST_CHECK_EQUAL(patch_list[1], "ra_off");
    BOOST_CHECK_EQUAL(patch_list[2], "radec_off");

    const std::vector<dp3::base::Patch::ConstPtr> patches =
        dp3::base::MakePatches(source_db, patch_list);
    CheckEqual(*patches[0], test_source_db::Expected[0]);
    CheckEqual(*patches[1], test_source_db::Expected[1]);
    CheckEqual(*patches[2], test_source_db::Expected[2]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
