// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "sky_model/ReadSkyModel.h"
#include "sky_model/SkyModelFunctions.h"

#include "base/test/LoggerFixture.h"
#include "common/test/unit/fixtures/fDirectory.h"
#include "common/test/unit/fixtures/fSkyModel.h"

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>

static const char* kSkyModelName = "unittest.skymodel";

BOOST_AUTO_TEST_SUITE(
    sky_model_functions,
    *boost::unit_test::fixture<dp3::base::test::LoggerFixture>() *
        boost::unit_test::fixture<dp3::common::test::FixtureDirectory>() *
        boost::unit_test::fixture<FixtureSkyModel>(FixtureSkyModel::Arguments{
            kSkyModelName}))

BOOST_AUTO_TEST_CASE(read_format) {
  std::string sky_model_filename = "read_format.skymodel";
  std::ofstream sky_model_file(sky_model_filename);

  std::string default_format =
      "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
  sky_model_file
      << "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \n";

  sky_model_file.close();

  std::string format = dp3::sky_model::ReadFormat("", sky_model_filename);

  BOOST_REQUIRE_EQUAL(format, default_format);
}

BOOST_AUTO_TEST_CASE(read_sky_model) {
  const std::string format = dp3::sky_model::ReadFormat("", kSkyModelName);

  BOOST_REQUIRE(!format.empty());
  const dp3::sky_model::SkyModel source_db =
      dp3::sky_model::ReadSkyModel(kSkyModelName, format);

  {
    const std::vector<dp3::sky_model::PatchInfo>& patches =
        source_db.GetPatches();

    BOOST_REQUIRE_EQUAL(patches.size(), 3u);
    BOOST_CHECK_EQUAL(patches[0].getName(), "center");
    BOOST_CHECK_EQUAL(patches[1].getName(), "ra_off");
    BOOST_CHECK_EQUAL(patches[2].getName(), "radec_off");
  }

  {
    const std::vector<std::string> patch_list =
        dp3::sky_model::MakePatchList(source_db, std::vector<std::string>{});
    BOOST_REQUIRE_EQUAL(patch_list.size(), 3u);
    BOOST_CHECK_EQUAL(patch_list[0], "center");
    BOOST_CHECK_EQUAL(patch_list[1], "ra_off");
    BOOST_CHECK_EQUAL(patch_list[2], "radec_off");

    const std::vector<std::shared_ptr<dp3::sky_model::Patch>> patches =
        dp3::sky_model::MakePatches(source_db, patch_list);
    CheckEqual(*patches[0], test_source_db::Expected[0]);
    CheckEqual(*patches[1], test_source_db::Expected[1]);
    CheckEqual(*patches[2], test_source_db::Expected[2]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
