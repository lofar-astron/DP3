// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SkymodelToSourceDB.h"
#include "../../SourceDB.h"

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>

BOOST_AUTO_TEST_SUITE(skymodeltosourcedb)

BOOST_AUTO_TEST_CASE(read_skymodel) {
  std::string skymodel_filename = "unittest.skymodel";
  std::ofstream skymodel_file(skymodel_filename);

  skymodel_file
      << "FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, "
         "PositionAngle, ReferenceFrequency='134e6', SpectralIndex='[0.0]'\n";
  skymodel_file
      << "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \n";
  skymodel_file
      << "ra_off, POINT, 16:58:28.205000, +63.44.34.314000, 1, , , , , \n";
  skymodel_file
      << "radec_off, POINT, 16:38:28.205000, +65.44.34.314000, 1, , , , , \n";

  skymodel_file.close();

  std::string format =
      dp3::parmdb::skymodel_to_source_db::ReadFormat("", skymodel_filename);
  if (format.empty()) {
    format = "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
  }

  dp3::parmdb::SourceDB sourceDB =
      dp3::parmdb::skymodel_to_source_db::MakeSourceDb(
          skymodel_filename, "", std::string(), format, "", "", false, true,
          false, dp3::parmdb::skymodel_to_source_db::GetSearchInfo("", "", ""),
          true);

  std::vector<std::string> patches = sourceDB.getPatches();

  BOOST_REQUIRE_EQUAL(patches.size(), 3u);
  BOOST_REQUIRE_EQUAL(patches[0], "center");
  BOOST_REQUIRE_EQUAL(patches[1], "ra_off");
  BOOST_REQUIRE_EQUAL(patches[2], "radec_off");

  remove(skymodel_filename.c_str());
}

BOOST_AUTO_TEST_CASE(read_format) {
  std::string skymodel_filename = "unittest.skymodel";
  std::ofstream skymodel_file(skymodel_filename);

  std::string default_format =
      "Name,Type,Ra,Dec,I,Q,U,V,MajorAxis,MinorAxis,Orientation";
  skymodel_file
      << "center, POINT, 16:38:28.205000, +63.44.34.314000, 1, , , , , \n";

  skymodel_file.close();

  std::string format =
      dp3::parmdb::skymodel_to_source_db::ReadFormat("", skymodel_filename);

  BOOST_REQUIRE_EQUAL(format, default_format);

  remove(skymodel_filename.c_str());
}

BOOST_AUTO_TEST_SUITE_END()