// tSourceDBUtil.cc: Test program for class SourceDBUtil
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SourceDBUtil.h"
#include "../../PointSource.h"

#include "../../../common/test/unit/fixtures/fDirectory.h"
#include "../../../common/test/unit/fixtures/fSkymodel.h"
#include "../../../parmdb/SourceDB.h"

#include <boost/test/unit_test.hpp>

using dp3::base::Direction;
using dp3::base::ModelComponent;
using dp3::base::Patch;
using dp3::base::PointSource;

static const std::string kSkymodelName = "unittest.skymodel";
static const std::string kSourceDBName = "unittest.sourcedb";

BOOST_AUTO_TEST_SUITE(
    source_db_util,
    // clang-format off
    *boost::unit_test::fixture<FixtureDirectory>()
    *boost::unit_test::fixture<FixtureSkymodel>(FixtureSkymodel::Arguments{
            kSkymodelName, kSourceDBName}))
// clang-format on

BOOST_AUTO_TEST_CASE(cluster_proximate_sources_empty) {
  std::vector<Patch::ConstPtr> input;
  std::vector<Patch::ConstPtr> output =
      dp3::base::clusterProximateSources(input, 1.0);
  BOOST_CHECK(output.empty());
}

BOOST_AUTO_TEST_CASE(cluster_proximate_sources_single) {
  std::vector<std::shared_ptr<const Patch>> input(1);
  std::vector<ModelComponent::ConstPtr> components{
      std::make_shared<PointSource>(Direction(0.2, 0.4)),
      std::make_shared<PointSource>(Direction(0.2, 0.5)),
      std::make_shared<PointSource>(Direction(0.2, 0.6)),
      std::make_shared<PointSource>(Direction(1.0, 1.1))};
  input[0] = std::make_shared<Patch>("a", components.begin(), components.end());

  std::vector<Patch::ConstPtr> output =
      dp3::base::clusterProximateSources(input, 0.5);

  BOOST_REQUIRE_EQUAL(output.size(), 2u);

  // Check patch data
  BOOST_REQUIRE_EQUAL(output[0]->nComponents(), 3u);
  BOOST_CHECK_CLOSE(output[0]->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->direction().dec, 0.5, 1e-5);

  BOOST_REQUIRE_EQUAL(output[1]->nComponents(), 1u);
  BOOST_CHECK_CLOSE(output[1]->direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->direction().dec, 1.1, 1e-5);

  // Check sources in patch
  BOOST_CHECK_CLOSE(output[0]->component(0)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(0)->direction().dec, 0.4, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(1)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(1)->direction().dec, 0.5, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(2)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(2)->direction().dec, 0.6, 1e-5);

  BOOST_CHECK_CLOSE(output[1]->component(0)->direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->component(0)->direction().dec, 1.1, 1e-5);
}

BOOST_AUTO_TEST_CASE(cluster_proximate_sources_multi) {
  std::vector<std::shared_ptr<const Patch>> input(2);
  std::vector<ModelComponent::ConstPtr> componentsA{
      std::make_shared<PointSource>(Direction(0.2, 0.4)),
      std::make_shared<PointSource>(Direction(0.2, 0.6)),
      std::make_shared<PointSource>(Direction(1.0, 1.1))};
  std::vector<ModelComponent::ConstPtr> componentsB{
      std::make_shared<PointSource>(Direction(0.2, 0.5))};
  input[0] =
      std::make_shared<Patch>("a", componentsA.begin(), componentsA.end());
  input[1] =
      std::make_shared<Patch>("b", componentsB.begin(), componentsB.end());

  std::vector<Patch::ConstPtr> output =
      dp3::base::clusterProximateSources(input, 0.5);

  BOOST_REQUIRE_EQUAL(output.size(), 2u);

  // Check patch data
  BOOST_REQUIRE_EQUAL(output[0]->nComponents(), 3u);
  BOOST_CHECK_CLOSE(output[0]->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->direction().dec, 0.5, 1e-5);

  BOOST_REQUIRE_EQUAL(output[1]->nComponents(), 1u);
  BOOST_CHECK_CLOSE(output[1]->direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->direction().dec, 1.1, 1e-5);

  // Check sources in patch
  BOOST_CHECK_CLOSE(output[0]->component(0)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(0)->direction().dec, 0.4, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(1)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(1)->direction().dec, 0.6, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(2)->direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->component(2)->direction().dec, 0.5, 1e-5);

  BOOST_CHECK_CLOSE(output[1]->component(0)->direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->component(0)->direction().dec, 1.1, 1e-5);
}

static void TestPatches(dp3::base::SourceDB&& source_db) {
  const std::vector<dp3::base::Patch::ConstPtr> patches =
      source_db.MakePatchList();
  BOOST_REQUIRE_EQUAL(patches.size(), 3);
  CheckEqual(*patches[0], test_source_db::Expected[0]);
  CheckEqual(*patches[1], test_source_db::Expected[1]);
  CheckEqual(*patches[2], test_source_db::Expected[2]);
}

BOOST_AUTO_TEST_CASE(source_db_make_patch_list) {
  const std::vector<std::string> source_patterns;
  TestPatches(dp3::base::SourceDB{kSkymodelName, source_patterns});
  TestPatches(dp3::base::SourceDB{kSourceDBName, source_patterns});
}

static void TestPatchesExplicit(dp3::base::SourceDB&& source_db) {
  const std::vector<dp3::base::Patch::ConstPtr> patches =
      source_db.MakePatchList();
  BOOST_REQUIRE_EQUAL(patches.size(), 2);
  CheckEqual(*patches[0],
             test_source_db::Patch{"HerA", 0.872665, -0.017453, 3.5, 1});
  CheckEqual(*patches[1],
             test_source_db::Patch{"VirAx", 0.366868, -0.174882, 10, 2});
}

// Using
// https://git.astron.nl/ro/lofar/-/blob/master/CEP/ParmDB/test/tmergesourcedb.in_2
BOOST_AUTO_TEST_CASE(source_db_make_patch_list_explicit_patch_list) {
  const std::string kSkymodel = "explicit.skymodel";
  const std::string kSourceDB = "explicit.sourcedb";
  FixtureSkymodel{FixtureSkymodel::Arguments{
      kSkymodel, kSourceDB, R"(# This is test input for tmergesourcedb
 #

format = name,   type, patch,   Ra,   Dec, I, Q, U, V, SpectralIndex, ReferenceFrequency

HerA, point, , 50deg, -1deg, 3.5
,,VirAx, 21deg, 10deg, 10
VirA1, point, VirAx , 21deg, 10deg, 8
VirA2, point, VirAx , 21.1deg, 10.1deg, 2)"}};

  const std::vector<std::string> source_patterns;

  TestPatchesExplicit(dp3::base::SourceDB{kSkymodel, source_patterns});
  TestPatchesExplicit(dp3::base::SourceDB{kSourceDB, source_patterns});
}

static void TestPolarized(dp3::base::SourceDB&& source_db) {
  BOOST_CHECK_EQUAL(source_db.CheckPolarized(), false);
}

BOOST_AUTO_TEST_CASE(source_db_check_polarised) {
  const std::vector<std::string> source_patterns;
  TestPolarized(dp3::base::SourceDB{kSkymodelName, source_patterns});
  TestPolarized(dp3::base::SourceDB{kSourceDBName, source_patterns});
}

BOOST_AUTO_TEST_CASE(source_db_empty_source_pattern_string) {
  const std::vector<std::string> source_patterns{""};
  BOOST_CHECK_THROW((dp3::base::SourceDB{kSkymodelName, source_patterns}),
                    std::runtime_error);
  BOOST_CHECK_THROW((dp3::base::SourceDB{kSourceDBName, source_patterns}),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
