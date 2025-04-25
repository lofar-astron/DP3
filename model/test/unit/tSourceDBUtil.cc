// tSourceDBUtil.cc: Test program for class SourceDBUtil
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SourceDBUtil.h"

#include "../../../base/PointSource.h"
#include "../../../common/test/unit/fixtures/fDirectory.h"
#include "../../../common/test/unit/fixtures/fSkymodel.h"
#include "../../../parmdb/SourceDB.h"

#include "../../../base/test/LoggerFixture.h"

#include <boost/test/unit_test.hpp>

using dp3::base::Direction;
using dp3::base::ModelComponent;
using dp3::base::PointSource;

using dp3::model::Patch;

static const std::string kSkymodelName = "unittest.skymodel";
static const std::string kSourceDBName = "unittest.sourcedb";

BOOST_AUTO_TEST_SUITE(
    source_db_util,
    // clang-format off
    *boost::unit_test::fixture<dp3::base::test::LoggerFixture>()
    *boost::unit_test::fixture<dp3::common::test::FixtureDirectory>()
    *boost::unit_test::fixture<FixtureSkymodel>(FixtureSkymodel::Arguments{
            kSkymodelName, kSourceDBName}))
// clang-format on

BOOST_AUTO_TEST_CASE(cluster_proximate_sources_empty) {
  std::vector<std::shared_ptr<Patch>> input;
  std::vector<std::shared_ptr<Patch>> output =
      dp3::model::clusterProximateSources(input, 1.0);
  BOOST_CHECK(output.empty());
}

BOOST_AUTO_TEST_CASE(cluster_proximate_sources_single) {
  std::vector<std::shared_ptr<Patch>> input(1);
  std::vector<std::shared_ptr<ModelComponent>> components{
      std::make_shared<PointSource>(Direction(0.2, 0.4)),
      std::make_shared<PointSource>(Direction(0.2, 0.5)),
      std::make_shared<PointSource>(Direction(0.2, 0.6)),
      std::make_shared<PointSource>(Direction(1.0, 1.1))};
  input[0] = std::make_shared<Patch>("a", components.begin(), components.end());

  std::vector<std::shared_ptr<Patch>> output =
      dp3::model::clusterProximateSources(input, 0.5);

  BOOST_REQUIRE_EQUAL(output.size(), 2u);

  // Check patch data
  BOOST_REQUIRE_EQUAL(output[0]->NComponents(), 3u);
  BOOST_CHECK_CLOSE(output[0]->Direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->Direction().dec, 0.5, 1e-5);

  BOOST_REQUIRE_EQUAL(output[1]->NComponents(), 1u);
  BOOST_CHECK_CLOSE(output[1]->Direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->Direction().dec, 1.1, 1e-5);

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
  std::vector<std::shared_ptr<Patch>> input(2);
  std::vector<std::shared_ptr<ModelComponent>> componentsA{
      std::make_shared<PointSource>(Direction(0.2, 0.4)),
      std::make_shared<PointSource>(Direction(0.2, 0.6)),
      std::make_shared<PointSource>(Direction(1.0, 1.1))};
  std::vector<std::shared_ptr<ModelComponent>> componentsB{
      std::make_shared<PointSource>(Direction(0.2, 0.5))};
  input[0] =
      std::make_shared<Patch>("a", componentsA.begin(), componentsA.end());
  input[1] =
      std::make_shared<Patch>("b", componentsB.begin(), componentsB.end());

  std::vector<std::shared_ptr<Patch>> output =
      dp3::model::clusterProximateSources(input, 0.5);

  BOOST_REQUIRE_EQUAL(output.size(), 2u);

  // Check patch data
  BOOST_REQUIRE_EQUAL(output[0]->NComponents(), 3u);
  BOOST_CHECK_CLOSE(output[0]->Direction().ra, 0.2, 1e-5);
  BOOST_CHECK_CLOSE(output[0]->Direction().dec, 0.5, 1e-5);

  BOOST_REQUIRE_EQUAL(output[1]->NComponents(), 1u);
  BOOST_CHECK_CLOSE(output[1]->Direction().ra, 1.0, 1e-5);
  BOOST_CHECK_CLOSE(output[1]->Direction().dec, 1.1, 1e-5);

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

static void TestPatches(dp3::model::SourceDBWrapper& source_db,
                        const std::vector<test_source_db::Patch>& expected) {
  const std::vector<std::shared_ptr<Patch>> patches = source_db.MakePatchList();
  BOOST_REQUIRE_EQUAL(patches.size(), expected.size());
  for (size_t i = 0; i < patches.size(); ++i)
    test_source_db::CheckEqual(*patches[i], expected[i]);
}

BOOST_AUTO_TEST_CASE(source_db_make_patch_list_empty_pattern) {
  const std::vector<std::string> filter;
  const std::vector<test_source_db::Patch> expected{
      test_source_db::Expected.begin(), test_source_db::Expected.end()};

  TestPatches(
      dp3::model::SourceDBWrapper(kSkymodelName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern),
      expected);
  TestPatches(
      dp3::model::SourceDBWrapper(kSourceDBName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern),
      expected);
}

BOOST_AUTO_TEST_CASE(source_db_make_patch_list_non_empty_pattern) {
  // The pattern is in non-alphabetic order, the results are expected to be
  // sorted.
  const std::vector<std::string> filter{"radec_off", "center"};
  const std::vector<test_source_db::Patch> expected{
      test_source_db::Expected[0], test_source_db::Expected[2]};

  TestPatches(
      dp3::model::SourceDBWrapper(kSkymodelName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern),
      expected);
  TestPatches(
      dp3::model::SourceDBWrapper(kSourceDBName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern),
      expected);
}

BOOST_AUTO_TEST_CASE(source_db_make_patch_list_empty_values) {
  const std::vector<std::string> filter;
  const std::vector<test_source_db::Patch> expected;
  TestPatches(
      dp3::model::SourceDBWrapper(kSkymodelName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kValue),
      expected);
  TestPatches(
      dp3::model::SourceDBWrapper(kSourceDBName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kValue),
      expected);
}

BOOST_AUTO_TEST_CASE(source_db_make_patch_list_non_empty_values) {
  // The pattern is in non-alphabetic order, the results are expected to be
  // in the same order.
  const std::vector<std::string> filter{"radec_off", "center"};
  const std::vector<test_source_db::Patch> expected{
      test_source_db::Expected[2], test_source_db::Expected[0]};

  TestPatches(
      dp3::model::SourceDBWrapper(kSkymodelName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kValue),
      expected);
  TestPatches(
      dp3::model::SourceDBWrapper(kSourceDBName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kValue),
      expected);
}

static void TestPatchesExplicit(dp3::model::SourceDBWrapper& source_db) {
  const std::vector<std::shared_ptr<Patch>> patches = source_db.MakePatchList();
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

  const std::vector<std::string> filter;

  TestPatchesExplicit(dp3::model::SourceDBWrapper(kSkymodel).Filter(
      filter, dp3::model::SourceDBWrapper::FilterMode::kPattern));
  TestPatchesExplicit(dp3::model::SourceDBWrapper(kSourceDB).Filter(
      filter, dp3::model::SourceDBWrapper::FilterMode::kPattern));
}

static void TestPolarized(dp3::model::SourceDBWrapper& source_db) {
  BOOST_CHECK_EQUAL(source_db.CheckPolarized(), false);
}

BOOST_AUTO_TEST_CASE(source_db_check_polarised) {
  const std::vector<std::string> filter;
  TestPolarized(
      dp3::model::SourceDBWrapper(kSkymodelName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern));
  TestPolarized(
      dp3::model::SourceDBWrapper(kSourceDBName)
          .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern));
}

BOOST_AUTO_TEST_CASE(source_db_empty_source_pattern_string) {
  const std::vector<std::string> filter{""};
  BOOST_CHECK_THROW(
      (dp3::model::SourceDBWrapper(kSkymodelName)
           .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern)),
      std::runtime_error);
  BOOST_CHECK_THROW(
      (dp3::model::SourceDBWrapper(kSourceDBName)
           .Filter(filter, dp3::model::SourceDBWrapper::FilterMode::kPattern)),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
