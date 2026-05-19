#include "base/PointSource.h"
#include "base/test/LoggerFixture.h"

#include "common/test/unit/fixtures/fDirectory.h"
#include "common/test/unit/fixtures/fSkyModel.h"

#include "sky_model/ReadSkyModel.h"
#include "sky_model/SkyModelSelection.h"

#include <boost/test/unit_test.hpp>

using dp3::base::Direction;
using dp3::base::ModelComponent;
using dp3::base::PointSource;

using dp3::sky_model::Patch;

static const std::string kSkyModelName = "unittest.skymodel";

static void TestPatches(dp3::sky_model::SkyModelSelection& source_db,
                        const std::vector<test_source_db::Patch>& expected) {
  const std::vector<std::shared_ptr<Patch>> patches = source_db.MakePatchList();
  BOOST_REQUIRE_EQUAL(patches.size(), expected.size());
  for (size_t i = 0; i < patches.size(); ++i)
    test_source_db::CheckEqual(*patches[i], expected[i]);
}

BOOST_AUTO_TEST_SUITE(
    sky_model_selection,
    *boost::unit_test::fixture<dp3::base::test::LoggerFixture>() *
        boost::unit_test::fixture<dp3::common::test::FixtureDirectory>() *
        boost::unit_test::fixture<FixtureSkyModel>(FixtureSkyModel::Arguments{
            kSkyModelName}))

BOOST_AUTO_TEST_CASE(select_matching_patches_with_empty_pattern) {
  const std::vector<std::string> filter;
  const std::vector<test_source_db::Patch> expected(
      test_source_db::Expected.begin(), test_source_db::Expected.end());

  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectMatchingPatches(filter);
  TestPatches(selection, expected);
}

BOOST_AUTO_TEST_CASE(select_all_patches) {
  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  const std::vector<test_source_db::Patch> expected(
      test_source_db::Expected.begin(), test_source_db::Expected.end());
  selection.SelectAllPatches();
  TestPatches(selection, expected);
}

BOOST_AUTO_TEST_CASE(select_matching_patches) {
  // The pattern is in non-alphabetic order, the results are expected to be
  // sorted.
  const std::vector<std::string> filter{"radec_off", "center"};
  const std::vector<test_source_db::Patch> expected{
      test_source_db::Expected[0], test_source_db::Expected[2]};

  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectMatchingPatches(filter);
  TestPatches(selection, expected);
}

BOOST_AUTO_TEST_CASE(select_patch_list_with_empty_filter) {
  const std::vector<std::string> filter;
  const std::vector<test_source_db::Patch> expected;
  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectPatchList(filter);
  TestPatches(selection, expected);
}

BOOST_AUTO_TEST_CASE(select_patch_list) {
  // The pattern is in non-alphabetic order, the results are expected to be
  // in the same order.
  const std::vector<std::string> filter{"radec_off", "center"};
  const std::vector<test_source_db::Patch> expected{
      test_source_db::Expected[2], test_source_db::Expected[0]};

  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectPatchList(filter);
  TestPatches(selection, expected);
}

static void TestPatchesExplicit(dp3::sky_model::SkyModelSelection& source_db) {
  const std::vector<std::shared_ptr<Patch>> patches = source_db.MakePatchList();
  BOOST_REQUIRE_EQUAL(patches.size(), 2);
  CheckEqual(*patches[0],
             test_source_db::Patch{"HerA", 0.872665, -0.017453, 3.5, 1});
  CheckEqual(*patches[1],
             test_source_db::Patch{"VirAx", 0.366868, -0.174882, 10, 2});
}

BOOST_AUTO_TEST_CASE(select_matching_patches_b) {
  const std::string kSkyModel = "explicit.skymodel";
  FixtureSkyModel{FixtureSkyModel::Arguments{
      kSkyModel, R"(# This is test input for tmergesourcedb
 #

format = name,   type, patch,   Ra,   Dec, I, Q, U, V, SpectralIndex, ReferenceFrequency

HerA, point, , 50deg, -1deg, 3.5
,,VirAx, 21deg, 10deg, 10
VirA1, point, VirAx , 21deg, 10deg, 8
VirA2, point, VirAx , 21.1deg, 10.1deg, 2)"}};

  const std::vector<std::string> filter;

  dp3::sky_model::SkyModel sky_model = dp3::sky_model::ReadSkyModel(kSkyModel);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectMatchingPatches(filter);
  TestPatchesExplicit(selection);
}

static void TestPolarized(dp3::sky_model::SkyModelSelection& source_db) {
  BOOST_CHECK_EQUAL(source_db.CheckPolarized(), false);
}

BOOST_AUTO_TEST_CASE(check_polarised) {
  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  selection.SelectMatchingPatches({});
  TestPolarized(selection);
}

BOOST_AUTO_TEST_CASE(empty_source_pattern_string) {
  dp3::sky_model::SkyModel sky_model =
      dp3::sky_model::ReadSkyModel(kSkyModelName);
  dp3::sky_model::SkyModelSelection selection(sky_model);
  BOOST_CHECK_THROW(selection.SelectMatchingPatches({""}), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
