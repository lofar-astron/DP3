// tSourceDBUtil.cc: Test program for class SourceDBUtil
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SourceDBUtil.h"
#include "../../PointSource.h"

#include "../../../parmdb/SourceDB.h"

#include <boost/test/unit_test.hpp>

using dp3::base::Direction;
using dp3::base::ModelComponent;
using dp3::base::Patch;
using dp3::base::PointSource;

BOOST_AUTO_TEST_SUITE(source_db_util)

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

BOOST_AUTO_TEST_SUITE_END()
