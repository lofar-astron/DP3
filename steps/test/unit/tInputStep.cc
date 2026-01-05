// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include "steps/InputStep.h"

#include <boost/test/unit_test.hpp>

#include "steps/MSBDAReader.h"
#include "steps/MsReader.h"
#include "steps/MultiMsReader.h"
#include "common/ParameterSet.h"
#include "common/test/unit/fixtures/fDirectory.h"

using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::MSBDAReader;
using dp3::steps::MsReader;
using dp3::steps::MultiMsReader;

BOOST_FIXTURE_TEST_SUITE(inputstep, dp3::common::test::FixtureDirectory)

BOOST_AUTO_TEST_CASE(reader_initialization_regular) {
  ExtractResource("tNDPPP.in_MS.tgz");
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MsReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_regular) {
  ExtractResource("tNDPPP.in_MS.tgz");
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MultiMsReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_bda) {
  ExtractResource("tNDPPP_bda.in_MS.tgz");
  ParameterSet parset;
  parset.add("msin", "tNDPPP_bda_tmp.MS");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MSBDAReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_bda) {
  ExtractResource("tNDPPP_bda.in_MS.tgz");
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_bda_tmp.MS, tNDPPP_bda_tmp.MS]");
  BOOST_CHECK_THROW(InputStep::CreateReader(parset), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_one_missing) {
  ExtractResource("tNDPPP.in_MS.tgz");
  ParameterSet parset;
  parset.add("msin", "[missing.ms, tNDPPP_tmp.MS]");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MultiMsReader*>(reader.get()));
}

BOOST_AUTO_TEST_SUITE_END()
