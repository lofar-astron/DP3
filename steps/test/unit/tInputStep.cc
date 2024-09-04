// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <boost/test/unit_test.hpp>

#include "../../MSBDAReader.h"
#include "../../MsReader.h"
#include "../../MultiMsReader.h"
#include "../../InputStep.h"
#include "../../../common/ParameterSet.h"

using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::MSBDAReader;
using dp3::steps::MsReader;
using dp3::steps::MultiMsReader;

BOOST_AUTO_TEST_SUITE(inputstep)

BOOST_AUTO_TEST_CASE(reader_initialization_regular) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MsReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_regular) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MultiMsReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_bda) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_bda_tmp.MS");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MSBDAReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_bda) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_bda_tmp.MS, tNDPPP_bda_tmp.MS]");
  BOOST_CHECK_THROW(InputStep::CreateReader(parset), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_one_missing) {
  ParameterSet parset;
  parset.add("msin", "[missing.ms, tNDPPP_tmp.MS]");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset);

  BOOST_TEST(dynamic_cast<MultiMsReader*>(reader.get()));
}

BOOST_AUTO_TEST_SUITE_END()
