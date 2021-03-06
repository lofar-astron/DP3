// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

#include "../../MSBDAReader.h"
#include "../../MSReader.h"
#include "../../MultiMSReader.h"
#include "../../InputStep.h"
#include "../../../common/ParameterSet.h"

using dp3::common::ParameterSet;
using dp3::steps::InputStep;
using dp3::steps::MSBDAReader;
using dp3::steps::MSReader;
using dp3::steps::MultiMSReader;

BOOST_AUTO_TEST_SUITE(dpinput)

BOOST_AUTO_TEST_CASE(reader_initialization_regular) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MSReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_regular) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MultiMSReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_bda) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  parset.add("bda", "true");
  std::unique_ptr<InputStep> reader = InputStep::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MSBDAReader*>(reader.get()));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_bda) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  parset.add("bda", "true");

  BOOST_CHECK_THROW(InputStep::CreateReader(parset, ""), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
