// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

#include "../../MSBDAReader.h"
#include "../../MSReader.h"
#include "../../MultiMSReader.h"
#include "../../DPInput.h"
#include "../../../Common/ParameterSet.h"

using DP3::ParameterSet;
using DP3::DPPP::DPInput;
using DP3::DPPP::MSBDAReader;
using DP3::DPPP::MSReader;
using DP3::DPPP::MultiMSReader;

BOOST_AUTO_TEST_SUITE(dpinput)

BOOST_AUTO_TEST_CASE(reader_initialization_regular) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  DPInput* reader = DPInput::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MSReader*>(reader));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_regular) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  DPInput* reader = DPInput::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MultiMSReader*>(reader));
}

BOOST_AUTO_TEST_CASE(reader_initialization_bda) {
  ParameterSet parset;
  parset.add("msin", "tNDPPP_tmp.MS");
  parset.add("bda", "true");
  DPInput* reader = DPInput::CreateReader(parset, "");

  BOOST_TEST(dynamic_cast<MSBDAReader*>(reader));
}

BOOST_AUTO_TEST_CASE(reader_initialization_multiple_bda) {
  ParameterSet parset;
  parset.add("msin", "[tNDPPP_tmp.MS, tNDPPP_tmp.MS]");
  parset.add("bda", "true");

  BOOST_CHECK_THROW(DPInput::CreateReader(parset, ""), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
