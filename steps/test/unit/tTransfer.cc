// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../Transfer.h"

#include <boost/test/unit_test.hpp>

#include <dp3/base/DPInfo.h>

#include "../../../common/ParameterSet.h"

using dp3::steps::Transfer;

BOOST_AUTO_TEST_SUITE(transfer)

BOOST_AUTO_TEST_CASE(constructor) {
  const std::string kSourceMs = "tNDPPP_tmp.MS";
  const std::string kPerformDataTransfer = "true";

  dp3::common::ParameterSet parset;
  parset.add("source_ms", kSourceMs);
  parset.add("data", kPerformDataTransfer);
  BOOST_CHECK_NO_THROW(Transfer(parset, ""));

  const std::string kWrongDataColumn = "TEMP_DATA";
  parset.add("datacolumn", kWrongDataColumn);
  BOOST_CHECK_THROW(Transfer(parset, ""), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
