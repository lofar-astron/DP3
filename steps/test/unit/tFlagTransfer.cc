// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../FlagTransfer.h"

#include <boost/test/unit_test.hpp>

#include <dp3/base/DPInfo.h>

#include "../../../common/ParameterSet.h"

using dp3::steps::FlagTransfer;

BOOST_AUTO_TEST_SUITE(flagtransfer)

BOOST_AUTO_TEST_CASE(constructor) {
  const std::string kSourceMs = "tNDPPP_tmp.MS";

  dp3::common::ParameterSet parset;
  parset.add("source_ms", kSourceMs);

  BOOST_CHECK_NO_THROW(FlagTransfer(parset, ""));
}

BOOST_AUTO_TEST_SUITE_END()