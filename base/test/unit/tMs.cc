// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../MS.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ms)

BOOST_AUTO_TEST_CASE(read_antenna_set) {
  // This MS has an antenna set field containing "HBA_DUAL_INNER".
  const std::string ms_name = "tNDPPP-generic.MS";
  const casacore::MeasurementSet ms(ms_name);
  const std::string antenna_set = dp3::base::ReadAntennaSet(ms);
  BOOST_TEST(antenna_set == "HBA_DUAL_INNER");
}

BOOST_AUTO_TEST_CASE(read_antenna_set_absent) {
  // This MS has no antenna set field.
  const std::string ms_name = "tNDPPP_tmp.MS";
  const casacore::MeasurementSet ms(ms_name);
  const std::string antenna_set = dp3::base::ReadAntennaSet(ms);
  BOOST_TEST(antenna_set.empty());
}

BOOST_AUTO_TEST_SUITE_END()
