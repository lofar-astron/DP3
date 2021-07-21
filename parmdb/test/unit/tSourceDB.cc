// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../ParmDBMeta.h"
#include "../../SourceDB.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(sourcedb)

BOOST_AUTO_TEST_CASE(open_non_existing) {
  dp3::parmdb::ParmDBMeta parmDBMeta("", "not/an/existing/file");
  BOOST_CHECK_THROW(dp3::parmdb::SourceDB(parmDBMeta, true, false),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
