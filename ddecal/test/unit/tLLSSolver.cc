// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../linear_solvers/LLSSolver.h"

#include <boost/test/unit_test.hpp>

using dp3::ddecal::LLSSolver;
using dp3::ddecal::LLSSolverType;

BOOST_AUTO_TEST_SUITE(lls_solver)

BOOST_AUTO_TEST_CASE(parsetype) {
  BOOST_CHECK(LLSSolver::ParseType("svd") == LLSSolverType::SVD);
  BOOST_CHECK(LLSSolver::ParseType("QR") == LLSSolverType::QR);
  BOOST_CHECK(LLSSolver::ParseType("nOrMaLeQuAtIoNs") ==
              LLSSolverType::NORMAL_EQUATIONS);
  BOOST_CHECK_THROW(LLSSolver::ParseType(""), std::runtime_error);
  BOOST_CHECK_THROW(LLSSolver::ParseType("foo"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
