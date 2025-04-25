// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <casacore/casa/BasicMath/Math.h>  // near

#include <vector>
#include <complex>

#include "../../constraints/RotationConstraint.h"
#include "../../constraints/RotationAndDiagonalConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using casacore::near;
using dp3::ddecal::Constraint;
using dp3::ddecal::RotationAndDiagonalConstraint;
using dp3::ddecal::RotationConstraint;
using std::complex;
using std::vector;

BOOST_AUTO_TEST_SUITE(rotation_constraint)

BOOST_AUTO_TEST_CASE(rotation_changes_input) {
  RotationConstraint constraint;
  constraint.Initialize(1, {1u}, {42.0e6});

  dp3::ddecal::SolutionTensor onesolution_tensor({1, 1, 1, 4});
  dp3::ddecal::SolutionSpan onesolution =
      aocommon::xt::CreateSpan(onesolution_tensor);

  onesolution(0, 0, 0, 0) = 1.0;
  onesolution(0, 0, 0, 1) = 2.0;
  onesolution(0, 0, 0, 2) = 3.0;
  onesolution(0, 0, 0, 3) = 4.0;

  constraint.Apply(onesolution, 0.0, nullptr);

  BOOST_CHECK_CLOSE(std::real(onesolution(0, 0, 0, 0)), 0.980581, 1.0e-3);
  BOOST_CHECK_CLOSE(std::real(onesolution(0, 0, 0, 1)), -0.196116, 1.0e-3);
  BOOST_CHECK_CLOSE(std::real(onesolution(0, 0, 0, 2)), 0.196116, 1.0e-3);
  BOOST_CHECK_CLOSE(std::real(onesolution(0, 0, 0, 3)), 0.980581, 1.0e-3);
}

BOOST_AUTO_TEST_CASE(rotation) {
  constexpr size_t kNDirections = 3;
  RotationConstraint constraint;
  constraint.Initialize(1, std::vector<uint32_t>(kNDirections, 1u), {42.0e6});

  dp3::ddecal::SolutionTensor onesolution_tensor({1, 1, kNDirections, 4});
  dp3::ddecal::SolutionSpan onesolution =
      aocommon::xt::CreateSpan(onesolution_tensor);
  for (double phi = -M_PI; phi + 0.01 < M_PI; phi += M_PI / 6) {
    for (size_t direction = 0; direction != kNDirections; ++direction) {
      onesolution(0, 0, direction, 0) = std::cos(phi + 0.3 * direction);
      onesolution(0, 0, direction, 1) = -std::sin(phi + 0.3 * direction);
      onesolution(0, 0, direction, 2) = std::sin(phi + 0.3 * direction);
      onesolution(0, 0, direction, 3) = std::cos(phi + 0.3 * direction);
    }
    vector<Constraint::Result> constraint_result =
        constraint.Apply(onesolution, 0.0, nullptr);

    BOOST_REQUIRE_EQUAL(constraint_result.size(), 1);
    BOOST_CHECK_EQUAL(constraint_result[0].axes, "ant,dir,freq");
    BOOST_CHECK_EQUAL(constraint_result[0].name, "rotation");
    BOOST_REQUIRE_EQUAL(constraint_result[0].dims.size(), 3);
    BOOST_CHECK_EQUAL(constraint_result[0].dims[0], 1);
    BOOST_CHECK_EQUAL(constraint_result[0].dims[1], kNDirections);
    BOOST_CHECK_EQUAL(constraint_result[0].dims[2], 1);
    const std::vector<double> &values = constraint_result[0].vals;
    BOOST_REQUIRE_EQUAL(values.size(), kNDirections);
    for (size_t direction = 0; direction != kNDirections; ++direction) {
      double expected_phi = phi + 0.3 * direction;
      if (expected_phi + 0.001 >= M_PI) expected_phi -= 2.0 * M_PI;
      BOOST_CHECK_CLOSE_FRACTION(values[direction], expected_phi, 1e-4);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
