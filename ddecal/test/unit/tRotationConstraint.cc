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

BOOST_AUTO_TEST_SUITE(rotationconstraint)

BOOST_AUTO_TEST_CASE(test_rotation_changes_input) {
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

BOOST_AUTO_TEST_CASE(test_rotation) {
  RotationConstraint constraint;
  constraint.Initialize(1, {1u}, {42.0e6});

  dp3::ddecal::SolutionTensor onesolution_tensor({1, 1, 1, 4});
  dp3::ddecal::SolutionSpan onesolution =
      aocommon::xt::CreateSpan(onesolution_tensor);
  double pi = 3.1415;
  for (double phi = -pi; phi < pi; phi += pi / 6) {
    /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
    */
    onesolution(0, 0, 0, 0) = cos(phi);
    onesolution(0, 0, 0, 1) = -sin(phi);
    onesolution(0, 0, 0, 2) = sin(phi);
    onesolution(0, 0, 0, 3) = cos(phi);
    vector<Constraint::Result> constraint_result;
    constraint_result = constraint.Apply(onesolution, 0., nullptr);

    BOOST_CHECK(constraint_result.size() == 1);
    BOOST_CHECK(constraint_result[0].axes == "ant,dir,freq");
    BOOST_CHECK(near(constraint_result[0].vals[0], phi));
    BOOST_CHECK(constraint_result[0].name == "rotation");
    BOOST_CHECK(constraint_result[0].dims.size() == 3);
    BOOST_CHECK(constraint_result[0].dims[0] == 1);
    BOOST_CHECK(constraint_result[0].dims[1] == 1);
    BOOST_CHECK(constraint_result[0].dims[2] == 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()
