// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <casacore/casa/BasicMath/Math.h>  // near

#include <vector>
#include <complex>

#include "../../constraints/RotationAndDiagonalConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace dp3::ddecal {

BOOST_AUTO_TEST_SUITE(rotation_and_diagonal_constraint)

BOOST_DATA_TEST_CASE(test_rotation_and_diagonal,
                     boost::unit_test::data::make({true, false}),
                     doRotationReference) {
  RotationAndDiagonalConstraint constraint;
  constraint.Initialize(1, {1u}, {42e6});
  constraint.SetDoRotationReference(doRotationReference);

  dp3::ddecal::SolutionTensor onesolution_tensor({1, 1, 1, 4});
  dp3::ddecal::SolutionSpan onesolution =
      aocommon::xt::CreateSpan(onesolution_tensor);
  double pi = 3.1415;
  double phi = pi / 6.0;
  const std::complex<double> i(0, 1.0);
  std::complex<double> a = 2.0 * exp(i * 0.3), b = 3.0 * exp(i * -0.2);

  /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
  */
  onesolution(0, 0, 0, 0) = a * std::cos(phi);
  onesolution(0, 0, 0, 1) = a * -std::sin(phi);
  onesolution(0, 0, 0, 2) = b * std::sin(phi);
  onesolution(0, 0, 0, 3) = b * std::cos(phi);

  std::vector<Constraint::Result> constraint_result;
  constraint_result = constraint.Apply(onesolution, 0.0, nullptr);
  BOOST_CHECK_EQUAL(constraint_result.size(), 3u);
  BOOST_CHECK_EQUAL(constraint_result[0].name, "rotation");
  BOOST_CHECK_EQUAL(constraint_result[0].axes, "ant,dir,freq");
  BOOST_CHECK_EQUAL(constraint_result[0].dims.size(), 3u);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[0], 1u);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[1], 1u);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[2], 1u);
  BOOST_CHECK_EQUAL(constraint_result[0].vals.size(), 1u);
  BOOST_CHECK(casacore::near(constraint_result[0].vals[0],
                             doRotationReference ? 0.0 : phi));

  BOOST_CHECK_EQUAL(constraint_result[1].name, "amplitude");
  BOOST_CHECK_EQUAL(constraint_result[1].axes, "ant,dir,freq,pol");
  BOOST_CHECK_EQUAL(constraint_result[1].dims.size(), 4u);
  BOOST_CHECK_EQUAL(constraint_result[1].dims[0], 1u);
  BOOST_CHECK_EQUAL(constraint_result[1].dims[1], 1u);
  BOOST_CHECK_EQUAL(constraint_result[1].dims[2], 1u);
  BOOST_CHECK_EQUAL(constraint_result[1].dims[3], 2u);
  BOOST_CHECK_EQUAL(constraint_result[1].vals.size(), 2u);
  BOOST_CHECK(casacore::near(constraint_result[1].vals[0], abs(a)));
  BOOST_CHECK(casacore::near(constraint_result[1].vals[1], abs(b)));

  BOOST_CHECK_EQUAL(constraint_result[2].name, "phase");
  BOOST_CHECK_EQUAL(constraint_result[2].axes, "ant,dir,freq,pol");
  BOOST_CHECK_EQUAL(constraint_result[2].dims.size(), 4u);
  BOOST_CHECK_EQUAL(constraint_result[2].dims[0], 1u);
  BOOST_CHECK_EQUAL(constraint_result[2].dims[1], 1u);
  BOOST_CHECK_EQUAL(constraint_result[2].dims[2], 1u);
  BOOST_CHECK_EQUAL(constraint_result[2].dims[3], 2u);
  BOOST_CHECK_EQUAL(constraint_result[2].vals.size(), 2u);
  BOOST_CHECK(casacore::near(constraint_result[2].vals[0], arg(a)));
  BOOST_CHECK(casacore::near(constraint_result[2].vals[1], arg(b)));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace dp3::ddecal
