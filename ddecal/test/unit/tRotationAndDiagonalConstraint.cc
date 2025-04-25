// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <vector>
#include <complex>

#include "../../constraints/RotationAndDiagonalConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace dp3::ddecal {

using base::CalType;

namespace {
constexpr double kReferencePhi = M_PI / 6.0;
constexpr size_t kNDirections = 3;
}  // namespace

BOOST_AUTO_TEST_SUITE(rotation_and_diagonal_constraint)

BOOST_AUTO_TEST_CASE(constrain_diagonal) {
  const std::array<std::complex<double>, 2> diagonal{std::complex{4.0, -3.0},
                                                     std::complex{21.0, -20.0}};
  std::array<std::complex<double>, 2> result = diagonal;
  ConstrainDiagonal(result, CalType::kDiagonal);
  BOOST_CHECK(diagonal == result);

  result = diagonal;
  ConstrainDiagonal(result, CalType::kDiagonalAmplitude);
  BOOST_CHECK_CLOSE_FRACTION(result[0].real(), 5.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[0].imag(), 0.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].real(), 29.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].imag(), 0.0, 1e-8);

  result = diagonal;
  ConstrainDiagonal(result, CalType::kDiagonalPhase);
  BOOST_CHECK_CLOSE_FRACTION(result[0].real(), 4.0 / 5.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[0].imag(), -3.0 / 5.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].real(), 21.0 / 29.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].imag(), -20.0 / 29.0, 1e-8);

  const std::complex<double> scalar((diagonal[0] + diagonal[1]) * 0.5);

  result = diagonal;
  ConstrainDiagonal(result, CalType::kScalar);
  BOOST_CHECK_CLOSE_FRACTION(result[0].real(), scalar.real(), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[0].imag(), scalar.imag(), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].real(), scalar.real(), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].imag(), scalar.imag(), 1e-8);

  result = diagonal;
  ConstrainDiagonal(result, CalType::kScalarAmplitude);
  BOOST_CHECK_CLOSE_FRACTION(result[0].real(), std::abs(scalar), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[0].imag(), 0.0, 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].real(), std::abs(scalar), 1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].imag(), 0.0, 1e-8);

  result = diagonal;
  ConstrainDiagonal(result, CalType::kScalarPhase);
  BOOST_CHECK_CLOSE_FRACTION(result[0].real(), scalar.real() / std::abs(scalar),
                             1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[0].imag(), scalar.imag() / std::abs(scalar),
                             1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].real(), scalar.real() / std::abs(scalar),
                             1e-8);
  BOOST_CHECK_CLOSE_FRACTION(result[1].imag(), scalar.imag() / std::abs(scalar),
                             1e-8);
}

std::vector<Constraint::Result> TestApplyConstraint(
    RotationAndDiagonalConstraint& constraint, bool enable_rotation_reference,
    std::complex<double> a, std::complex<double> b) {
  constraint.Initialize(1, std::vector<uint32_t>(kNDirections, 1u), {42e6});
  constraint.SetDoRotationReference(enable_rotation_reference);

  dp3::ddecal::SolutionTensor solutions_tensor({1, 1, kNDirections, 4});
  dp3::ddecal::SolutionSpan solutions =
      aocommon::xt::CreateSpan(solutions_tensor);
  /* Solution is of the form ((a,0),(0,b))*rot(phi)
     with rot(phi) = ((cos(phi),-sin(phi)),(sin(phi),cos(phi)))
  */
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    solutions(0, 0, direction, 0) =
        a * std::cos(kReferencePhi + 0.1 * direction);
    solutions(0, 0, direction, 1) =
        a * -std::sin(kReferencePhi + 0.1 * direction);
    solutions(0, 0, direction, 2) =
        b * std::sin(kReferencePhi + 0.1 * direction);
    solutions(0, 0, direction, 3) =
        b * std::cos(kReferencePhi + 0.1 * direction);
  }

  std::vector<Constraint::Result> constraint_result;
  constraint_result = constraint.Apply(solutions, 0.0, nullptr);
  BOOST_CHECK_GE(constraint_result.size(), 2u);
  BOOST_CHECK_EQUAL(constraint_result[0].name, "rotation");
  BOOST_CHECK_EQUAL(constraint_result[0].axes, "ant,dir,freq");
  BOOST_REQUIRE_EQUAL(constraint_result[0].dims.size(), 3u);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[0], 1u);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[1], kNDirections);
  BOOST_CHECK_EQUAL(constraint_result[0].dims[2], 1u);
  BOOST_REQUIRE_EQUAL(constraint_result[0].vals.size(), kNDirections);
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    BOOST_CHECK_CLOSE_FRACTION(
        constraint_result[0].vals[direction],
        enable_rotation_reference ? 0.0 : kReferencePhi + 0.1 * direction,
        1e-6);
  }
  return constraint_result;
}

void CheckDiagonalAmplitude(const Constraint::Result& result,
                            std::complex<double> a, std::complex<double> b) {
  BOOST_CHECK_EQUAL(result.name, "amplitude");
  BOOST_CHECK_EQUAL(result.axes, "ant,dir,freq,pol");
  BOOST_REQUIRE_EQUAL(result.dims.size(), 4u);
  BOOST_CHECK_EQUAL(result.dims[0], 1u);
  BOOST_CHECK_EQUAL(result.dims[1], kNDirections);
  BOOST_CHECK_EQUAL(result.dims[2], 1u);
  BOOST_CHECK_EQUAL(result.dims[3], 2u);
  BOOST_REQUIRE_EQUAL(result.vals.size(), kNDirections * 2);
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction * 2], std::abs(a), 1e-4);
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction * 2 + 1], std::abs(b),
                               1e-4);
  }
}

void CheckDiagonalPhase(const Constraint::Result& result,
                        std::complex<double> a, std::complex<double> b) {
  BOOST_CHECK_EQUAL(result.name, "phase");
  BOOST_CHECK_EQUAL(result.axes, "ant,dir,freq,pol");
  BOOST_REQUIRE_EQUAL(result.dims.size(), 4u);
  BOOST_CHECK_EQUAL(result.dims[0], 1u);
  BOOST_CHECK_EQUAL(result.dims[1], kNDirections);
  BOOST_CHECK_EQUAL(result.dims[2], 1u);
  BOOST_CHECK_EQUAL(result.dims[3], 2u);
  BOOST_REQUIRE_EQUAL(result.vals.size(), kNDirections * 2);
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction * 2], std::arg(a), 1e-4);
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction * 2 + 1], std::arg(b),
                               1e-4);
  }
}

void CheckScalarAmplitude(const Constraint::Result& result,
                          std::complex<double> a) {
  BOOST_CHECK_EQUAL(result.name, "amplitude");
  BOOST_CHECK_EQUAL(result.axes, "ant,dir,freq");
  BOOST_REQUIRE_EQUAL(result.dims.size(), 3u);
  BOOST_CHECK_EQUAL(result.dims[0], 1u);
  BOOST_CHECK_EQUAL(result.dims[1], kNDirections);
  BOOST_CHECK_EQUAL(result.dims[2], 1u);
  BOOST_REQUIRE_EQUAL(result.vals.size(), kNDirections);
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction], std::abs(a), 1e-4);
  }
}

void CheckScalarPhase(const Constraint::Result& result,
                      std::complex<double> a) {
  BOOST_CHECK_EQUAL(result.name, "phase");
  BOOST_CHECK_EQUAL(result.axes, "ant,dir,freq");
  BOOST_REQUIRE_EQUAL(result.dims.size(), 3u);
  BOOST_CHECK_EQUAL(result.dims[0], 1u);
  BOOST_CHECK_EQUAL(result.dims[1], kNDirections);
  BOOST_CHECK_EQUAL(result.dims[2], 1u);
  for (size_t direction = 0; direction != kNDirections; ++direction) {
    BOOST_CHECK_CLOSE_FRACTION(result.vals[direction], std::arg(a), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal) {
  RotationAndDiagonalConstraint constraint(CalType::kDiagonal);
  const std::complex<double> kA = std::polar(2.0, 0.3);
  const std::complex<double> kB = std::polar(3.0, -0.2);
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kB);
  BOOST_REQUIRE_EQUAL(results.size(), 3u);
  CheckDiagonalAmplitude(results[1], kA, kB);
  CheckDiagonalPhase(results[2], kA, kB);
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal_with_reference) {
  RotationAndDiagonalConstraint constraint(CalType::kDiagonal);
  const std::complex<double> kA = std::polar(2.0, 0.3);
  const std::complex<double> kB = std::polar(3.0, -0.2);
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, true, kA, kB);
  BOOST_REQUIRE_EQUAL(results.size(), 3u);
  CheckDiagonalAmplitude(results[1], kA, kB);
  CheckDiagonalPhase(results[2], kA, kB);
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal_amplitude) {
  RotationAndDiagonalConstraint constraint(CalType::kDiagonalAmplitude);
  const std::complex<double> kA = 2.0;
  const std::complex<double> kB = 3.0;
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kB);
  BOOST_REQUIRE_EQUAL(results.size(), 2u);
  CheckDiagonalAmplitude(results[1], kA, kB);
}

BOOST_AUTO_TEST_CASE(rotation_and_diagonal_phase) {
  RotationAndDiagonalConstraint constraint(CalType::kDiagonalPhase);
  const std::complex<double> kA = std::polar(1.0, 0.3);
  const std::complex<double> kB = std::polar(1.0, -0.2);
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kB);
  BOOST_REQUIRE_EQUAL(results.size(), 2u);
  CheckDiagonalPhase(results[1], kA, kB);
}

BOOST_AUTO_TEST_CASE(rotation_and_scalar) {
  RotationAndDiagonalConstraint constraint(CalType::kScalar);
  const std::complex<double> kA = std::polar(2.0, 0.3);
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kA);
  BOOST_REQUIRE_EQUAL(results.size(), 3u);
  CheckScalarAmplitude(results[1], kA);
  CheckScalarPhase(results[2], kA);
}

BOOST_AUTO_TEST_CASE(rotation_and_scalar_amplitude) {
  RotationAndDiagonalConstraint constraint(CalType::kScalarAmplitude);
  const std::complex<double> kA = 2.0;
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kA);
  BOOST_REQUIRE_EQUAL(results.size(), 2u);
  CheckScalarAmplitude(results[1], kA);
}

BOOST_AUTO_TEST_CASE(rotation_and_scalar_phase) {
  RotationAndDiagonalConstraint constraint(CalType::kScalarPhase);
  const std::complex<double> kA = std::polar(1.0, 0.3);
  std::vector<Constraint::Result> results =
      TestApplyConstraint(constraint, false, kA, kA);
  BOOST_REQUIRE_EQUAL(results.size(), 2u);
  CheckScalarPhase(results[1], kA);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace dp3::ddecal
