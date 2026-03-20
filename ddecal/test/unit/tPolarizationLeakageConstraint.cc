#include <complex>

#include "ddecal/constraints/PolarizationLeakageConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

using dp3::ddecal::Constraint;
using dp3::ddecal::ConstraintResult;
using dp3::ddecal::PolarizationLeakageConstraint;
using std::complex;
using std::vector;

BOOST_AUTO_TEST_SUITE(polarization_leakage_constraint)

BOOST_AUTO_TEST_CASE(apply_constraint) {
  constexpr size_t kNCorrelations = 4;
  constexpr size_t kNAntennas = 2;
  constexpr size_t kNSubSolutions = 3;
  constexpr size_t kNChannels = 2;
  const std::vector<uint32_t> solutions_per_direction{1, 2};
  const std::vector<double> frequencies{123.0, 124.0};

  PolarizationLeakageConstraint constraint;
  constraint.Initialize(kNAntennas, solutions_per_direction, frequencies);

  dp3::ddecal::SolutionTensor solutions_tensor(
      {kNChannels, kNAntennas, kNSubSolutions, kNCorrelations});
  dp3::ddecal::SolutionSpan solutions =
      aocommon::xt::CreateSpan(solutions_tensor);

  size_t counter = 1;
  for (std::complex<double>& solution : solutions_tensor) {
    solution = std::complex<double>(counter, -0.5);
    ++counter;
  }

  constexpr double kTime = 42.0;
  constraint.Apply(solutions, kTime, nullptr);

  counter = 1;
  for (size_t ch = 0; ch < kNChannels; ++ch) {
    for (size_t ant = 0; ant < kNAntennas; ++ant) {
      for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
           ++sub_solution) {
        const std::complex<double>* matrix =
            &solutions(ch, ant, sub_solution, 0);
        BOOST_CHECK_EQUAL(matrix[0].real(), 1.0);
        BOOST_CHECK_EQUAL(matrix[0].imag(), 0.0);
        BOOST_CHECK_CLOSE(matrix[1].real(), counter + 1, 1e-6);
        BOOST_CHECK_EQUAL(matrix[1].imag(), -0.5);
        BOOST_CHECK_CLOSE(matrix[2].real(), counter + 2, 1e-6);
        BOOST_CHECK_EQUAL(matrix[2].imag(), -0.5);
        BOOST_CHECK_EQUAL(matrix[3].real(), 1.0);
        BOOST_CHECK_EQUAL(matrix[3].imag(), 0.0);
        counter += 4;
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
