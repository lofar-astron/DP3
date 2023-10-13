// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../constraints/TECConstraint.h"

#include <vector>
#include <cmath>
#include <complex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <xtensor/xview.hpp>
#include <xtensor/xbuilder.hpp>  // linspace
#include <xtensor/xadapt.hpp>

using dp3::ddecal::ApproximateTECConstraint;
using dp3::ddecal::Constraint;
using dp3::ddecal::TECConstraint;

namespace {
const double kTecConstant = -8.44797245e9;

const size_t kNAntennas = 10;
const size_t kNChannels = 142;
const size_t kNSolutions = 1;
const size_t kNPolarizations = 1;

const size_t kApproximatingIterations = 4;

std::unique_ptr<TECConstraint> CreateConstraint(TECConstraint::Mode mode,
                                                bool approximate_tec) {
  if (approximate_tec) {
    auto constraint = std::make_unique<ApproximateTECConstraint>(mode);
    constraint->SetMaxApproximatingIterations(kApproximatingIterations);
    return constraint;
  } else {
    return std::make_unique<TECConstraint>(mode);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(tecconstraint)

BOOST_DATA_TEST_CASE(tec_only, boost::unit_test::data::make({false, true}),
                     approximate_tec) {
  std::unique_ptr<TECConstraint> constraint =
      CreateConstraint(TECConstraint::Mode::kTecOnly, approximate_tec);

  std::vector<double> channel_frequencies(kNChannels);
  xt::adapt(channel_frequencies) = xt::linspace(120.0e6, 150.0e6, kNChannels);
  constraint->Initialize(kNAntennas, {1u}, channel_frequencies);

  dp3::ddecal::SolutionTensor onesolution(
      {kNChannels, kNAntennas, kNSolutions, kNPolarizations});

  const xt::xtensor<double, 1> tec_values = xt::linspace(0.0, 3.0, kNAntennas);

  xt::view(onesolution, xt::all(), xt::all(), 0, 0) = xt::exp(
      std::complex<double>(0, 1.0) * tec_values * kTecConstant /
      xt::view(xt::adapt(channel_frequencies), xt::all(), xt::newaxis()));

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(onesolution);
  std::vector<Constraint::Result> constraint_result;

  if (approximate_tec) {
    // Do an approximation, which yields no results.
    constraint->PrepareIteration(false, 0, false);
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_CHECK(constraint_result.empty());

    // Tell the constraint the final iteration is next.
    constraint->PrepareIteration(false, 0, true);
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_CHECK_EQUAL(constraint_result.size(), 2);

    // TODO: Verify the result for ApproximateTECConstraint.
  } else {
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_REQUIRE_EQUAL(constraint_result.size(), 2);

    const Constraint::Result& tec_result = constraint_result[0];
    const Constraint::Result& error_result = constraint_result[1];

    std::vector<size_t> expected_dims = {kNAntennas, kNSolutions, 1};

    BOOST_CHECK_EQUAL(tec_result.name, "tec");
    BOOST_CHECK_EQUAL(tec_result.axes, "ant,dir,freq");
    BOOST_CHECK_EQUAL_COLLECTIONS(tec_result.dims.begin(),
                                  tec_result.dims.end(), expected_dims.begin(),
                                  expected_dims.end());
    BOOST_REQUIRE_EQUAL(tec_result.vals.size(), kNAntennas);
    BOOST_CHECK_SMALL(tec_result.vals[0], 1.0e-6);

    BOOST_CHECK_EQUAL(error_result.name, "error");
    BOOST_CHECK_EQUAL(error_result.axes, "ant,dir,freq");
    BOOST_CHECK_EQUAL_COLLECTIONS(error_result.dims.begin(),
                                  error_result.dims.end(),
                                  expected_dims.begin(), expected_dims.end());
    BOOST_REQUIRE_EQUAL(error_result.vals.size(), kNAntennas);
    BOOST_CHECK_SMALL(error_result.vals[0], 1.0e-6);

    // Only test that the weights are all equal, the exact value doesn't matter
    const double weight0 = tec_result.weights[0];
    BOOST_CHECK_CLOSE(error_result.weights[0], weight0, 1.0e-6);
    for (size_t ant = 1; ant < kNAntennas; ++ant) {
      BOOST_CHECK_CLOSE(tec_result.vals[ant], tec_values[ant], 1.0e-3);
      BOOST_CHECK_SMALL(error_result.vals[ant], 1.0e-6);
      BOOST_CHECK_CLOSE(error_result.weights[ant], weight0, 1.0e-6);
      BOOST_CHECK_CLOSE(tec_result.weights[ant], weight0, 1.0e-6);
    }
  }
}

BOOST_DATA_TEST_CASE(tec_and_phase, boost::unit_test::data::make({false, true}),
                     approximate_tec) {
  std::unique_ptr<TECConstraint> constraint = CreateConstraint(
      TECConstraint::Mode::kTecAndCommonScalar, approximate_tec);

  std::vector<double> channel_frequencies(kNChannels);
  xt::adapt(channel_frequencies) = xt::linspace(120.0e6, 150.0e6, kNChannels);
  constraint->Initialize(kNAntennas, {1u}, channel_frequencies);

  dp3::ddecal::SolutionTensor onesolution(
      {kNChannels, kNAntennas, kNSolutions, kNPolarizations});

  const double kScalarPhase0 = 1.23;
  const xt::xtensor<double, 1> tec_values = xt::linspace(0.0, 3.0, kNAntennas);
  const xt::xtensor<double, 1> scalar_phases =
      xt::linspace(kScalarPhase0, 10.0, kNAntennas);

  xt::view(onesolution, xt::all(), xt::all(), 0, 0) =
      xt::exp(std::complex<double>(0, 1.0) *
              (scalar_phases + (tec_values * kTecConstant /
                                xt::view(xt::adapt(channel_frequencies),
                                         xt::all(), xt::newaxis()))));

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(onesolution);
  std::vector<Constraint::Result> constraint_result;

  if (approximate_tec) {
    // Do an approximation, which yields no results.
    constraint->PrepareIteration(false, 0, false);
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_CHECK(constraint_result.empty());

    // Tell the constraint the final iteration is next.
    constraint->PrepareIteration(false, 0, true);
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_CHECK_EQUAL(constraint_result.size(), 3);

    // TODO: Verify the result for ApproximateTECConstraint.
  } else {
    constraint_result = constraint->Apply(onesolution_span, 0.0, nullptr);
    BOOST_REQUIRE_EQUAL(constraint_result.size(), 3);

    const Constraint::Result& tec_result = constraint_result[0];
    const Constraint::Result& phase_result = constraint_result[1];
    const Constraint::Result& error_result = constraint_result[2];

    std::vector<size_t> expected_dims = {kNAntennas, kNSolutions, 1};

    BOOST_CHECK(tec_result.name == "tec");
    BOOST_CHECK_EQUAL(tec_result.axes, "ant,dir,freq");
    BOOST_CHECK_EQUAL_COLLECTIONS(tec_result.dims.begin(),
                                  tec_result.dims.end(), expected_dims.begin(),
                                  expected_dims.end());
    BOOST_CHECK(tec_result.vals.size() == kNAntennas);
    BOOST_CHECK(phase_result.name == "phase");
    BOOST_CHECK_EQUAL(phase_result.axes, "ant,dir,freq");
    BOOST_CHECK_EQUAL_COLLECTIONS(phase_result.dims.begin(),
                                  phase_result.dims.end(),
                                  expected_dims.begin(), expected_dims.end());
    BOOST_CHECK(phase_result.vals.size() == kNAntennas);
    BOOST_CHECK_SMALL(tec_result.vals[0], 1.0e-6);
    // Only test that the weights are all equal, the exact value doesn't matter
    const double weight0 = tec_result.weights[0];
    BOOST_CHECK_CLOSE(phase_result.weights[0], weight0, 1.0e-3);
    BOOST_CHECK_CLOSE(error_result.weights[0], weight0, 1.0e-3);

    BOOST_CHECK_SMALL(phase_result.vals[0], 1.0e-6);
    BOOST_CHECK_SMALL(error_result.vals[0], 1.0e-6);

    for (size_t ant = 1; ant < kNAntennas; ++ant) {
      BOOST_CHECK_CLOSE(tec_result.vals[ant], tec_values[ant], 1.0e-3);
      // Add a big number of 2pi to make sure that phase_result is positive
      // Put number between -pi and pi
      // Use station 0 as phase reference
      const double result_phase_referenced = scalar_phases[ant] - kScalarPhase0;
      BOOST_CHECK_SMALL(
          std::fmod(phase_result.vals[ant] - result_phase_referenced +
                        100 * M_PI + M_PI,
                    2.0 * M_PI) -
              M_PI,
          1.0e-3);
      BOOST_CHECK_SMALL(error_result.vals[ant], 1.0e-6);
      BOOST_CHECK_CLOSE(tec_result.weights[ant], weight0, 1.0e-3);
      BOOST_CHECK_CLOSE(phase_result.weights[ant], weight0, 1.0e-3);
      BOOST_CHECK_CLOSE(error_result.weights[ant], weight0, 1.0e-3);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
