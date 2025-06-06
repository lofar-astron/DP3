// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../constraints/AntennaConstraint.h"

#include <vector>
#include <complex>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <xtensor/xview.hpp>

using dp3::ddecal::AntennaConstraint;
using dp3::ddecal::Constraint;

BOOST_AUTO_TEST_SUITE(antennaconstraint)

BOOST_AUTO_TEST_CASE(test_antenna_constraint) {
  AntennaConstraint constraint;

  const size_t kNAntennas = 10;
  const size_t kNChannels = 2;
  const size_t kNSubSolutions = 1;
  const size_t kNPolarizations = 4;

  constraint.Initialize(kNAntennas, {1u}, {42.0e6, 84.0e6});
  const std::set<size_t> antenna_set1{0, 1, 2, 3, 4};
  const std::set<size_t> antenna_set2{5, 6};
  constraint.SetAntennaSets({antenna_set1, antenna_set2});

  dp3::ddecal::SolutionTensor onesolution(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});

  // Set one antenna in the first antenna group to a value (a separate one for
  // both channels and two polarizations)
  onesolution(0, 0, 0, 0) = std::complex(1.0, 0.0);
  onesolution(0, 0, 0, 3) = std::complex(17.0, 0.0);
  xt::view(onesolution, 1, 0, 0, xt::all()) = std::complex(3.0, 4.0);

  // Set one antenna in the second antenna group to a value (a separate one for
  // both channels)
  xt::view(onesolution, 0, 5, 0, xt::all()) = std::complex(1.0, 0.0);
  xt::view(onesolution, 1, 6, 0, xt::all()) = std::complex(3.0, 4.0);

  xt::view(onesolution, 0, 7, 0, xt::all()) = std::complex(1.0e8, 0.0);

  dp3::ddecal::SolutionSpan onesolution_span =
      aocommon::xt::CreateSpan(onesolution);
  std::vector<Constraint::Result> constraint_result =
      constraint.Apply(onesolution_span, 0.0, nullptr);

  BOOST_CHECK(constraint_result.empty());

  // Since only one antenna in the antenna group is nonzero, the mean is the
  // reciprocal of the one value (for each channel, for each polarization)
  for (size_t ant_nr : antenna_set1) {
    BOOST_CHECK_CLOSE(std::abs(onesolution(0, ant_nr, 0, 0)),
                      1.0 / antenna_set1.size(), 1.0e-8);
    BOOST_CHECK_CLOSE(std::abs(onesolution(0, ant_nr, 0, 3)),
                      17.0 / antenna_set1.size(), 1.0e-8);
    BOOST_CHECK_CLOSE(onesolution(1, ant_nr, 0, 0).real(),
                      3.0 / antenna_set1.size(), 1.0e-8);
    BOOST_CHECK_CLOSE(onesolution(1, ant_nr, 0, 0).imag(),
                      4.0 / antenna_set1.size(), 1.0e-8);
  }

  BOOST_CHECK_CLOSE(onesolution(0, 7, 0, 0).real(), 1.0e8, 1.0e-8);

  for (size_t ant_nr : antenna_set2) {
    BOOST_CHECK_CLOSE(std::abs(onesolution(0, ant_nr, 0, 0)),
                      1.0 / antenna_set2.size(), 1.0e-8);
    BOOST_CHECK_CLOSE(onesolution(1, ant_nr, 0, 0).real(),
                      3.0 / antenna_set2.size(), 1.0e-8);
    BOOST_CHECK_CLOSE(onesolution(1, ant_nr, 0, 0).imag(),
                      4.0 / antenna_set2.size(), 1.0e-8);
  }
}

BOOST_AUTO_TEST_SUITE_END()
