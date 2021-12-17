// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SolutionResampler.h"

#include <iostream>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <boost/test/execution_monitor.hpp>

BOOST_AUTO_TEST_SUITE(resample_solution)

using dp3::ddecal::SolutionResampler;

BOOST_AUTO_TEST_CASE(sub_steps) {
  const size_t solution_interval = 6;

  {
    const std::vector<size_t> n_solutions_per_direction = {1, 3};
    SolutionResampler resample(n_solutions_per_direction, 1, 1,
                               solution_interval);
    BOOST_CHECK_EQUAL(resample.GetNrSubSteps(), 3u);
  }

  {
    const std::vector<size_t> n_solutions_per_direction = {2, 3};
    SolutionResampler resample(n_solutions_per_direction, 1, 1,
                               solution_interval);
    BOOST_CHECK_EQUAL(resample.GetNrSubSteps(), solution_interval);
  }
}

BOOST_AUTO_TEST_CASE(upsample) {
  const size_t solution_interval = 6;
  const size_t n_time_steps = 1;
  const size_t n_channel_blocks = 1;
  const size_t n_antenna = 2;
  // Values must be integer divisor of solution_interval
  const std::vector<size_t> n_solutions_per_direction = {2, 3, 1};
  const size_t n_pol = 1;

  SolutionResampler resample(n_solutions_per_direction, n_antenna, n_pol,
                             solution_interval);

  const size_t n_directions = n_solutions_per_direction.size();
  const size_t n_solutions = std::accumulate(
      n_solutions_per_direction.begin(), n_solutions_per_direction.end(), 0u);

  std::vector<std::vector<std::vector<std::complex<double>>>> original_solution;
  original_solution.resize(n_time_steps);
  for (auto& time_solution : original_solution) {
    time_solution.resize(n_channel_blocks);
    for (auto& channel_solution : time_solution) {
      for (size_t index = 0; index != n_antenna * n_solutions * n_pol;
           ++index) {
        // Use +1 so that solution vector does not contain zeros.
        channel_solution.push_back(index + 1);
      }
    }
  }

  std::vector<std::vector<std::vector<std::complex<double>>>>
      upsampled_solution = resample.Upsample(original_solution);

  // Least common multiple of n_solutions_per_direction
  const size_t n_substeps = 6u;
  BOOST_CHECK_EQUAL(upsampled_solution.size(), n_substeps);
  for (const auto& time_solution : upsampled_solution) {
    BOOST_CHECK_EQUAL(time_solution.size(), n_channel_blocks);
    for (const auto& channel_solution : time_solution) {
      BOOST_CHECK_EQUAL(channel_solution.size(),
                        n_antenna * n_directions * n_pol);
    }
  }

  const std::vector<std::vector<std::complex<double>>> ref_solution = {
      {1.0, 3.0, 6.0, 7.0, 9.0, 12.0},  {1.0, 3.0, 6.0, 7.0, 9.0, 12.0},
      {1.0, 4.0, 6.0, 7.0, 10.0, 12.0}, {2.0, 4.0, 6.0, 8.0, 10.0, 12.0},
      {2.0, 5.0, 6.0, 8.0, 11.0, 12.0}, {2.0, 5.0, 6.0, 8.0, 11.0, 12.0},
  };

  for (size_t substep = 0; substep != n_substeps; ++substep) {
    BOOST_CHECK_EQUAL_COLLECTIONS(upsampled_solution[substep][0].begin(),
                                  upsampled_solution[substep][0].end(),
                                  ref_solution[substep].begin(),
                                  ref_solution[substep].end());
  }
}

BOOST_AUTO_TEST_SUITE_END()
