// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SolutionResampler.h"

#include <iostream>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <boost/test/execution_monitor.hpp>

BOOST_AUTO_TEST_SUITE(solution_resampler)

using dp3::ddecal::ConstraintResult;
using dp3::ddecal::SolutionResampler;

BOOST_AUTO_TEST_CASE(sub_steps) {
  const size_t solution_interval = 6;

  {
    const std::vector<size_t> n_solutions_per_direction = {1, 3};
    SolutionResampler resample(n_solutions_per_direction, 1, solution_interval);
    BOOST_CHECK_EQUAL(resample.GetNrSubSteps(), 3u);
  }

  {
    const std::vector<size_t> n_solutions_per_direction = {2, 3};
    SolutionResampler resample(n_solutions_per_direction, 1, solution_interval);
    BOOST_CHECK_EQUAL(resample.GetNrSubSteps(), solution_interval);
  }
}

BOOST_AUTO_TEST_CASE(upsample_normal_solutions) {
  const size_t solution_interval = 6;
  const size_t n_time_steps = 1;
  const size_t n_channel_blocks = 1;
  const size_t n_antenna = 2;
  // Values must be integer divisor of solution_interval
  const std::vector<size_t> n_solutions_per_direction = {2, 3, 1};
  const size_t n_pol = 1;

  SolutionResampler resample(n_solutions_per_direction, n_antenna,
                             solution_interval);

  const size_t n_directions = n_solutions_per_direction.size();
  const size_t n_solutions = std::accumulate(
      n_solutions_per_direction.begin(), n_solutions_per_direction.end(), 0u);

  std::vector<std::vector<std::vector<std::complex<double>>>> original_solution(
      n_time_steps);
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

  const std::vector<std::vector<std::vector<std::complex<double>>>>
      upsampled_solution = resample.Upsample(original_solution, n_pol);

  // 6 is divisible by all values in n_solutions_per_direction
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

BOOST_AUTO_TEST_CASE(upsample_constraint_results) {
  // Least common multiple of n_solutions_per_direction
  constexpr size_t solution_interval = 15;
  constexpr size_t n_intervals = 3;
  constexpr size_t n_channel_blocks = 4;
  constexpr size_t n_antenna = 2;
  constexpr size_t n_constraints = 2;
  constexpr size_t n_constraint_sub_types = 2;
  // Values must be integer divisor of solution_interval
  const std::vector<size_t> n_solutions_per_direction = {3, 1, 5};

  SolutionResampler resample(n_solutions_per_direction, n_antenna,
                             solution_interval);

  const size_t n_directions = n_solutions_per_direction.size();
  const size_t n_sub_solutions = std::accumulate(
      n_solutions_per_direction.begin(), n_solutions_per_direction.end(), 0u);

  std::vector<std::vector<std::vector<ConstraintResult>>> input(n_intervals);
  for (size_t interval = 0; interval != n_intervals; ++interval) {
    input[interval].resize(n_constraints);
    for (size_t constraint = 0; constraint != n_constraints; ++constraint) {
      input[interval][constraint].resize(n_constraint_sub_types);
      for (size_t subtype = 0; subtype != n_constraint_sub_types; ++subtype) {
        ConstraintResult& input_result = input[interval][constraint][subtype];
        input_result.dims = {n_antenna, n_sub_solutions, n_channel_blocks};
        input_result.vals.resize(n_antenna * n_sub_solutions *
                                 n_channel_blocks);
        input_result.name = "the name";
        input_result.axes = "ant,dir,freq";
        std::iota(input_result.vals.begin(), input_result.vals.end(), 0.0);
        input_result.weights.assign(
            n_antenna * n_sub_solutions * n_channel_blocks, 1.0);
      }
    }
  }

  const std::vector<std::vector<std::vector<ConstraintResult>>>
      upsampled_solution = resample.Upsample(input);

  // Least common multiple of n_solutions_per_direction
  BOOST_REQUIRE_EQUAL(upsampled_solution.size(),
                      solution_interval * n_intervals);
  for (const auto& time_solution : upsampled_solution) {
    BOOST_REQUIRE_EQUAL(time_solution.size(), n_constraints);
    for (const auto& constraint_solution : time_solution) {
      BOOST_REQUIRE_EQUAL(constraint_solution.size(), n_constraint_sub_types);
      for (const ConstraintResult& result : constraint_solution) {
        BOOST_REQUIRE_EQUAL(result.name, "the name");
        BOOST_REQUIRE_EQUAL(result.axes, "ant,dir,freq");
        BOOST_REQUIRE_EQUAL(result.dims.size(), 3);
        const size_t count = n_antenna * n_directions * n_channel_blocks;
        BOOST_REQUIRE_EQUAL(result.vals.size(), count);
        BOOST_REQUIRE_EQUAL(result.weights.size(), count);
        BOOST_CHECK_EQUAL(result.weights[0], 1.0);
      }
    }
  }

  const std::vector<double> solutions0_to_3 = {
      0.0,  1.0,
      2.0,  3.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      16.0, 17.0,
      18.0, 19.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      36.0, 37.0,
      38.0, 39.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      52.0, 53.0,
      54.0, 55.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 0; i != 3; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions0_to_3.begin(),
                                  solutions0_to_3.end());
  }

  // At timestep 3, the fastest direction, direction 2, changes because it
  // reduces the 15 timesteps per intervals by a factor of 5.
  const std::vector<double> solutions3_to_5 = {
      0.0,  1.0,
      2.0,  3.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      20.0, 21.0,
      22.0, 23.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      36.0, 37.0,
      38.0, 39.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      56.0, 57.0,
      58.0, 59.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 3; i != 5; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions3_to_5.begin(),
                                  solutions3_to_5.end());
  }
  // At timestep 5, direction 0 changes because it reduces the
  // 15 timesteps per intervals by a factor of 3.
  const std::vector<double> solutions5_to_6 = {
      4.0,  5.0,
      6.0,  7.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      20.0, 21.0,
      22.0, 23.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      40.0, 41.0,
      42.0, 43.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      56.0, 57.0,
      58.0, 59.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 5; i != 6; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions5_to_6.begin(),
                                  solutions5_to_6.end());
  }
  // At timestep 6, direction 2 changes again.
  const std::vector<double> solutions6_to_9 = {
      4.0,  5.0,
      6.0,  7.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      24.0, 25.0,
      26.0, 27.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      40.0, 41.0,
      42.0, 43.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      60.0, 61.0,
      62.0, 63.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 6; i != 9; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions6_to_9.begin(),
                                  solutions6_to_9.end());
  }
  // At timestep 9, direction 2 changes again.
  const std::vector<double> solutions9_to_10 = {
      4.0,  5.0,
      6.0,  7.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      28.0, 29.0,
      30.0, 31.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      40.0, 41.0,
      42.0, 43.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      64.0, 65.0,
      66.0, 67.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 9; i != 10; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions9_to_10.begin(),
                                  solutions9_to_10.end());
  }
  // At timestep 10, direction 0 changes
  const std::vector<double> solutions10_to_12 = {
      8.0,  9.0,
      10.0, 11.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      28.0, 29.0,
      30.0, 31.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      44.0, 45.0,
      46.0, 47.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      64.0, 65.0,
      66.0, 67.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 10; i != 12; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions10_to_12.begin(),
                                  solutions10_to_12.end());
  }
  // At timestep 12, direction 2 changes
  const std::vector<double> solutions12_to_15 = {
      8.0,  9.0,
      10.0, 11.0,  // ant 0 direction 0 with 3 subsolutions and 4 channels
      12.0, 13.0,
      14.0, 15.0,  // ant 0 direction 1 with 1 subsolution and 4 channels
      32.0, 33.0,
      34.0, 35.0,  // ant 0 direction 2 with 5 subsolutions and 4 channels
      44.0, 45.0,
      46.0, 47.0,  // ant 1 direction 0 with 3 subsolutions and 4 channels
      48.0, 49.0,
      50.0, 51.0,  // ant 1 direction 1 with 1 subsolution and 4 channels
      68.0, 69.0,
      70.0, 71.0  // ant 1 direction 2 with 5 subsolutions and 4 channels
  };
  for (size_t i = 12; i != 15; ++i) {
    const ConstraintResult& result = upsampled_solution[i][0][0];
    BOOST_CHECK_EQUAL_COLLECTIONS(result.vals.begin(), result.vals.end(),
                                  solutions12_to_15.begin(),
                                  solutions12_to_15.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()
