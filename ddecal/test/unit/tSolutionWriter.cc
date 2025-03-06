// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../../SolutionWriter.h"

#include <iomanip>
#include <iostream>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <boost/test/execution_monitor.hpp>

#include <schaapcommon/h5parm/h5parm.h>

BOOST_AUTO_TEST_SUITE(write_solution)

using dp3::ddecal::SolutionWriter;

BOOST_AUTO_TEST_CASE(solution_time_size) {
  /**
   * The solutions container has a length of 200 time slots.
   * The measurement set time interval is 5 seconds.
   * The solution time interval is 50 seconds.
   * The solutions' time span is 10000 seconds.
   * However, the difference between the end and start time is 5000 seconds.
   * Accordingly, the H5 file time axis must include only the times that fit
   * into that time range and discard the rest.
   **/

  const int kNChannelBlocks = 10;
  const int kNSolutions = 200;
  const int kNAntennas = 3;
  const int kNDirections = 10;
  const int kBlockSize = kNAntennas * kNDirections;
  const int kExpectedTimeSize = 100;

  {
    // Scalar arguments
    const double kStartTime = 5000000000;
    const double kEndTime = 5000005000;
    const double kMsInterval = 5.0;
    const size_t kTimestepsPerSolutionInterval = 10;
    const dp3::base::CalType mode = dp3::base::CalType::kScalar;

    // Empty container argument
    std::vector<std::vector<std::vector<dp3::ddecal::Constraint::Result>>>
        constraint_solutions;

    // Comment in H5 file
    std::string history = "DP3 testing -- solution_time_size";

    // Prepare several container arguments
    // and write the solution into an H5 file.

    std::vector<std::vector<std::vector<std::complex<double>>>> solutions(
        kNSolutions);
    for (size_t i = 0; i < solutions.size(); i++) {
      solutions[i].resize(kNChannelBlocks);
      for (size_t j = 0; j < solutions[i].size(); j++) {
        solutions[i][j].resize(kBlockSize);
        for (size_t k = 0; k < solutions[i][j].size(); k++) {
          solutions[i][j][k] = 0.2;
        }
      }
    }

    std::vector<std::string> antenna_names(kNAntennas);
    for (size_t i = 0; i < antenna_names.size(); i++) {
      std::stringstream text;
      text << "CS" << std::setfill('0') << std::setw(7) << i;
      antenna_names[i] = text.str();
    }

    std::vector<dp3::base::Direction> source_directions(kNDirections);
    for (size_t i = 0; i < source_directions.size(); i++) {
      source_directions[i] = dp3::base::Direction(0.1, 0.1);
    }

    std::vector<std::vector<std::string>> directions(kNDirections);
    for (size_t i = 0; i < directions.size(); i++) {
      directions[i].resize(1);
    }

    std::vector<double> chan_freqs(kNChannelBlocks);
    for (size_t i = 0; i < chan_freqs.size(); i++) {
      chan_freqs[i] = 1.200e+08 + i * 1.0e6;
    }

    std::vector<double> chan_block_freqs(kNChannelBlocks);
    for (size_t i = 0; i < chan_block_freqs.size(); i++) {
      chan_block_freqs[i] = chan_freqs[i];
    }

    const std::vector<size_t> solutions_per_direction(kNDirections, 1);

    SolutionWriter solution_writer("output.h5");
    solution_writer.Write(solutions, constraint_solutions, kStartTime, kEndTime,
                          kMsInterval, kTimestepsPerSolutionInterval,
                          solutions_per_direction, mode, antenna_names,
                          source_directions, directions, chan_freqs,
                          chan_block_freqs, history);
  }

  schaapcommon::h5parm::H5Parm h5_parm("output.h5");
  schaapcommon::h5parm::AxisInfo axis_info =
      h5_parm.GetSolTab("amplitude000").GetAxis("time");

  BOOST_CHECK_EQUAL(kExpectedTimeSize, axis_info.size);
}

BOOST_AUTO_TEST_SUITE_END()
