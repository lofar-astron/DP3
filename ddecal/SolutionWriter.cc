// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionWriter.h"

#include "../common/StreamUtil.h"
#include "../common/StringTools.h"

#include <cassert>
#include <numeric>
#include <algorithm>

using dp3::base::CalType;
using dp3::common::operator<<;

namespace {
std::vector<std::string> GetDirectionNames(
    const std::vector<std::vector<std::string>>& directions) {
  std::vector<std::string> result;
  result.reserve(directions.size());

  for (const std::vector<std::string>& direction : directions) {
    std::stringstream ss;
    ss << direction;
    result.emplace_back(ss.str());
  }

  return result;
}
}  // namespace

namespace dp3 {
namespace ddecal {

SolutionWriter::SolutionWriter(const std::string& filename)
    : h5parm_(filename, true) {}

void SolutionWriter::AddAntennas(
    const std::vector<std::string>& all_antenna_names,
    const std::vector<std::array<double, 3>>& all_antenna_positions) {
  h5parm_.AddAntennas(all_antenna_names, all_antenna_positions);
}

void SolutionWriter::Write(
    const std::vector<std::vector<std::vector<std::complex<double>>>>&
        solutions,
    const std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>&
        constraint_solutions,
    const double start_time, const double solution_interval,
    const base::CalType mode,
    const std::vector<std::string>& used_antenna_names,
    const std::vector<base::Direction>& source_directions,
    const std::vector<std::vector<std::string>>& directions,
    const std::vector<double>& chan_freqs,
    const std::vector<double>& chan_block_freqs, const std::string& history) {
  const size_t n_times = solutions.size();
  const size_t n_channel_blocks = chan_block_freqs.size();
  const size_t n_antennas = used_antenna_names.size();
  const size_t n_directions = directions.size();
  const std::vector<std::string> direction_names =
      GetDirectionNames(directions);

  std::vector<std::pair<double, double>> h5_source_directions;
  h5_source_directions.reserve(source_directions.size());
  for (const base::Direction& direction : source_directions) {
    h5_source_directions.emplace_back(direction.ra, direction.dec);
  }
  h5parm_.AddSources(direction_names, h5_source_directions);

  std::vector<double> sol_times(n_times);
  for (size_t t = 0; t < n_times; ++t) {
    sol_times[t] = start_time + (t + 0.5) * solution_interval;
  }

  if (constraint_solutions.empty() || constraint_solutions.front().empty()) {
    // Record the actual iterands of the solver, not constraint results
    unsigned int n_pol;

    std::vector<std::string> polarizations;
    if (mode == CalType::kDiagonal || mode == CalType::kDiagonalPhase ||
        mode == CalType::kDiagonalAmplitude) {
      n_pol = 2;
      polarizations.emplace_back("XX");
      polarizations.emplace_back("YY");
    } else if (mode == CalType::kFullJones) {
      polarizations.emplace_back("XX");
      polarizations.emplace_back("XY");
      polarizations.emplace_back("YX");
      polarizations.emplace_back("YY");
      n_pol = 4;
    } else {
      n_pol = 1;
    }

    // Put solutions in a contiguous piece of memory
    std::vector<std::complex<double>> contiguous_solutions;
    const size_t n_solutions =
        n_times * n_channel_blocks * n_antennas * n_directions * n_pol;
    contiguous_solutions.reserve(n_solutions);
    for (const auto& time_solution : solutions) {
      for (const auto& block_solution : time_solution) {
        contiguous_solutions.insert(contiguous_solutions.end(),
                                    block_solution.begin(),
                                    block_solution.end());
      }
    }
    assert(contiguous_solutions.size() == n_solutions);

    std::vector<schaapcommon::h5parm::AxisInfo> axes;
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("time", n_times));
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("freq", n_channel_blocks));
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("ant", n_antennas));
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("dir", n_directions));
    if (n_pol > 1) {
      axes.emplace_back(schaapcommon::h5parm::AxisInfo("pol", n_pol));
    }

    int n_soltabs = 1;
    // For [scalar]complexgain, store two soltabs: phase and amplitude.
    if (mode == CalType::kScalar || mode == CalType::kDiagonal ||
        mode == CalType::kFullJones) {
      n_soltabs = 2;
    }
    for (int soltab_index = 0; soltab_index < n_soltabs; ++soltab_index) {
      bool store_phase = true;  // false means store amplitude.
      switch (mode) {
        case CalType::kScalar:
        case CalType::kDiagonal:
        case CalType::kFullJones:
          store_phase = soltab_index == 0;
          break;
        case CalType::kScalarPhase:
        case CalType::kDiagonalPhase:
          store_phase = true;
          break;
        case CalType::kScalarAmplitude:
        case CalType::kDiagonalAmplitude:
          store_phase = false;
          break;
        default:
          throw std::runtime_error("Constraint should have produced output");
      }

      schaapcommon::h5parm::SolTab& soltab =
          store_phase ? h5parm_.CreateSolTab("phase000", "phase", axes)
                      : h5parm_.CreateSolTab("amplitude000", "amplitude", axes);

      soltab.SetComplexValues(contiguous_solutions, std::vector<double>(),
                              !store_phase, history);
      soltab.SetAntennas(used_antenna_names);
      soltab.SetSources(GetDirectionNames(directions));
      if (n_pol > 1) {
        soltab.SetPolarizations(polarizations);
      }
      soltab.SetFreqs(chan_block_freqs);
      soltab.SetTimes(sol_times);
    }  // solnums loop
  } else {
    // Record the Constraint::Result in the H5Parm
    size_t n_constraints = constraint_solutions.front().size();

    for (size_t constraint_index = 0; constraint_index < n_constraints;
         ++constraint_index) {
      // Number of solution names, e.g. 2 for "TEC" and "ScalarPhase"
      size_t n_names = constraint_solutions[0][constraint_index].size();
      for (size_t name_index = 0; name_index < n_names; ++name_index) {
        // Get the result of the constraint solution at first time to get
        // metadata
        ddecal::Constraint::Result first_result =
            constraint_solutions[0][constraint_index][name_index];

        std::vector<hsize_t> dims(first_result.dims.size() +
                                  1);  // Add time dimension at beginning
        dims[0] = constraint_solutions.size();  // Number of times
        size_t n_solutions = dims[0];
        for (size_t i = 1; i < dims.size(); ++i) {
          dims[i] = first_result.dims[i - 1];
          n_solutions *= dims[i];
        }

        std::vector<std::string> firstaxesnames =
            common::stringtools::tokenize(first_result.axes, ",");

        std::vector<schaapcommon::h5parm::AxisInfo> axes;
        axes.emplace_back(schaapcommon::h5parm::AxisInfo(
            "time", constraint_solutions.size()));
        for (size_t axis_index = 0; axis_index < firstaxesnames.size();
             ++axis_index) {
          axes.emplace_back(schaapcommon::h5parm::AxisInfo(
              firstaxesnames[axis_index], first_result.dims[axis_index]));
        }

        // Put solutions in a contiguous piece of memory
        std::vector<double> contiguous_solutions;
        contiguous_solutions.reserve(n_solutions);
        for (size_t time = 0; time < n_times; ++time) {
          if (constraint_solutions[time].size() !=
              constraint_solutions.front().size())
            throw std::runtime_error(
                "Constraint " + std::to_string(constraint_index) +
                " did not produce a correct output at time step " +
                std::to_string(time) + ": got " +
                std::to_string(constraint_solutions[time].size()) +
                " results, expecting " +
                std::to_string(constraint_solutions.front().size()));
          const std::vector<double>& cs_vals =
              constraint_solutions[time][constraint_index][name_index].vals;
          contiguous_solutions.insert(contiguous_solutions.end(),
                                      cs_vals.begin(), cs_vals.end());
        }

        // Put solution weights in a contiguous piece of memory
        std::vector<double> weights;
        if (!constraint_solutions.front()[constraint_index][name_index]
                 .weights.empty()) {
          weights.reserve(n_solutions);
          for (unsigned int time = 0; time < solutions.size(); ++time) {
            const std::vector<double>& cs_weights =
                constraint_solutions[time][constraint_index][name_index]
                    .weights;
            weights.insert(weights.end(), cs_weights.begin(), cs_weights.end());
          }
        }

        std::string solTabName = first_result.name + "000";
        schaapcommon::h5parm::SolTab& soltab =
            h5parm_.CreateSolTab(solTabName, first_result.name, axes);
        soltab.SetValues(contiguous_solutions, weights, history);
        soltab.SetAntennas(used_antenna_names);
        soltab.SetSources(GetDirectionNames(directions));

        if (soltab.HasAxis("pol")) {
          std::vector<std::string> polarizations;
          switch (soltab.GetAxis("pol").size) {
            case 2:
              polarizations.emplace_back("XX");
              polarizations.emplace_back("YY");
              break;
            case 4:
              polarizations.emplace_back("XX");
              polarizations.emplace_back("XY");
              polarizations.emplace_back("YX");
              polarizations.emplace_back("YY");
              break;
            default:
              throw std::runtime_error(
                  "No metadata for numpolarizations = " +
                  std::to_string(soltab.GetAxis("pol").size));
          }
          soltab.SetPolarizations(polarizations);
        }

        // Set channel block frequencies. Do not use chan_block_freqs, because
        // constraint may have changed size.
        unsigned int n_channel_blocks = 1;
        if (soltab.HasAxis("freq")) {
          n_channel_blocks = soltab.GetAxis("freq").size;
        }
        std::vector<double> freqs(n_channel_blocks);
        size_t channel_index_start = 0;
        const size_t n_channels = chan_freqs.size();
        for (size_t block = 0; block != n_channel_blocks; ++block) {
          const size_t channel_index_end =
              (block + 1) * n_channels / n_channel_blocks;
          const size_t block_size = channel_index_end - channel_index_start;
          const double mean_freq =
              std::accumulate(chan_freqs.begin() + channel_index_start,
                              chan_freqs.begin() + channel_index_end, 0.0) /
              block_size;
          freqs[block] = mean_freq;
          channel_index_start = channel_index_end;
        }
        soltab.SetFreqs(freqs);

        soltab.SetTimes(sol_times);
      }
    }
  }
}

}  // namespace ddecal
}  // namespace dp3