// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolutionWriter.h"

#include "../common/StreamUtil.h"
#include "../common/StringTools.h"

#include <cassert>
#include <numeric>
#include <algorithm>

#include "SolutionResampler.h"

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

/// Allows constructing h5parm_ directly in the SolutionWriter constructor
/// and catching any exception. Also avoids using the deprecated implicitly
//  declared assignment operator of H5Parm, by using the constructor instead.
schaapcommon::h5parm::H5Parm ConstructH5Parm(const std::string& filename) {
  try {
    return schaapcommon::h5parm::H5Parm(filename, true);
  } catch (H5::FileIException& e) {
    throw std::runtime_error("Error opening '" + filename +
                             "' for writing: " + e.getCDetailMsg());
  }
}

std::vector<double> GetSolutionTimes(size_t n_times, double start_time,
                                     double end_time, double solution_interval,
                                     double ms_timestep_duration) {
  // Filling and trimming the time axis to have the same time range of the
  // input data. sol_times keeps the center of the time interval.
  // If an integral number of solution intervals fit in the ms, then the first
  // interval that is entirely outside the ms has a start time that equals the
  // end time of the ms. Whether this inequality holds depends on the rounding.
  // The extra half ms_interval margin makes sure that the outlying interval
  // is always excluded.
  std::vector<double> sol_times(n_times);
  for (size_t t = 0; t < n_times; ++t) {
    sol_times[t] = start_time + (t + 0.5) * solution_interval;
    if ((start_time + t * solution_interval) >
        (end_time - 0.5 * ms_timestep_duration)) {
      sol_times.resize(t);
      break;
    }
  }
  return sol_times;
}

}  // namespace

namespace dp3 {
namespace ddecal {

SolutionWriter::SolutionWriter(const std::string& filename)
    : h5parm_(ConstructH5Parm(filename)) {}

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
    double start_time, double end_time, double ms_timestep_duration,
    size_t n_interval_timesteps,
    const std::vector<size_t>& solutions_per_direction, base::CalType mode,
    const std::vector<std::string>& used_antenna_names,
    const std::vector<base::Direction>& source_directions,
    const std::vector<std::vector<std::string>>& directions,
    const std::vector<double>& chan_freqs,
    const std::vector<double>& chan_block_freqs, const std::string& history) {
  const size_t n_directions = solutions_per_direction.size();
  const size_t n_sub_solutions = std::accumulate(
      solutions_per_direction.begin(), solutions_per_direction.end(), 0u);
  if (n_sub_solutions == n_directions) {
    const double solution_duration =
        ms_timestep_duration * n_interval_timesteps;
    WriteDirect(solutions, constraint_solutions, start_time, end_time,
                ms_timestep_duration, solution_duration, mode,
                used_antenna_names, source_directions, directions, chan_freqs,
                chan_block_freqs, history);
  } else {
    const size_t n_antennas = used_antenna_names.size();
    const SolutionResampler resampler(solutions_per_direction, n_antennas,
                                      GetNPolarizations(mode),
                                      n_interval_timesteps);
    const size_t n_resampled_interval_timesteps =
        n_interval_timesteps / resampler.GetNrSubSteps();

    std::vector<std::vector<std::vector<std::complex<double>>>>
        upsampled_solutions = resampler.Upsample(solutions);
    WriteDirect(upsampled_solutions, constraint_solutions, start_time, end_time,
                ms_timestep_duration,
                ms_timestep_duration * n_resampled_interval_timesteps, mode,
                used_antenna_names, source_directions, directions, chan_freqs,
                chan_block_freqs, history);
  }
}

void SolutionWriter::WriteDirect(
    const std::vector<std::vector<std::vector<std::complex<double>>>>&
        solutions,
    const std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>&
        constraint_solutions,
    const double start_time, const double end_time,
    const double ms_timestep_duration, const double solution_interval,
    const base::CalType mode,
    const std::vector<std::string>& used_antenna_names,
    const std::vector<base::Direction>& source_directions,
    const std::vector<std::vector<std::string>>& directions,
    const std::vector<double>& chan_freqs,
    const std::vector<double>& chan_block_freqs, const std::string& history) {
  const std::vector<std::string> direction_names =
      GetDirectionNames(directions);

  std::vector<std::pair<double, double>> h5_source_directions;
  h5_source_directions.reserve(source_directions.size());
  for (const base::Direction& direction : source_directions) {
    h5_source_directions.emplace_back(direction.ra, direction.dec);
  }
  h5parm_.AddSources(direction_names, h5_source_directions);

  const std::vector<double> solution_times =
      GetSolutionTimes(solutions.size(), start_time, end_time,
                       solution_interval, ms_timestep_duration);
  if (constraint_solutions.empty() || constraint_solutions.front().empty()) {
    WriteSolverResults(solutions, mode, used_antenna_names, source_directions,
                       directions, chan_freqs, chan_block_freqs, history,
                       solution_times);
  } else {
    WriteConstraintResults(constraint_solutions, mode, used_antenna_names,
                           source_directions, directions, chan_freqs,
                           chan_block_freqs, history, solution_times);
  }
}

void SolutionWriter::WriteSolverResults(
    const std::vector<std::vector<std::vector<std::complex<double>>>>&
        solutions,
    base::CalType mode, const std::vector<std::string>& used_antenna_names,
    const std::vector<base::Direction>& source_directions,
    const std::vector<std::vector<std::string>>& directions,
    const std::vector<double>& chan_freqs,
    const std::vector<double>& chan_block_freqs, const std::string& history,
    const std::vector<double>& solution_times) {
  const size_t n_channel_blocks = chan_block_freqs.size();
  const size_t n_antennas = used_antenna_names.size();
  const size_t n_directions = directions.size();
  const size_t n_times = solution_times.size();

  unsigned int n_pol;
  std::vector<std::string> polarizations;
  if (mode == CalType::kDiagonal || mode == CalType::kDiagonalPhase ||
      mode == CalType::kDiagonalAmplitude) {
    polarizations = {"XX", "YY"};
    n_pol = 2;
  } else if (mode == CalType::kFullJones) {
    polarizations = {"XX", "XY", "YX", "YY"};
    n_pol = 4;
  } else {
    n_pol = 1;
  }

  // Put solutions in a contiguous piece of memory
  std::vector<std::complex<double>> contiguous_solutions;
  const size_t n_solutions =
      n_times * n_channel_blocks * n_antennas * n_directions * n_pol;
  contiguous_solutions.reserve(n_solutions);
  for (size_t i = 0; i < n_times; i++) {
    for (const auto& block_solution : solutions[i]) {
      contiguous_solutions.insert(contiguous_solutions.end(),
                                  block_solution.begin(), block_solution.end());
    }
  }
  assert(contiguous_solutions.size() == n_solutions);

  std::vector<schaapcommon::h5parm::AxisInfo> axes;
  axes.push_back({"time", static_cast<unsigned>(n_times)});
  axes.push_back({"freq", static_cast<unsigned>(n_channel_blocks)});
  axes.push_back({"ant", static_cast<unsigned>(n_antennas)});
  axes.push_back({"dir", static_cast<unsigned>(n_directions)});
  if (n_pol > 1) {
    axes.push_back({"pol", static_cast<unsigned>(n_pol)});
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
    soltab.SetTimes(solution_times);
  }  // solnums loop
}

void SolutionWriter::WriteConstraintResults(
    const std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>&
        constraint_solutions,
    base::CalType mode, const std::vector<std::string>& used_antenna_names,
    const std::vector<base::Direction>& source_directions,
    const std::vector<std::vector<std::string>>& directions,
    const std::vector<double>& chan_freqs,
    const std::vector<double>& chan_block_freqs, const std::string& history,
    const std::vector<double>& solution_times) {
  const size_t n_constraints = constraint_solutions.front().size();
  const size_t n_times = solution_times.size();

  for (size_t constraint_index = 0; constraint_index < n_constraints;
       ++constraint_index) {
    // Number of solution names, e.g. 2 for "TEC" and "ScalarPhase"
    const size_t n_names = constraint_solutions[0][constraint_index].size();
    for (size_t name_index = 0; name_index < n_names; ++name_index) {
      // Get the result of the constraint solution at first time to get
      // metadata
      const ddecal::Constraint::Result& first_result =
          constraint_solutions[0][constraint_index][name_index];

      std::vector<hsize_t> dims(first_result.dims.size() +
                                1);           // Add time dimension at beginning
      dims[0] = constraint_solutions.size();  // Number of times
      size_t n_solutions = dims[0];
      for (size_t i = 1; i < dims.size(); ++i) {
        dims[i] = first_result.dims[i - 1];
        n_solutions *= dims[i];
      }

      std::vector<std::string> first_axes_names =
          common::stringtools::tokenize(first_result.axes, ",");

      std::vector<schaapcommon::h5parm::AxisInfo> axes;
      axes.push_back(
          {"time", static_cast<unsigned>(constraint_solutions.size())});
      for (size_t axis_index = 0; axis_index < first_axes_names.size();
           ++axis_index) {
        axes.push_back({first_axes_names[axis_index],
                        static_cast<unsigned>(first_result.dims[axis_index])});
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
        contiguous_solutions.insert(contiguous_solutions.end(), cs_vals.begin(),
                                    cs_vals.end());
      }

      // Put solution weights in a contiguous piece of memory
      std::vector<double> weights;
      if (!constraint_solutions.front()[constraint_index][name_index]
               .weights.empty()) {
        weights.reserve(n_solutions);
        for (unsigned int time = 0; time < constraint_solutions.size();
             ++time) {
          const std::vector<double>& cs_weights =
              constraint_solutions[time][constraint_index][name_index].weights;
          weights.insert(weights.end(), cs_weights.begin(), cs_weights.end());
        }
      }

      const std::string sol_tab_name = first_result.name + "000";
      schaapcommon::h5parm::SolTab& soltab =
          h5parm_.CreateSolTab(sol_tab_name, first_result.name, axes);
      soltab.SetValues(contiguous_solutions, weights, history);
      soltab.SetAntennas(used_antenna_names);
      soltab.SetSources(GetDirectionNames(directions));

      if (soltab.HasAxis("pol")) {
        std::vector<std::string> polarizations;
        switch (soltab.GetAxis("pol").size) {
          case 2:
            polarizations = {"XX", "YY"};
            break;
          case 4:
            polarizations = {"XX", "XY", "YX", "YY"};
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

      soltab.SetTimes(solution_times);
    }
  }
}

}  // namespace ddecal
}  // namespace dp3
