// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTester.h"

#include "../../../ddecal/gain_solvers/SolverBase.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/threadpool.h>

#include <boost/test/unit_test.hpp>

#include <random>
#include <utility>

using aocommon::MC2x2;
using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;

namespace dp3 {
namespace ddecal {
namespace test {

namespace {

std::pair<size_t, size_t> GetAveragingFactors(size_t baseline_index) {
  size_t time_averaging_factor = 1;
  size_t channel_averaging_factor = 1;
  // Divide the baselines in blocks of 8 baselines, where each baseline
  // in the block uses a different averaging strategy.
  switch (baseline_index % 8) {
    case 0:  // Do not use averaging / keep using the factor 1.
      break;
    case 1:
      time_averaging_factor = 2;
      break;
    case 2:
      time_averaging_factor = 3;
      break;
    case 3:
      channel_averaging_factor = 2;
      break;
    case 4:
      channel_averaging_factor = 3;
      break;
    case 5:
      time_averaging_factor = 2;
      channel_averaging_factor = 3;
      break;
    case 6:
      time_averaging_factor = 3;
      channel_averaging_factor = 2;
      break;
    case 7:
      time_averaging_factor = 4;
      channel_averaging_factor = 4;
      break;
  }
  return std::make_pair(time_averaging_factor, channel_averaging_factor);
}

}  // namespace

constexpr size_t SolverTester::kDDSolutionsPerDirection[];

SolverTester::SolverTester()
    : antennas1_(),
      antennas2_(),
      solver_solutions_(kNChannelBlocks),
      n_solutions_(0),

      // Set the solution interval so that all data is in the first interval.
      bda_solver_buffer_(-1.0, kNBDATimes + 1.0, kNBaselines) {
  antennas1_.reserve(kNBaselines);
  antennas2_.reserve(kNBaselines);
  for (size_t ant1 = 0; ant1 != kNAntennas; ++ant1) {
    for (size_t ant2 = ant1 + 1; ant2 != kNAntennas; ++ant2) {
      antennas1_.push_back(ant1);
      antennas2_.push_back(ant2);
    }
  }

  // Use 4 threads in the solver tests. This value corresponds to the
  // PROCESSORS setting in CMakeLists.txt for the solver tests.
  aocommon::ThreadPool::GetInstance().SetNThreads(4);
}

std::vector<std::string> SolverTester::CreateDirectionNames() {
  std::vector<std::string> keys;
  for (size_t i = 0; i < kNDirections; ++i) {
    keys.emplace_back(std::to_string(i));
  }
  return keys;
}

std::vector<dp3::base::DPBuffer> SolverTester::FillDdIntervalData() {
  std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  std::mt19937 mt(0);

  std::vector<std::vector<size_t>> solution_indices(kNDirections);
  size_t first_solution_index = 0;
  for (size_t d = 0; d != kNDirections; ++d) {
    solution_indices[d].resize(kNRegularTimes);
    for (size_t t = 0; t != kNRegularTimes; ++t) {
      solution_indices[d][t] =
          first_solution_index +
          t * n_solutions_per_direction_[d] / kNRegularTimes;
    }
    first_solution_index += n_solutions_per_direction_[d];
  }

  const std::vector<std::string> direction_names = CreateDirectionNames();

  std::vector<std::unique_ptr<dp3::base::DPBuffer>> unweighted_buffers;
  std::vector<dp3::base::DPBuffer> weighted_buffers;
  for (size_t timestep = 0; timestep != kNRegularTimes; ++timestep) {
    unweighted_buffers.emplace_back(std::make_unique<base::DPBuffer>());
    weighted_buffers.emplace_back();

    unweighted_buffers.back()->GetData().resize(kShape);
    unweighted_buffers.back()->GetWeights().resize(kShape);
    unweighted_buffers.back()->GetFlags().resize(kShape);
    unweighted_buffers.back()->GetWeights().fill(1.0f);
    unweighted_buffers.back()->GetFlags().fill(false);

    base::DPBuffer::DataType& time_data = unweighted_buffers.back()->GetData();
    time_data.fill(std::complex<float>{0.0f, 0.0f});

    for (size_t d = 0; d != kNDirections; ++d) {
      const std::string& name = direction_names[d];
      unweighted_buffers.back()->AddData(name);
      std::complex<float>* this_direction =
          unweighted_buffers.back()->GetData(name).data();

      for (size_t bl = 0; bl != kNBaselines; ++bl) {
        for (size_t ch = 0; ch != kNChannels; ++ch) {
          const size_t matrix_index = (bl * kNChannels + ch) * 4;
          this_direction[matrix_index + 0] =
              std::complex<float>(uniform_data(mt), uniform_data(mt));
          this_direction[matrix_index + 1] =
              std::complex<float>(uniform_data(mt), uniform_data(mt)) * 0.1f;
          this_direction[matrix_index + 2] =
              std::complex<float>(uniform_data(mt), uniform_data(mt)) * 0.1f;
          this_direction[matrix_index + 3] =
              std::complex<float>(uniform_data(mt), uniform_data(mt)) * 1.5f;
        }
      }
    }
    size_t baseline_index = 0;
    for (size_t a1 = 0; a1 != kNAntennas; ++a1) {
      for (size_t a2 = a1 + 1; a2 != kNAntennas; ++a2) {
        for (size_t ch = 0; ch != kNChannels; ++ch) {
          MC2x2 perturbed_model = MC2x2::Zero();
          for (size_t d = 0; d != kNDirections; ++d) {
            const std::string& name = direction_names[d];
            const size_t solution_index = solution_indices[d][timestep];
            const MC2x2 val(&unweighted_buffers.back()->GetData(name)(
                baseline_index, ch, 0));
            MC2x2 left(
                input_solutions_[(a1 * n_solutions_ + solution_index) * 2 + 0],
                0.0, 0.0,
                input_solutions_[(a1 * n_solutions_ + solution_index) * 2 + 1]);
            const MC2x2 right(
                input_solutions_[(a2 * n_solutions_ + solution_index) * 2 + 0],
                0.0, 0.0,
                input_solutions_[(a2 * n_solutions_ + solution_index) * 2 + 1]);
            perturbed_model += left * val.MultiplyHerm(right);
          }
          for (size_t p = 0; p != 4; ++p) {
            time_data(baseline_index, ch, p) = perturbed_model.Get(p);
          }
        }
        ++baseline_index;
      }
    }
  }

  dp3::ddecal::AssignAndWeight(unweighted_buffers, direction_names,
                               weighted_buffers, false, false);
  return weighted_buffers;
}

const BdaSolverBuffer& SolverTester::FillBDAData() {
  using dp3::common::Fields;

  std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  std::mt19937 mt(0);

  const std::vector<float> weights(kNChannels * kNPolarizations, 1.0f);

  // For each block of 8 baselines, the averaged data size corresponds to
  // 1 + 1/2 + 1/2 + 1/3 + 1/3 + 1/6 + 1/6 + 1/16 = 3 + 1/16
  // non-averaged baselines. -> Divide the non-averaged data by two.
  const size_t buffer_size =
      kNBDATimes * kNBaselines * kNChannels * kNPolarizations / 2;

  // Initialize the data buffers.
  const Fields fields = Fields(Fields::Single::kData) |
                        Fields(Fields::Single::kWeights) |
                        Fields(Fields::Single::kFlags);

  const std::vector<std::string> direction_names = CreateDirectionNames();

  auto buffer = std::make_unique<BdaBuffer>(buffer_size, fields);
  for (const std::string& name : direction_names) buffer->AddData(name);

  // Do the outer loop over time, since the BDA rows should be ordered.
  for (size_t time = 0; time < kNBDATimes; ++time) {
    for (size_t bl = 0; bl < kNBaselines; ++bl) {
      const std::pair<size_t, size_t> averaging_factors =
          GetAveragingFactors(bl);
      const size_t n_averaged_times = averaging_factors.first;
      const size_t n_averaged_channels =
          (kNChannels + averaging_factors.second - 1) /
          averaging_factors.second;

      if ((time % n_averaged_times) == 0) {
        const std::size_t row = buffer->GetRows().size();

        BOOST_REQUIRE(buffer->AddRow(
            time, n_averaged_times * kBdaUnitTimeInterval,
            n_averaged_times * kBdaUnitTimeInterval, bl, n_averaged_channels,
            kNPolarizations, nullptr, nullptr, weights.data()));

        // Generate model data for all directions.
        for (const std::string& name : direction_names) {
          std::complex<float>* model_data = buffer->GetData(row, name);
          for (size_t ch = 0; ch < n_averaged_channels; ++ch) {
            model_data[0] =
                std::complex<float>(uniform_data(mt), uniform_data(mt));
            model_data[1] = std::complex<float>(uniform_data(mt) * 0.1,
                                                uniform_data(mt) * 0.1);
            model_data[2] = std::complex<float>(uniform_data(mt) * 0.1,
                                                uniform_data(mt) * 0.1);
            model_data[3] = std::complex<float>(uniform_data(mt) * 1.5,
                                                uniform_data(mt) * 1.5);
            model_data += kNPolarizations;
          }
        }

        // Generate measurement data using the model data.
        const size_t ant1 = antennas1_[bl];
        const size_t ant2 = antennas2_[bl];
        for (size_t ch = 0; ch < n_averaged_channels; ++ch) {
          MC2x2 perturbed_model = MC2x2::Zero();

          for (size_t dir = 0; dir < kNDirections; ++dir) {
            const size_t ant1_index = (ant1 * kNDirections + dir) * 2;
            const size_t ant2_index = (ant2 * kNDirections + dir) * 2;

            MC2x2 val(&buffer->GetData(row, direction_names[dir])[ch * 4]);
            MC2x2 left(input_solutions_[ant1_index], 0.0, 0.0,
                       input_solutions_[ant1_index + 1]);
            MC2x2 right(input_solutions_[ant2_index], 0.0, 0.0,
                        input_solutions_[ant2_index + 1]);
            perturbed_model += left * val.MultiplyHerm(right);
          }
          perturbed_model.AssignTo(&buffer->GetData(row)[ch * 4]);
        }
      }
    }
  }

  const size_t n_bda_buffer_rows = buffer->GetRows().size();
  bda_solver_buffer_.AppendAndWeight(std::move(buffer), direction_names, false);
  BOOST_REQUIRE_EQUAL(bda_solver_buffer_.GetIntervalRows().size(),
                      n_bda_buffer_rows);

  return bda_solver_buffer_;
}

void SolverTester::InitializeSolver(SolverBase& solver) const {
  InitializeSolverSettings(solver);
  solver.Initialize(kNAntennas, n_solutions_per_direction_, kNChannelBlocks);
}

void SolverTester::InitializeSolverSettings(SolverBase& solver) const {
  std::vector<SolverBase*> solvers = solver.ConstraintSolvers();
  // If this is a hybrid solver, it won't list itself, so add it manually
  if (std::find(solvers.begin(), solvers.end(), &solver) == solvers.end())
    solvers.push_back(&solver);
  for (SolverBase* s : solvers) {
    s->SetAccuracy(kAccuracy);
    s->SetStepSize(kStepSize);
    s->SetPhaseOnly(kPhaseOnly);
  }
  // Because the maximum number of iterations need to be set before adding
  // solvers to a hybrid solver, we only set the max of the hybrid solver:
  solver.SetMaxIterations(kMaxIterations);
}

void SolverTester::InitializeNSolutions(bool use_dd_intervals) {
  if (use_dd_intervals) {
    n_solutions_per_direction_.assign(std::begin(kDDSolutionsPerDirection),
                                      std::end(kDDSolutionsPerDirection));
    n_solutions_ = std::accumulate(std::begin(kDDSolutionsPerDirection),
                                   std::end(kDDSolutionsPerDirection), 0.0);
  } else {
    n_solutions_per_direction_.assign(kNDirections, 1);
    n_solutions_ = kNDirections;
  }
}

void SolverTester::SetScalarSolutions(bool use_dd_intervals) {
  InitializeNSolutions(use_dd_intervals);
  input_solutions_.resize(kNAntennas * n_solutions_ * 2);
  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t ant = 0; ant != kNAntennas; ++ant) {
    for (size_t s = 0; s != n_solutions_; ++s) {
      const size_t index = (ant * n_solutions_ + s) * 2;
      input_solutions_[index] = {uniform_sols(mt), uniform_sols(mt)};
      if (s != 0) input_solutions_[index] *= 0.5f;
      input_solutions_[index + 1] = input_solutions_[index];
    }
  }

  // Initialize unit-matrices as initial solver solution values.
  for (auto& vec : solver_solutions_) {
    vec.assign(n_solutions_ * kNAntennas, 1.0);
  }
}

void SolverTester::SetDiagonalSolutions(bool use_dd_intervals) {
  InitializeNSolutions(use_dd_intervals);
  input_solutions_.resize(kNAntennas * n_solutions_ * 2);
  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t ant = 0; ant != kNAntennas; ++ant) {
    for (size_t pol = 0; pol != 2; ++pol) {
      for (size_t s = 0; s != n_solutions_; ++s) {
        const size_t index = (ant * n_solutions_ + s) * 2 + pol;
        input_solutions_[index] = {uniform_sols(mt), uniform_sols(mt)};
        if (s != 0) input_solutions_[index] *= 0.5f;
      }
    }
  }

  // Initialize unit-matrices as initial solver solution values.
  for (auto& vec : solver_solutions_) {
    vec.assign(n_solutions_ * kNAntennas * 2, 1.0);
  }
}

void SolverTester::CheckScalarResults(double tolerance) {
  for (size_t ch = 0; ch != kNChannelBlocks; ++ch) {
    for (size_t ant = 0; ant != kNAntennas; ++ant) {
      for (size_t s = 0; s != n_solutions_; ++s) {
        const std::complex<double> sol0 = solver_solutions_[ch][s];
        const std::complex<double> inp0 = input_solutions_[s * 2];

        const std::complex<double> sol =
            solver_solutions_[ch][s + ant * n_solutions_];
        const std::complex<double> inp =
            input_solutions_[(s + ant * n_solutions_) * 2];

        // Compare the squared quantities, because the phase has an ambiguity
        BOOST_CHECK_CLOSE(std::norm(sol), std::norm(inp), tolerance);

        // Reference to antenna0 to check if relative phase is correct
        BOOST_CHECK_LT(std::abs((sol * std::conj(sol0)).real() -
                                (inp * std::conj(inp0)).real()),
                       tolerance);
        BOOST_CHECK_LT(std::abs((sol * std::conj(sol0)).imag() -
                                (inp * std::conj(inp0)).imag()),
                       tolerance);
      }
    }
  }
}

void SolverTester::CheckDiagonalResults(double tolerance) {
  for (size_t ch = 0; ch != kNChannelBlocks; ++ch) {
    for (size_t ant = 0; ant != kNAntennas; ++ant) {
      for (size_t s = 0; s != n_solutions_; ++s) {
        const std::complex<double> solX0 = solver_solutions_[ch][s * 2];
        const std::complex<double> solY0 = solver_solutions_[ch][s * 2 + 1];
        const std::complex<double> inpX0 = input_solutions_[s * 2];
        const std::complex<double> inpY0 = input_solutions_[s * 2 + 1];

        const std::complex<double> solX =
            solver_solutions_[ch][(s + ant * n_solutions_) * 2];
        const std::complex<double> solY =
            solver_solutions_[ch][(s + ant * n_solutions_) * 2 + 1];
        const std::complex<double> inpX =
            input_solutions_[(s + ant * n_solutions_) * 2];
        const std::complex<double> inpY =
            input_solutions_[(s + ant * n_solutions_) * 2 + 1];

        // Compare the squared quantities, because the phase has an ambiguity
        BOOST_CHECK_CLOSE(std::norm(solX), std::norm(inpX), tolerance);
        BOOST_CHECK_CLOSE(std::norm(solY), std::norm(inpY), tolerance);

        // Reference to antenna0 to check if relative phase is correct
        BOOST_CHECK_CLOSE((solX * std::conj(solX0)).real(),
                          (inpX * std::conj(inpX0)).real(), tolerance);
        BOOST_CHECK_CLOSE((solY * std::conj(solY0)).real(),
                          (inpY * std::conj(inpY0)).real(), tolerance);
      }
    }
  }
}

}  // namespace test
}  // namespace ddecal
}  // namespace dp3
