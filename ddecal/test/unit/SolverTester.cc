// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTester.h"

#include "../../../ddecal/gain_solvers/RegularSolverBase.h"
#include "../../../ddecal/gain_solvers/BdaSolverBase.h"

#include <aocommon/matrix2x2.h>

#include <boost/make_unique.hpp>
#include <boost/test/unit_test.hpp>

#include <random>
#include <utility>

using aocommon::MC2x2;

namespace dp3 {
namespace base {
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

SolverTester::SolverTester()
    : antennas1_(),
      antennas2_(),
      input_solutions_(kNAntennas * kNDirections * 2),
      solver_solutions_(kNChannelBlocks),

      data_store_(kNRegularTimes),
      weight_store_(kNRegularTimes),
      solver_buffer_(),

      bda_data_buffer_(),
      bda_model_buffers_(),
      bda_solver_buffer_() {
  antennas1_.reserve(kNBaselines);
  antennas2_.reserve(kNBaselines);
  for (size_t ant1 = 0; ant1 != kNAntennas; ++ant1) {
    for (size_t ant2 = ant1 + 1; ant2 != kNAntennas; ++ant2) {
      antennas1_.push_back(ant1);
      antennas2_.push_back(ant2);
    }
  }
}

const SolverBuffer& SolverTester::FillData() {
  std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  std::mt19937 mt(0);
  std::vector<std::vector<DPBuffer>> model_buffers;

  for (size_t timestep = 0; timestep != kNRegularTimes; ++timestep) {
    data_buffers_.emplace_back();
    data_buffers_.back().setData(casacore::Cube<std::complex<float>>(
        kNPolarizations, kNChannels, kNBaselines, 0));
    data_buffers_.back().setWeights(
        casacore::Cube<float>(kNPolarizations, kNChannels, kNBaselines, 0));

    casacore::Cube<std::complex<float>>& time_data =
        data_buffers_.back().getData();
    casacore::Cube<float> time_weights = data_buffers_.back().getWeights();

    model_buffers.emplace_back();
    std::vector<DPBuffer>& model_time_buffers = model_buffers.back();

    for (size_t d = 0; d != kNDirections; ++d) {
      model_time_buffers.emplace_back();
      model_time_buffers.back().setData(casacore::Cube<std::complex<float>>(
          kNPolarizations, kNChannels, kNBaselines));
      std::complex<float>* this_direction =
          model_time_buffers.back().getData().data();

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
            MC2x2 val(&model_time_buffers[d].getData()(0, ch, baseline_index));
            MC2x2 left(input_solutions_[(a1 * kNDirections + d) * 2 + 0], 0.0,
                       0.0, input_solutions_[(a1 * kNDirections + d) * 2 + 1]);
            MC2x2 right(input_solutions_[(a2 * kNDirections + d) * 2 + 0], 0.0,
                        0.0, input_solutions_[(a2 * kNDirections + d) * 2 + 1]);
            MC2x2 res;
            MC2x2::ATimesB(res, left, val);
            // Use 'left' as scratch for the result.
            MC2x2::ATimesHermB(left, res, right);
            perturbed_model += left;
          }
          for (size_t p = 0; p != 4; ++p) {
            time_data(p, ch, baseline_index) = perturbed_model[p];
            time_weights(p, ch, baseline_index) = 1.0;
          }
        }
        ++baseline_index;
      }
    }
  }

  solver_buffer_.AssignAndWeight(data_buffers_, std::move(model_buffers));

  return solver_buffer_;
}

const BDASolverBuffer& SolverTester::FillBDAData() {
  assert(bda_model_buffers_.empty());

  std::uniform_real_distribution<float> uniform_data(-1.0, 1.0);
  std::mt19937 mt(0);

  // Reusable vector for measurement data and model data.
  // Declaring it outside the loops below allows memory reuse.
  std::vector<std::complex<float>> data;

  // Initialize the data buffers. The solvers only need the data field.
  BDABuffer::Fields bda_fields(false);
  bda_fields.data = true;
  bda_data_buffer_ = boost::make_unique<BDABuffer>(kBDABufferSize, bda_fields);
  bda_model_buffers_.reserve(kNDirections);
  for (size_t dir = 0; dir < kNDirections; ++dir) {
    bda_model_buffers_.push_back(
        boost::make_unique<BDABuffer>(kBDABufferSize, bda_fields));
  }

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
        // Generate model data for all directions.
        for (size_t dir = 0; dir < kNDirections; ++dir) {
          data.clear();
          for (size_t ch = 0; ch < n_averaged_channels; ++ch) {
            data.emplace_back(uniform_data(mt), uniform_data(mt));
            data.emplace_back(uniform_data(mt) * 0.1, uniform_data(mt) * 0.1);
            data.emplace_back(uniform_data(mt) * 0.1, uniform_data(mt) * 0.1);
            data.emplace_back(uniform_data(mt) * 1.5, uniform_data(mt) * 1.5);
          }
          BOOST_REQUIRE(bda_model_buffers_[dir]->AddRow(
              time, averaging_factors.first, n_averaged_times, bl,
              n_averaged_channels, kNPolarizations, data.data()));
        }

        // Generate measurement data using the model data.
        const size_t ant1 = antennas1_[bl];
        const size_t ant2 = antennas2_[bl];
        data.clear();
        for (size_t ch = 0; ch < n_averaged_channels; ++ch) {
          MC2x2 perturbed_model = MC2x2::Zero();
          for (size_t dir = 0; dir < kNDirections; ++dir) {
            const size_t ant1_index = (ant1 * kNDirections + dir) * 2;
            const size_t ant2_index = (ant2 * kNDirections + dir) * 2;

            MC2x2 val(&bda_model_buffers_[dir]->GetRows().back().data[ch * 4]);
            MC2x2 left(input_solutions_[ant1_index], 0.0, 0.0,
                       input_solutions_[ant1_index + 1]);
            MC2x2 right(input_solutions_[ant2_index], 0.0, 0.0,
                        input_solutions_[ant2_index + 1]);
            MC2x2 res;
            MC2x2::ATimesB(res, left, val);
            // Use 'left' as scratch for the result.
            MC2x2::ATimesHermB(left, res, right);
            perturbed_model += left;
          }
          for (size_t p = 0; p != kNPolarizations; ++p) {
            data.emplace_back(perturbed_model[p]);
          }
        }
        BOOST_REQUIRE(bda_data_buffer_->AddRow(
            time, averaging_factors.first, n_averaged_times, bl,
            n_averaged_channels, kNPolarizations, data.data()));
      }
    }
  }

  bda_solver_buffer_.SetDataRows(*bda_data_buffer_);
  bda_solver_buffer_.SetModelDataRows(bda_model_buffers_);

  return bda_solver_buffer_;
}

void SolverTester::InitializeSolver(
    dp3::base::RegularSolverBase& solver) const {
  solver.SetMaxIterations(kMaxIterations);
  solver.SetAccuracy(kAccuracy);
  solver.SetStepSize(kStepSize);
  solver.SetNThreads(kNThreads);
  solver.SetPhaseOnly(kPhaseOnly);
  solver.Initialize(kNAntennas, kNDirections, kNChannels, kNChannelBlocks,
                    antennas1_, antennas2_);
}

void SolverTester::InitializeSolver(dp3::base::BdaSolverBase& solver) const {
  solver.SetMaxIterations(kMaxIterations);
  solver.SetAccuracy(kAccuracy);
  solver.SetStepSize(kStepSize);
  solver.SetNThreads(kNThreads);
  solver.SetPhaseOnly(kPhaseOnly);
  solver.Initialize(kNAntennas, kNDirections, kNChannelBlocks);
}

void SolverTester::SetScalarSolutions() {
  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t ant = 0; ant != kNAntennas; ++ant) {
    for (size_t dir = 0; dir != kNDirections; ++dir) {
      const size_t index = (ant * kNDirections + dir) * 2;
      input_solutions_[index] = {uniform_sols(mt), uniform_sols(mt)};
      if (dir != 0) input_solutions_[index] *= 0.5f;
      input_solutions_[index + 1] = input_solutions_[index];
    }
  }

  // Initialize unit-matrices as initial solver solution values.
  for (auto& vec : solver_solutions_) {
    vec.assign(kNDirections * kNAntennas, 1.0);
  }
}

void SolverTester::SetDiagonalSolutions() {
  std::mt19937 mt;
  std::uniform_real_distribution<float> uniform_sols(1.0, 2.0);
  for (size_t ant = 0; ant != kNAntennas; ++ant) {
    for (size_t pol = 0; pol != 2; ++pol) {
      for (size_t dir = 0; dir != kNDirections; ++dir) {
        const size_t index = (ant * kNDirections + dir) * 2 + pol;
        input_solutions_[index] = {uniform_sols(mt), uniform_sols(mt)};
        if (dir != 0) input_solutions_[index] *= 0.5f;
      }
    }
  }

  // Initialize unit-matrices as initial solver solution values.
  for (auto& vec : solver_solutions_) {
    vec.assign(kNDirections * kNAntennas * 2, 1.0);
  }
}

void SolverTester::CheckScalarResults(double tolerance) {
  for (size_t ch = 0; ch != kNChannelBlocks; ++ch) {
    for (size_t ant = 0; ant != kNAntennas; ++ant) {
      for (size_t d = 0; d != kNDirections; ++d) {
        const std::complex<double> sol0 = solver_solutions_[ch][d];
        const std::complex<double> inp0 = input_solutions_[d * 2];

        const std::complex<double> sol =
            solver_solutions_[ch][d + ant * kNDirections];
        const std::complex<double> inp =
            input_solutions_[(d + ant * kNDirections) * 2];

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
      for (size_t d = 0; d != kNDirections; ++d) {
        std::complex<double> solX0 = solver_solutions_[ch][d * 2];
        std::complex<double> solY0 = solver_solutions_[ch][d * 2 + 1];
        std::complex<double> inpX0 = input_solutions_[d * 2];
        std::complex<double> inpY0 = input_solutions_[d * 2 + 1];

        std::complex<double> solX =
            solver_solutions_[ch][(d + ant * kNDirections) * 2];
        std::complex<double> solY =
            solver_solutions_[ch][(d + ant * kNDirections) * 2 + 1];
        std::complex<double> inpX =
            input_solutions_[(d + ant * kNDirections) * 2];
        std::complex<double> inpY =
            input_solutions_[(d + ant * kNDirections) * 2 + 1];

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
}  // namespace base
}  // namespace dp3
