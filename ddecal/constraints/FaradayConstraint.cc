#include "FaradayConstraint.h"

#include <cassert>
#include <cmath>

#include <aocommon/constants.h>
#include <aocommon/staticfor.h>

#include "RotationAndDiagonalConstraint.h"
#include "RotationConstraint.h"

namespace dp3::ddecal {

using common::phase_fitting::FitSample;

void FaradayConstraint::Initialize(
    size_t nAntennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, solutions_per_direction, frequencies);

  frequencies_ = frequencies;

  std::vector<FitSample> grid;
  // The x values of the fitter must be ordered lowest to highest, so add
  // them in reverse order of frequency.
  for (auto iter = frequencies.rbegin(); iter != frequencies.rend(); ++iter) {
    const double wavelength = aocommon::kSpeedOfLight / *iter;
    grid.emplace_back(FitSample{wavelength * wavelength, 0.0, 0.0});
  }
  fit_range_ = common::phase_fitting::GetRange(grid);
  if (max_rotation_value_) {
    assert(*max_rotation_value_ > 0.0);
    if (fit_range_.wrap_step * fit_range_.max_wraps > *max_rotation_value_ &&
        fit_range_.wrap_step != 0.0) {
      fit_range_.max_wraps =
          std::ceil(*max_rotation_value_ / fit_range_.wrap_step);
    }
  }

  Result& rm_result = results_.emplace_back();
  rm_result.vals.resize(NAntennas() * NSubSolutions());
  rm_result.axes = "ant,dir";
  rm_result.dims.resize(2);
  rm_result.dims[0] = NAntennas();
  rm_result.dims[1] = NSubSolutions();
  rm_result.name = "rotationmeasure";

  std::vector<Result> diagonal_result = MakeDiagonalResults(
      NAntennas(), NSubSolutions(), NChannelBlocks(), diagonal_solution_type_);
  for (Result& result : diagonal_result)
    results_.emplace_back(std::move(result));
}

void FaradayConstraint::SetWeights(const std::vector<double>& weights) {
  assert(weights.size() == NAntennas() * NChannelBlocks());
  sub_solution_weights_.clear();
  for (size_t i = 0; i != NSubSolutions(); ++i) {
    sub_solution_weights_.emplace_back(weights);
  }
  for (Result& result : results_) result.weights.clear();
  Result& rm_result = results_.front();
  const size_t n_polarizations = GetNPolarizations(diagonal_solution_type_);
  for (size_t antenna = 0; antenna != NAntennas(); ++antenna) {
    const double* antenna_weights = &weights[antenna * NChannelBlocks()];
    double weight_sum = 0.0;
    for (size_t channel = 0; channel != NChannelBlocks(); ++channel) {
      weight_sum += antenna_weights[channel];
    }
    for (size_t sub_solution = 0; sub_solution != NSubSolutions();
         ++sub_solution) {
      rm_result.weights.emplace_back(weight_sum);
      for (size_t result_index = 1; result_index < results_.size();
           ++result_index) {
        for (size_t channel = 0; channel != NChannelBlocks(); ++channel) {
          for (size_t p = 0; p != n_polarizations; ++p) {
            results_[result_index].weights.emplace_back(
                antenna_weights[channel]);
          }
        }
      }
    }
  }
}

void FaradayConstraint::SetSubSolutionWeights(
    const std::vector<std::vector<double>>& sub_solution_weights) {
  sub_solution_weights_ = sub_solution_weights;
  for (Result& result : results_) result.weights.clear();
  Result& rm_result = results_.front();
  for (size_t antenna = 0; antenna != NAntennas(); ++antenna) {
    for (const std::vector<double>& weights : sub_solution_weights) {
      const double* antenna_weights = &weights[antenna * NChannelBlocks()];
      double weight_sum = 0.0;
      for (size_t channel = 0; channel != NChannelBlocks(); ++channel) {
        weight_sum += antenna_weights[channel];
        for (size_t result_index = 1; result_index < results_.size();
             ++result_index)
          results_[result_index].weights.emplace_back(antenna_weights[channel]);
      }
      rm_result.weights.emplace_back(weight_sum);
    }
  }
}

std::vector<Constraint::Result> FaradayConstraint::Apply(
    SolutionSpan& solutions, double,
    [[maybe_unused]] std::ostream* statStream) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NDirections());
  assert(solutions.shape(3) == 4);  // 2x2 full jones solutions
  assert(sub_solution_weights_.size() == NSubSolutions());
  const size_t n_fits = NAntennas() * NSubSolutions();
  aocommon::StaticFor<size_t> loop;
  loop.Run(0, n_fits, [&](size_t begin_index, size_t end_index, size_t thread) {
    std::vector<FitSample> scratch_space;
    for (size_t i = begin_index; i != end_index; ++i) {
      const size_t sub_solution = i % NSubSolutions();
      const size_t antenna = i / NSubSolutions();
      PerformFit(solutions, sub_solution, antenna, scratch_space);
    }
  });
  return results_;
}

void FaradayConstraint::PerformFit(SolutionSpan& solutions, size_t sub_solution,
                                   size_t antenna,
                                   std::vector<FitSample>& scratch_space) {
  const size_t n_channels = NChannelBlocks();
  scratch_space.clear();
  for (size_t i = 0; i != n_channels; ++i) {
    // The x values of the fitter must be ordered lowest to highest in
    // wavelength^2
    const size_t channel = n_channels - i - 1;
    const std::complex<double>* data =
        &solutions(channel, antenna, sub_solution, 0);
    double angle = RotationConstraint::FitRotation(data);
    float weight =
        sub_solution_weights_[sub_solution][antenna * n_channels + channel];
    const double frequency = frequencies_[channel];
    const double wavelength = aocommon::kSpeedOfLight / frequency;
    const double unit_faraday_depth = wavelength * wavelength;
    if (!std::isfinite(angle)) {
      angle = 0.0;
      weight = 0.0;
    }

    // The angle has an ambiguity in [-pi, pi] of pi (180 degrees), because a
    // rotation matrix for a rotation of pi negates the values, but
    // so does a matrix with values of -1 on the diagonal. Therefore,
    // the angle is multiplied by 2 so that it has a 2 pi ambiguity,
    // for which the fitter is written.
    scratch_space.emplace_back(
        FitSample{unit_faraday_depth, angle * 2.0, weight});
  }

  const double slope =
      0.5 * common::phase_fitting::FitSlope(scratch_space, fit_range_);
  results_.front().vals[antenna * NSubSolutions() + sub_solution] = slope;

  for (size_t i = 0; i != n_channels; ++i) {
    const size_t channel = n_channels - i - 1;
    std::complex<double>* data = &solutions(channel, antenna, sub_solution, 0);
    const double angle = scratch_space[i].x * slope;
    // Derotate the solutions for the found Faraday rotation
    // This calculates diagonal = diag(D R^1). For the X value this is therefore
    // X = D(0,0) R^-1(0,0) + D(1,0) R^-1(0,1)
    // X = D(0,0) cos(angle) + D(1,0) -sin(angle)
    // and thus
    // Y = D(0,1) sin(angle) + D(1,1) cos(angle)
    const double derotate = angle;
    const double sin_angle = std::sin(derotate);
    const double cos_angle = std::cos(derotate);
    std::array<std::complex<double>, 2> diagonal{
        data[0] * cos_angle - data[1] * sin_angle,
        data[2] * sin_angle + data[3] * cos_angle};
    ConstrainDiagonal(diagonal, diagonal_solution_type_);

    if (diagonal_solution_type_ != base::CalType::kRotation)
      StoreDiagonal(&results_[1], diagonal, channel, antenna, sub_solution,
                    NChannelBlocks(), NSubSolutions(), diagonal_solution_type_);
    RotationConstraint::SetRotation(data, angle);
    data[0] *= diagonal[0];
    data[1] *= diagonal[1];
    data[2] *= diagonal[0];
    data[3] *= diagonal[1];
  }
}

}  // namespace dp3::ddecal
