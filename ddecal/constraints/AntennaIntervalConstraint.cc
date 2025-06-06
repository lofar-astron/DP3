#include "AntennaIntervalConstraint.h"

#include <algorithm>
#include <stdexcept>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmasked_view.hpp>
#include <xtensor/xview.hpp>

namespace dp3::ddecal {

void AntennaIntervalConstraint::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  assert(n_antennas != 0);
  assert(!solutions_per_direction.empty());
  if (antenna_averaging_factors_.size() != n_antennas)
    throw std::runtime_error(
        "The 'antenna_averaging_factors' setting for the antenna interval "
        "constraint had an incorrect number of values: found " +
        std::to_string(antenna_averaging_factors_.size()) +
        " values, but the data has " + std::to_string(n_antennas) +
        " antennas");

  const uint32_t max_solutions_in_direction = *std::max_element(
      solutions_per_direction.begin(), solutions_per_direction.end());
  const std::pair min_max_antenna_interval = std::minmax_element(
      antenna_averaging_factors_.begin(), antenna_averaging_factors_.end());
  if (*min_max_antenna_interval.first == 0)
    throw std::runtime_error(
        "One of the values in the 'antenna_averaging_factors' setting is zero "
        "in "
        "the antenna interval constraint");
  // An antenna averaging factor is not allowed to exceed the max nr of
  // direction subsolutions.
  if (*min_max_antenna_interval.second > max_solutions_in_direction)
    throw std::runtime_error(
        "The maximum number of solutions in a direction is " +
        std::to_string(max_solutions_in_direction) +
        ", whereas the maximum averaging factor in the antenna interval "
        "constraint is set to a higher value of " +
        std::to_string(*min_max_antenna_interval.second));

  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);
}

std::vector<Constraint::Result> AntennaIntervalConstraint::Apply(
    SolutionSpan& solutions, [[maybe_unused]] double time,
    [[maybe_unused]] std::ostream* stat_stream) {
  const size_t n_channels = solutions.shape(0);
  const size_t n_antennas = solutions.shape(1);
  const size_t n_polarizations = solutions.shape(3);
  const size_t n_directions = NDirections();
  // All n_channels/polarizations are going to be added at once into
  // vectors. Declaring these variables outside the loop allows reusing
  // memory.
  xt::xtensor<dcomplex, 2> solution_set;
  xt::xtensor<size_t, 2> solution_set_is_finite;
  std::array<size_t, 2> shape = {n_channels, n_polarizations};
  xt::xtensor<dcomplex, 2> averages(shape);
  xt::xtensor<size_t, 2> counts(shape);
  for (size_t antenna = 0; antenna != n_antennas; ++antenna) {
    const size_t n_averaged_intervals = antenna_averaging_factors_[antenna];
    size_t sub_solution_index = 0;
    for (size_t direction = 0; direction != n_directions; ++direction) {
      size_t interval_start = 0;
      while (interval_start < GetSubSolutions(direction)) {
        averages.fill(0.0);
        counts.fill(0);

        // Calculate the sum of solutions over the intervals
        const size_t end = std::min<size_t>(
            n_averaged_intervals + interval_start, GetSubSolutions(direction));
        const size_t interval_size = end - interval_start;
        for (size_t interval = 0; interval != interval_size; ++interval) {
          solution_set = xt::view(solutions, xt::all(), antenna,
                                  sub_solution_index + interval, xt::all());
          solution_set_is_finite = xt::isfinite(solution_set);
          xt::masked_view(solution_set, !solution_set_is_finite)
              .fill(std::complex<double>(0.0, 0.0));

          averages += solution_set;
          counts += solution_set_is_finite;
        }
        averages /= xt::cast<double>(counts);

        // Assign average to involved intervals
        for (size_t interval = 0; interval != interval_size; ++interval) {
          xt::view(solutions, xt::all(), antenna, sub_solution_index + interval,
                   xt::all()) = averages;
        }
        sub_solution_index += interval_size;
        interval_start += interval_size;
      }
    }
  }
  return std::vector<Constraint::Result>();
}

}  // namespace dp3::ddecal
