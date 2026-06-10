#include <vector>
#include <complex>

#include "ddecal/constraints/FaradayConstraint.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <aocommon/constants.h>
#include <aocommon/matrix2x2.h>

namespace dp3::ddecal {

using aocommon::MC2x2;

namespace {
constexpr size_t kNPolarizations = 4;
constexpr size_t kNAntennas = 7;
constexpr size_t kNSubSolutions = 3;
constexpr size_t kNChannels = 100;
constexpr double kStartFrequency = 30e6;
constexpr double kBandwidth = 45e6;

std::vector<ConstraintResult> Constrain(FaradayConstraint& constraint,
                                        SolutionTensor& solutions_tensor) {
  const std::vector<uint32_t> solutions_per_directdion(kNSubSolutions, 1);
  std::vector<double> frequencies(kNChannels);
  for (size_t i = 0; i != frequencies.size(); ++i) {
    frequencies[i] = kStartFrequency + i * kBandwidth / frequencies.size();
  }
  constraint.Initialize(kNAntennas, solutions_per_directdion, frequencies);
  constraint.SetWeights(std::vector<double>(kNAntennas * kNChannels, 1.0));
  dp3::ddecal::SolutionSpan solutions =
      aocommon::xt::CreateSpan(solutions_tensor);
  constraint.Apply(solutions, 0.0);
  return constraint.GetResult();
}

void Fill(SolutionTensor& solutions_tensor, const MC2x2 value) {
  constexpr size_t n_values = kNAntennas * kNChannels * kNSubSolutions;
  for (size_t i = 0; i != n_values; ++i) {
    value.AssignTo(&solutions_tensor.flat(i * 4));
  }
}

void ApplyFaradayRotation(SolutionTensor& solutions_tensor, size_t antenna,
                          size_t sub_solution, double rotation_value) {
  for (size_t channel = 0; channel != kNChannels; ++channel) {
    const double frequency =
        kStartFrequency + channel * kBandwidth / kNChannels;
    const double wavelength = aocommon::kSpeedOfLight / frequency;
    const double rotation = wavelength * wavelength * rotation_value;
    const double cos_rotation = std::cos(rotation);
    const double sin_rotation = std::sin(rotation);
    MC2x2 m(&solutions_tensor(channel, antenna, sub_solution, 0));
    m *= MC2x2(cos_rotation, -sin_rotation, sin_rotation, cos_rotation);
    m.AssignTo(&solutions_tensor(channel, antenna, sub_solution, 0));
  }
}

/**
 * Check if two values with pi-ambiguity are close. Note that this is
 * purposely not 2*pi, because Faraday values have a pi-ambiguity, not
 * 2*pi like phase values have.
 */
void CheckPiWrappedPhase(double result, double expected) {
  const double wrapped_result = std::fmod(result + 2.0 * M_PI, M_PI);
  double wrapped_expected = std::fmod(expected + 2.0 * M_PI, M_PI);
  // Bring expected to the same wrap as the wrapped result
  if (wrapped_expected - wrapped_result < -0.5 * M_PI)
    wrapped_expected += M_PI;
  else if (wrapped_expected - wrapped_result > 0.5 * M_PI)
    wrapped_expected -= M_PI;
  BOOST_CHECK_CLOSE_FRACTION(wrapped_result, wrapped_expected, 1.0e-3);
}

void CheckDiagonal(const std::vector<ConstraintResult>& result,
                   std::complex<double> diagonal_value_x,
                   std::complex<double> diagonal_value_y) {
  const ConstraintResult& amplitude = result[1];
  BOOST_CHECK_EQUAL(amplitude.vals.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);
  const ConstraintResult& phase = result[2];
  BOOST_CHECK_EQUAL(phase.vals.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);

  const double kExpectedPhaseX = std::arg(diagonal_value_x);
  const double kExpectedPhaseY = std::arg(diagonal_value_y);
  for (size_t antenna = 0; antenna != kNAntennas; ++antenna) {
    for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
         ++sub_solution) {
      const size_t index =
          (antenna * kNSubSolutions + sub_solution) * kNChannels * 2;
      const double* amplitude_ptr = &amplitude.vals[index];
      const double* phase_ptr = &phase.vals[index];
      for (size_t i = 0; i != kNChannels; ++i) {
        const double* amplitudes = &amplitude_ptr[i * 2];
        const double* phases = &phase_ptr[i * 2];
        BOOST_CHECK_CLOSE_FRACTION(amplitudes[0], std::abs(diagonal_value_x),
                                   1.0e-3);
        CheckPiWrappedPhase(phases[0], kExpectedPhaseX);
        BOOST_CHECK_CLOSE_FRACTION(amplitudes[1], std::abs(diagonal_value_y),
                                   1.0e-3);
        CheckPiWrappedPhase(phases[1], kExpectedPhaseY);
      }
    }
  }
}

void CheckFaraday(const std::vector<ConstraintResult>& result,
                  double antenna_rotation_step,
                  double sub_solution_rotation_step) {
  BOOST_REQUIRE_EQUAL(result.size(), 3u);
  const ConstraintResult& faraday = result[0];
  for (size_t antenna = 0; antenna != kNAntennas; ++antenna) {
    for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
         ++sub_solution) {
      const double rotation_value = antenna * antenna_rotation_step +
                                    sub_solution * sub_solution_rotation_step;
      const double result_value =
          faraday.vals[antenna * kNSubSolutions + sub_solution];
      BOOST_CHECK_CLOSE_FRACTION(rotation_value, result_value, 1.0e-3);
    }
  }
}

dp3::ddecal::SolutionTensor MakeSolutions(
    double antenna_rotation_step, double sub_solution_rotation_step,
    std::complex<double> diagonal_value_x,
    std::complex<double> diagonal_value_y) {
  dp3::ddecal::SolutionTensor solutions_tensor(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});
  Fill(solutions_tensor, MC2x2(diagonal_value_x, 0.0, 0.0, diagonal_value_y));
  for (size_t antenna = 0; antenna != kNAntennas; ++antenna) {
    for (size_t sub_solution = 0; sub_solution != kNSubSolutions;
         ++sub_solution) {
      const double rotation_value = antenna * antenna_rotation_step +
                                    sub_solution * sub_solution_rotation_step;
      ApplyFaradayRotation(solutions_tensor, antenna, sub_solution,
                           rotation_value);
    }
  }
  return solutions_tensor;
}

}  // namespace

using base::CalType;

BOOST_AUTO_TEST_SUITE(faraday_constraint)

BOOST_AUTO_TEST_CASE(unity_faraday_only) {
  dp3::ddecal::SolutionTensor solutions_tensor(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});
  Fill(solutions_tensor, MC2x2::Unity());

  FaradayConstraint constraint(base::CalType::kRotation, {});
  std::vector<ConstraintResult> result =
      Constrain(constraint, solutions_tensor);

  BOOST_REQUIRE_EQUAL(result.size(), 1u);
  const ConstraintResult& faraday = result[0];
  BOOST_CHECK_EQUAL(faraday.name, "rotationmeasure");
  BOOST_CHECK_EQUAL(faraday.axes, "ant,dir");
  BOOST_REQUIRE_EQUAL(faraday.dims.size(), 2u);
  BOOST_CHECK_EQUAL(faraday.dims[0], kNAntennas);
  BOOST_CHECK_EQUAL(faraday.dims[1], kNSubSolutions);
  BOOST_CHECK_EQUAL(faraday.vals.size(), kNAntennas * kNSubSolutions);
  for (double f : faraday.vals) BOOST_CHECK_LT(std::abs(f), 1e-6);
}

BOOST_AUTO_TEST_CASE(zero_rotation_with_diagonal) {
  dp3::ddecal::SolutionTensor solutions_tensor(
      {kNChannels, kNAntennas, kNSubSolutions, kNPolarizations});
  const std::complex<double> diagonal_value(3.0, 4.0);
  Fill(solutions_tensor, MC2x2(diagonal_value, 0.0, 0.0, diagonal_value));

  FaradayConstraint constraint(base::CalType::kDiagonal, {});
  std::vector<ConstraintResult> result =
      Constrain(constraint, solutions_tensor);

  BOOST_REQUIRE_EQUAL(result.size(), 3u);

  const ConstraintResult& faraday = result[0];
  BOOST_CHECK_EQUAL(faraday.name, "rotationmeasure");
  BOOST_CHECK_EQUAL(faraday.axes, "ant,dir");
  BOOST_REQUIRE_EQUAL(faraday.dims.size(), 2u);
  BOOST_CHECK_EQUAL(faraday.dims[0], kNAntennas);
  BOOST_CHECK_EQUAL(faraday.dims[1], kNSubSolutions);
  BOOST_CHECK_EQUAL(faraday.weights.size(), kNAntennas * kNSubSolutions);
  BOOST_CHECK_EQUAL(faraday.vals.size(), kNAntennas * kNSubSolutions);
  for (double f : faraday.vals) BOOST_CHECK_LT(std::abs(f), 1e-6);

  const ConstraintResult& amplitude = result[1];
  BOOST_CHECK_EQUAL(amplitude.name, "amplitude");
  BOOST_CHECK_EQUAL(amplitude.axes, "ant,dir,freq,pol");
  BOOST_REQUIRE_EQUAL(amplitude.dims.size(), 4u);
  BOOST_CHECK_EQUAL(amplitude.dims[0], kNAntennas);
  BOOST_CHECK_EQUAL(amplitude.dims[1], kNSubSolutions);
  BOOST_CHECK_EQUAL(amplitude.vals.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);
  BOOST_CHECK_EQUAL(amplitude.weights.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);
  for (double a : amplitude.vals)
    BOOST_CHECK_CLOSE_FRACTION(a, std::abs(diagonal_value), 1e-6);

  const ConstraintResult& phase = result[2];
  BOOST_CHECK_EQUAL(phase.name, "phase");
  BOOST_CHECK_EQUAL(phase.axes, "ant,dir,freq,pol");
  BOOST_REQUIRE_EQUAL(phase.dims.size(), 4u);
  BOOST_CHECK_EQUAL(phase.dims[0], kNAntennas);
  BOOST_CHECK_EQUAL(phase.dims[1], kNSubSolutions);
  BOOST_CHECK_EQUAL(phase.vals.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);
  BOOST_CHECK_EQUAL(phase.weights.size(),
                    kNAntennas * kNSubSolutions * kNChannels * 2);

  const double kExpectedPhase = std::arg(diagonal_value);
  for (double p : phase.vals)
    BOOST_CHECK_CLOSE_FRACTION(p, kExpectedPhase, 1e-6);
}

BOOST_AUTO_TEST_CASE(rotation) {
  constexpr std::complex<double> diagonal_value_x(3.0, 4.0);
  constexpr std::complex<double> diagonal_value_y(7.0, -2.0);
  constexpr double antenna_rotation_step = 0.1;
  constexpr double sub_solution_rotation_step = 0.07;
  dp3::ddecal::SolutionTensor solutions_tensor =
      MakeSolutions(antenna_rotation_step, sub_solution_rotation_step,
                    diagonal_value_x, diagonal_value_y);
  FaradayConstraint constraint(base::CalType::kDiagonal, {});
  std::vector<ConstraintResult> result =
      Constrain(constraint, solutions_tensor);

  CheckFaraday(result, antenna_rotation_step, sub_solution_rotation_step);
  CheckDiagonal(result, diagonal_value_x, diagonal_value_y);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace dp3::ddecal
