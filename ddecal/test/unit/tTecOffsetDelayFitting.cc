#include "ddecal/constraints/TecOffsetDelayFitting.h"

#include <cassert>  // can be removed once aocommon is updated
#include <complex>
#include <cmath>
#include <random>
#include <vector>

#include <aocommon/logger.h>
#include <aocommon/multiarray.h>

#include <boost/test/unit_test.hpp>

using aocommon::Logger;

namespace {

double Model(double x, double a, double b, double c) {
  double f = a / x + b + c * x;
  double wrap = std::fmod(f, 2.0 * M_PI);
  if (wrap >= M_PI)
    wrap -= 2.0 * M_PI;
  else if (wrap < -M_PI)
    wrap += 2.0 * M_PI;
  return wrap;
}

double WrapOffset(double x) {
  if (x >= M_PI)
    x -= 2.0 * M_PI;
  else if (x < -M_PI)
    x += 2.0 * M_PI;
  return x;
}

const std::array<TecOffsetDelayFittingMethod, 2> kMethods{
    TecOffsetDelayFittingMethod::LeastSquares,
    TecOffsetDelayFittingMethod::VonMises};

}  // namespace

BOOST_AUTO_TEST_SUITE(tec_offset_delay_fitting)

BOOST_AUTO_TEST_CASE(refinement_fit) {
  std::vector<double> x;
  for (size_t i = 0; i != 50; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);
  std::vector<double> y(x.size());
  const TecOffsetDelayValues true_values{.a = 0.7, .b = 0.31, .c = -0.03};
  EvaluateLinearTecOffsetValues(true_values, x, y);

  std::optional<TecOffsetDelayValues> result =
      LinearTecOffsetDelaySolve(x, y, w);
  BOOST_REQUIRE(result.has_value());
  BOOST_CHECK_CLOSE_FRACTION(result->a, true_values.a, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(WrapOffset(result->b), true_values.b, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result->c, true_values.c, 1e-6);

  result = GradientTecOffsetDelaySolve(x, y, w);
  BOOST_REQUIRE(result.has_value());
  BOOST_CHECK_CLOSE_FRACTION(result->a, true_values.a, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(WrapOffset(result->b), true_values.b, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result->c, true_values.c, 1e-6);
}

BOOST_AUTO_TEST_CASE(refinement_fit_with_noise) {
  std::vector<double> x;
  for (size_t i = 0; i != 50; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);

  std::mt19937 rnd;
  std::normal_distribution<double> gaus(0.0, 0.1);

  for (size_t repeat = 0; repeat != 5; ++repeat) {
    std::vector<double> y;
    const double true_a = 0.2, true_b = 1.0, true_c = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
      y.push_back(Model(x[i], true_a, true_b, true_c) + gaus(rnd));
    }

    for (TecOffsetDelayFittingMethod method : kMethods) {
      std::optional<TecOffsetDelayValues> result =
          method == TecOffsetDelayFittingMethod::LeastSquares
              ? LinearTecOffsetDelaySolve(x, y, w)
              : GradientTecOffsetDelaySolve(x, y, w);
      BOOST_REQUIRE(result.has_value());
      BOOST_CHECK_LT(result->a, true_a + 0.4);
      BOOST_CHECK_GT(result->a, true_a - 0.4);
      BOOST_CHECK_LT(WrapOffset(result->b), true_b + 0.1);
      BOOST_CHECK_GT(WrapOffset(result->b), true_b - 0.1);

      BOOST_CHECK_LT(result->c, true_c + 0.01);
      BOOST_CHECK_GT(result->c, true_c - 0.01);
    }
  }
}

BOOST_AUTO_TEST_CASE(demonstrate_von_mises, *boost::unit_test::disabled()) {
  std::vector<double> x;
  for (size_t i = 0; i != 100; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);

  std::mt19937 rnd;
  std::uniform_real_distribution<double> uniform_dist(0.0, 0.5);

  double a_err[2] = {0.0, 0.0};
  double b_err[2] = {0.0, 0.0};
  double c_err[2] = {0.0, 0.0};
  double cost[2] = {0.0, 0.0};
  const size_t kNRepeats = 100000;
  size_t wrap_count = 0;
  for (size_t repeat = 0; repeat != kNRepeats; ++repeat) {
    const double uniform = /*uniform_dist(rnd) * 1.0 + */ 2.0;
    std::normal_distribution<double> gaus(0.0, uniform * uniform);
    std::vector<double> y;
    const double true_a = uniform_dist(rnd) * 0.5;
    const double true_b = uniform_dist(rnd) * 2.0 - 1.0;
    const double true_c = -uniform_dist(rnd) * 0.01;
    for (size_t i = 0; i < x.size(); ++i) {
      const double f = true_a / x[i] + true_b + true_c * x[i];
      if (f < -M_PI || f >= M_PI) ++wrap_count;
      // Add Gaussian noise in real space, then convert back to phase space
      const std::complex<double> z(std::polar(1.0, f) +
                                   std::complex(gaus(rnd), gaus(rnd)));
      const double y_sample = std::arg(z);
      y.push_back(WrapOffset(y_sample));
    }

    for (size_t m = 0; m != 2; ++m) {
      std::optional<TecOffsetDelayValues> result =
          m == 0 ? LinearTecOffsetDelaySolve(x, y, w)
                 : GradientTecOffsetDelaySolve(x, y, w);
      BOOST_REQUIRE(result.has_value());
      cost[m] += TecOffsetDelayCost(x, y, w, *result);
      a_err[m] += (result->a - true_a) * (result->a - true_a);
      b_err[m] +=
          (WrapOffset(result->b) - true_b) * (WrapOffset(result->b) - true_b);
      c_err[m] += (result->c - true_c) * (result->c - true_c);
    }
  }

  Logger::Info << "Average nr of samples that wrapped: "
               << double(wrap_count) / kNRepeats << " / " << x.size() << '\n';
  for (size_t m = 0; m != 2; ++m) {
    cost[m] = cost[m] / (kNRepeats * x.size());
    a_err[m] = std::sqrt(a_err[m] / kNRepeats);
    b_err[m] = std::sqrt(b_err[m] / kNRepeats);
    c_err[m] = std::sqrt(c_err[m] / kNRepeats);
  }
  Logger::Info << "              Least sq\tVon Mises\n"
               << "Cost/sample:  " << cost[0] << '\t' << cost[1] << '\n'
               << "TEC error:    " << a_err[0] << '\t' << a_err[1] << '\n'
               << "Offset error: " << b_err[0] << '\t' << b_err[1] << '\n'
               << "Delay error:  " << c_err[0] << '\t' << c_err[1] << '\n';
}

BOOST_AUTO_TEST_CASE(fit_exact_with_b) {
  std::vector<double> x;
  for (size_t i = 0; i != 50; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);
  std::vector<double> y(x.size());
  const TecOffsetDelayValues true_values{.a = 12, .b = 3.1, .c = -0.3};
  EvaluateLinearTecOffsetValues(true_values, x, y);

  for (TecOffsetDelayFittingMethod method : kMethods) {
    const TecOffsetDelayValues result =
        TecOffsetDelayGridSearch(x, y, w, true, 10, method);
    BOOST_CHECK_CLOSE_FRACTION(result.a, true_values.a, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(WrapOffset(result.b), true_values.b, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(result.c, true_values.c, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(fit_exact_without_b) {
  std::vector<double> x;
  for (size_t i = 0; i != 50; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);
  std::vector<double> y(x.size());
  const TecOffsetDelayValues true_values{.a = 11, .b = 0.0, .c = 0.25};
  EvaluateLinearTecOffsetValues(true_values, x, y);

  for (TecOffsetDelayFittingMethod method : kMethods) {
    const TecOffsetDelayValues result =
        TecOffsetDelayGridSearch(x, y, w, false, 10, method);
    BOOST_CHECK_CLOSE_FRACTION(result.a, true_values.a, 1e-6);
    BOOST_CHECK_EQUAL(result.b, 0.0);
    BOOST_CHECK_CLOSE_FRACTION(result.c, true_values.c, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(fit_with_noise) {
  std::vector<double> x;
  for (size_t i = 0; i != 100; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);

  std::mt19937 rnd;
  std::normal_distribution<double> gaus(0.0, 0.1);

  for (size_t repeat = 0; repeat != 5; ++repeat) {
    std::vector<double> y;
    const double true_a = -10, true_b = -2.0, true_c = 0.2;
    for (size_t i = 0; i < x.size(); ++i) {
      y.push_back(Model(x[i], true_a, true_b, true_c) + gaus(rnd));
    }

    for (TecOffsetDelayFittingMethod method : kMethods) {
      const TecOffsetDelayValues result =
          TecOffsetDelayGridSearch(x, y, w, true, 4, method);
      BOOST_CHECK_LT(result.a, true_a + 0.4);
      BOOST_CHECK_GT(result.a, true_a - 0.4);
      BOOST_CHECK_LT(WrapOffset(result.b), true_b + 0.05);
      BOOST_CHECK_GT(WrapOffset(result.b), true_b - 0.05);

      BOOST_CHECK_LT(result.c, true_c + 0.01);
      BOOST_CHECK_GT(result.c, true_c - 0.01);
    }
  }
}

/**
 * Test set to make sure the interface supports using a MultiArray.
 */
BOOST_AUTO_TEST_CASE(demonstrate_multiarray) {
  aocommon::MultiArray<double, double, double> data(50, false);
  for (size_t i = 0; i != data.Size(); ++i) data.Get<0>()[i] = i + 1;
  std::fill_n(data.Get<2>(), data.Size(), 1.0);
  const TecOffsetDelayValues true_values{.a = 0.3, .b = -1.8, .c = -0.08};
  EvaluateLinearTecOffsetValues(true_values, data.Span<0>(), data.Span<1>());

  const TecOffsetDelayValues result = TecOffsetDelayGridSearch(
      data.Span<0>(), data.Span<1>(), data.Span<2>(), true, 10);
  BOOST_CHECK_CLOSE_FRACTION(result.a, true_values.a, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(WrapOffset(result.b), true_values.b, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(result.c, true_values.c, 1e-6);
}

/**
 * This test can be used to produce a file ('cost.txt') that can be used to plot
 * the cost as function of a and c.
 */
BOOST_AUTO_TEST_CASE(plot, *boost::unit_test::disabled()) {
  std::vector<double> x;
  for (size_t i = 0; i != 100; ++i) x.emplace_back(i + 1);
  std::vector<double> w(x.size(), 1.0);

  std::mt19937 rnd;
  std::normal_distribution<double> gaus(0.0, 0.5);

  std::vector<double> y;
  const double true_a = -10, true_b = 0.0, true_c = 0.2;
  for (size_t i = 0; i < x.size(); ++i) {
    y.push_back(Model(x[i], true_a, true_b, true_c) + gaus(rnd));
  }
  PlotCostValues("cost.txt", x, y, w, 10);
}

BOOST_AUTO_TEST_SUITE_END()
