#ifndef TEC_OFFSET_DELAY_FITTING_H_
#define TEC_OFFSET_DELAY_FITTING_H_

#include <optional>
#include <span>
#include <string>

/**
 * @file
 * Contains functionality for fitting TEC, offset and delay to phase data.
 */

struct TecOffsetDelayValues {
  /** The 1/x coefficient that corresponds with the TEC value (unnormalized) */
  double a;
  /** The phase offset in radians */
  double b;
  /** The delay (from clock or otherwise) */
  double c;
};

/**
 * Perform a linear least-squares fit to solve the function y(x) = a / x + b + c
 * * x, where a is the TEC, b is the offset and c the delay.
 *
 * This function does not work well in the presence of wrapping: use the
 * @ref TecOffsetDelayGridSearch() for that.
 *
 * @returns an empty optional if the fit fails, otherwise
 */
std::optional<TecOffsetDelayValues> LinearTecOffsetDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights);

/**
 * Same as @ref LinearTecOffsetDelaySolve(), but assumes the delay offset b=0.
 */
std::optional<TecOffsetDelayValues> LinearTecDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights);

/**
 * Minimizes with a cost function of the sum of -cos(y - model), i.e. assuming a
 * von Mises distribution.
 */
std::optional<TecOffsetDelayValues> GradientTecOffsetDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights);

/**
 * Same as @ref GradientTecOffsetDelaySolve(), but assumes the delay offset b=0.
 */
std::optional<TecOffsetDelayValues> GradientTecDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights);

/**
 * Evaluates the function y(x) = a/x + b + c x and writes the result in @p y.
 * All values will be wrapped to the range [-pi, pi).
 */
void EvaluateLinearTecOffsetValues(const TecOffsetDelayValues& fit,
                                   std::span<const double> x,
                                   std::span<double> y);

enum class TecOffsetDelayFittingMethod { LeastSquares, VonMises };

/**
 * Perform a grid search for the TEC, offset and delay parameters. It optimizes
 * a, b, c in the function y(x) = a/x + b + c x, modulus 2pi, to minimize the
 * cost function.
 *
 * The cost function during the grid search is the sum of -cos(y_i - model_i),
 * i.e. it assumes the values follow the Van Mises distribution.
 * During refinement, the cost depends on the specified method:
 * least squares or Von Mises. When using least squares, the problem is linear
 * and therefore fast. The Van Mises uses a non-linear minimization function and
 * is therefore considerably slower.
 *
 * This function scales as max_wraps^2 x n (n being the data size). The
 * advantage of this method is that it is guaranteed to find the global minimum
 * in the noiseless case, but it is slow for large values of max_wraps. If this
 * becomes a bottleneck, a random sampling could be considered instead.
 *
 * @param x_data The frequency values, must be positive and in strictly
 * increasing order. At least 3 values must be specified.
 * @param y_data The phase values to be fitted (normally between 0 and 2pi, but
 * this is not required).
 * @param include_b true to fit b, false to assume b=0.
 * @param max_wraps The maximum number of wraps expected. The grid searched will
 * have dimensions of 2 x max_wraps, spread in such a way that the global
 * minimum is guaranteed to be found.
 */
TecOffsetDelayValues TecOffsetDelayGridSearch(
    std::span<const double> x_data, std::span<const double> y_data,
    std::span<const double> weights, bool include_b, size_t max_wraps,
    TecOffsetDelayFittingMethod method =
        TecOffsetDelayFittingMethod::LeastSquares);

double TecOffsetDelayCost(std::span<const double> x_data,
                          std::span<const double> y_data,
                          std::span<const double> weights,
                          const TecOffsetDelayValues& fit);

void PlotCostValues(const std::string& filename, std::span<const double> x_data,
                    std::span<const double> y_data,
                    std::span<const double> weights, size_t max_wraps);
#endif
