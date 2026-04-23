#include "TecOffsetDelayFitting.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

namespace {
struct LinearFitWorkspace {
  ~LinearFitWorkspace() noexcept {
    // Even if errors occur during/before allocation it is ok to call the
    // gsl_*_free functions, as they do not fail when a nullptr is given.
    gsl_multifit_linear_free(work);
    gsl_matrix_free(x);
    gsl_vector_free(y);
    gsl_vector_free(coefficients);
    gsl_matrix_free(covariance);
    gsl_vector_free(weights);
  }

  gsl_matrix* x = nullptr;
  gsl_vector* y = nullptr;
  gsl_vector* coefficients = nullptr;
  gsl_matrix* covariance = nullptr;
  gsl_vector* weights = nullptr;
  gsl_multifit_linear_workspace* work = nullptr;
};

struct NonLinearFitWorkspace {
  ~NonLinearFitWorkspace() noexcept {
    gsl_vector_free(initial_values);
    gsl_multimin_fdfminimizer_free(minimizer);
  }

  LinearFitWorkspace linear_ws;
  std::span<const double> x;
  std::vector<double> y;
  std::span<const double> weights;
  size_t dim = 3;
  gsl_vector* initial_values = nullptr;
  gsl_multimin_fdfminimizer* minimizer = nullptr;
};

void AllocateLinearSolver(LinearFitWorkspace& work, std::span<const double> x,
                          std::span<const double> weights, bool include_b) {
  const size_t n = x.size();
  const size_t p = include_b ? 3 : 2;
  // Allocate design matrix x (n rows, p columns)
  work.x = gsl_matrix_alloc(n, p);
  work.y = gsl_vector_alloc(n);
  work.coefficients = gsl_vector_alloc(p);
  work.covariance = gsl_matrix_alloc(p, p);
  work.weights = gsl_vector_alloc(n);  // uniform weights = 1
  work.work = gsl_multifit_linear_alloc(n, p);

  if (work.x == nullptr || work.y == nullptr || work.coefficients == nullptr ||
      work.covariance == nullptr || work.weights == nullptr ||
      work.work == nullptr)
    throw std::bad_alloc();

  if (include_b) {
    for (size_t i = 0; i < n; ++i) {
      assert(x[i] != 0.0);
      const double inv_x = 1.0 / x[i];

      gsl_matrix_set(work.x, i, 0, inv_x);  // coefficient for a
      gsl_matrix_set(work.x, i, 1, 1.0);    // coefficient for b
      gsl_matrix_set(work.x, i, 2, x[i]);   // coefficient for c
      gsl_vector_set(work.weights, i, weights[i]);
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      assert(x[i] != 0.0);
      const double inv_x = 1.0 / x[i];

      gsl_matrix_set(work.x, i, 0, inv_x);  // coefficient for a
      gsl_matrix_set(work.x, i, 1, x[i]);   // coefficient for c
      gsl_vector_set(work.weights, i, weights[i]);
    }
  }
}

void SetLinearSolverYValues(LinearFitWorkspace& work,
                            std::span<const double> y) {
  const size_t n = y.size();
  for (size_t i = 0; i < n; ++i) {
    gsl_vector_set(work.y, i, y[i]);
  }
}

void AllocateGradientSolver(NonLinearFitWorkspace& ws,
                            std::span<const double> x,
                            std::span<const double> y,
                            std::span<const double> weights, bool include_b) {
  AllocateLinearSolver(ws.linear_ws, x, weights, include_b);
  ws.x = x;
  ws.y = std::vector<double>(y.begin(), y.end());
  ws.weights = weights;
  ws.dim = include_b ? 3 : 2;

  ws.initial_values = gsl_vector_alloc(ws.dim);

  const gsl_multimin_fdfminimizer_type* T =
      gsl_multimin_fdfminimizer_vector_bfgs2;
  ws.minimizer = gsl_multimin_fdfminimizer_alloc(T, ws.dim);

  if (!ws.minimizer || !ws.initial_values) throw std::bad_alloc();
}

/**
 * Perform a linear solve from an already initialized and assigned workspace.
 *
 * The linear solving is separated in allocation, setting and the solving
 * functions so that the grid search (which does repeated linear solves of the
 * same size) can reuse the allocations.
 * @returns an unset optional on failure, otherwise a TecOffsetDelayFitResult
 * with the found coefficients.
 */
std::optional<TecOffsetDelayValues> LinearSolve(LinearFitWorkspace& ws,
                                                bool include_b) {
  // Solve the linear system
  double chisq;
  const int status = gsl_multifit_wlinear(
      ws.x, ws.weights, ws.y, ws.coefficients, ws.covariance, &chisq, ws.work);

  const bool success = (status == GSL_SUCCESS);
  if (success) {
    TecOffsetDelayValues result;
    result.a = gsl_vector_get(ws.coefficients, 0);
    if (include_b) {
      result.b = gsl_vector_get(ws.coefficients, 1);
      result.c = gsl_vector_get(ws.coefficients, 2);
    } else {
      result.b = 0.0;
      result.c = gsl_vector_get(ws.coefficients, 1);
    }
    return result;
  } else {
    return {};
  }
}

inline double Wrap(double x) {
  constexpr double twopi = 2.0 * M_PI;
  double wrap = std::fmod(x, twopi);
  if (wrap >= M_PI)
    wrap -= twopi;
  else if (wrap < -M_PI)
    wrap += twopi;
  return wrap;
}

inline constexpr double UnwrappedModel(double x, double a, double b, double c) {
  return a / x + b + c * x;
}

double VonMisesCost(const gsl_vector* v, void* params) {
  const NonLinearFitWorkspace& data =
      *static_cast<NonLinearFitWorkspace*>(params);
  const TecOffsetDelayValues values{
      .a = gsl_vector_get(v, 0),
      .b = data.dim == 3 ? gsl_vector_get(v, 1) : 0.0,
      .c = gsl_vector_get(v, data.dim - 1)};
  return TecOffsetDelayCost(data.x, data.y, data.weights, values);
}

void VonMisesCostDerivative(const gsl_vector* v, void* params, gsl_vector* df) {
  const NonLinearFitWorkspace& data =
      *static_cast<NonLinearFitWorkspace*>(params);

  const double a = gsl_vector_get(v, 0);
  const double b = data.dim == 3 ? gsl_vector_get(v, 1) : 0.0;
  const double c = gsl_vector_get(v, data.dim - 1);

  double da = 0.0, db = 0.0, dc = 0.0;
  for (size_t i = 0; i < data.x.size(); ++i) {
    const double xi = data.x[i];
    const double r = data.y[i] - (a / xi + b + c * xi);
    const double s = -std::sin(r) * data.weights[i];

    da += s / xi;
    db += s;
    dc += s * xi;
  }

  gsl_vector_set(df, 0, da);
  if (data.dim == 3) gsl_vector_set(df, 1, db);
  gsl_vector_set(df, data.dim - 1, dc);
}

void VonMisesCostFdf(const gsl_vector* v, void* params, double* f,
                     gsl_vector* df) {
  *f = VonMisesCost(v, params);
  VonMisesCostDerivative(v, params, df);
}

}  // namespace

double TecOffsetDelayCost(std::span<const double> x_data,
                          std::span<const double> y_data,
                          std::span<const double> weights,
                          const TecOffsetDelayValues& fit) {
  double cost = 0.0;
  for (size_t i = 0; i != x_data.size(); ++i) {
    const double x = x_data[i];
    const double f = fit.a / x + fit.b + fit.c * x;
    const double r = y_data[i] - f;
    cost += 1.0 - std::cos(r) * weights[i];
  }
  return cost;
}

std::optional<TecOffsetDelayValues> LinearTecOffsetDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights) {
  assert(x.size() == y.size() && x.size() == weights.size());
  assert(x.size() >= 3);
  assert(x.front() != x.back());

  LinearFitWorkspace work;
  AllocateLinearSolver(work, x, weights, true);
  SetLinearSolverYValues(work, y);
  return LinearSolve(work, true);
}

std::optional<TecOffsetDelayValues> LinearTecDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights) {
  assert(x.size() == y.size() && x.size() == weights.size());
  assert(x.size() >= 3);
  assert(x.front() != x.back());

  LinearFitWorkspace work;
  AllocateLinearSolver(work, x, weights, false);
  SetLinearSolverYValues(work, y);
  return LinearSolve(work, false);
}

std::optional<TecOffsetDelayValues> GradientTecOffsetDelaySolve(
    NonLinearFitWorkspace& workspace) {
  gsl_multimin_function_fdf func;
  func.n = workspace.dim;
  func.f = VonMisesCost;
  func.df = VonMisesCostDerivative;
  func.fdf = VonMisesCostFdf;
  func.params = &workspace;

  const bool include_b = workspace.dim == 3;

  // Initialize from linear solution
  SetLinearSolverYValues(workspace.linear_ws, workspace.y);
  const std::optional<TecOffsetDelayValues> initial_values =
      LinearSolve(workspace.linear_ws, include_b);
  if (!initial_values) return {};
  gsl_vector_set(workspace.initial_values, 0, initial_values->a);
  if (include_b) {
    gsl_vector_set(workspace.initial_values, 1, initial_values->b);
    gsl_vector_set(workspace.initial_values, 2, initial_values->c);
  } else {
    gsl_vector_set(workspace.initial_values, 1, initial_values->c);
  }

  const double step_size = 1e-3;
  const double tolerance = 1e-5;
  gsl_multimin_fdfminimizer_set(workspace.minimizer, &func,
                                workspace.initial_values, step_size, tolerance);

  size_t iter = 0;
  int status;

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(workspace.minimizer);
    if (status) break;

    status = gsl_multimin_test_gradient(workspace.minimizer->gradient, 1e-5);
  } while (status == GSL_CONTINUE && iter < 200);

  // Ignore status: if the linear fit was perfect, GSL will report an error, but
  // the values will be accurate.
  TecOffsetDelayValues fit;
  fit.a = gsl_vector_get(workspace.minimizer->x, 0);
  if (include_b) {
    fit.b = gsl_vector_get(workspace.minimizer->x, 1);
    fit.c = gsl_vector_get(workspace.minimizer->x, 2);
  } else {
    fit.b = 0.0;
    fit.c = gsl_vector_get(workspace.minimizer->x, 1);
  }
  return fit;
}

std::optional<TecOffsetDelayValues> GradientTecOffsetDelaySolve(
    std::span<const double> x, std::span<const double> y,
    std::span<const double> weights) {
  NonLinearFitWorkspace ws;
  AllocateGradientSolver(ws, x, y, weights, true);
  return GradientTecOffsetDelaySolve(ws);
}

TecOffsetDelayValues TecOffsetDelayGridSearch(
    std::span<const double> x_data, std::span<const double> y_data,
    std::span<const double> weights, bool include_b, size_t max_wraps,
    TecOffsetDelayFittingMethod method) {
  assert(x_data.size() == y_data.size() && x_data.size() == weights.size());
  assert(x_data.size() >= 3);
  assert(x_data.front() != x_data.back());

  NonLinearFitWorkspace gradient_workspace;
  LinearFitWorkspace linear_workspace;
  if (method == TecOffsetDelayFittingMethod::LeastSquares) {
    AllocateLinearSolver(linear_workspace, x_data, weights, include_b);
  } else {
    AllocateGradientSolver(gradient_workspace, x_data, y_data, weights,
                           include_b);
  }

  const double x_range = x_data.back() - x_data.front();
  const double max_y_range = 2 * M_PI * max_wraps;
  const double max_a =
      max_y_range / (1.0 / x_data.front() - 1.0 / x_data.back());
  const double max_c = max_y_range / x_range;

  std::vector<double> residual_data(y_data.size());
  double best_cost = std::numeric_limits<double>::max();
  TecOffsetDelayValues best_fit{.a = 0.0, .b = 0.0, .c = 0.0};

  // Perform the grid search
  for (size_t a_index = 0; a_index != max_wraps * 4 + 2; ++a_index) {
    // Calculate position in range of [-1, 1]
    const double a_position = double(a_index) / double(max_wraps * 2) - 1.0;
    const double a = a_position * max_a;

    for (size_t c_index = 0; c_index != max_wraps * 4 + 2; ++c_index) {
      const double c_position = double(c_index) / double(max_wraps * 2) - 1.0;
      const double c = c_position * max_c;

      // The a and c values are used as starting guess by subtracting the
      // corresponding model from the data and wrapping them all to [-pi, pi].
      // The linear fit is then done on the residual to determine coefficient
      // offsets. The b values are solved for so we don't need to unwrap for b.
      for (size_t i = 0; i != x_data.size(); ++i)
        residual_data[i] =
            Wrap(y_data[i] - UnwrappedModel(x_data[i], a, 0.0, c));

      std::optional<TecOffsetDelayValues> residual_fit;
      if (method == TecOffsetDelayFittingMethod::LeastSquares) {
        SetLinearSolverYValues(linear_workspace, residual_data);
        residual_fit = LinearSolve(linear_workspace, include_b);
      } else {
        std::swap(residual_data, gradient_workspace.y);
        residual_fit = GradientTecOffsetDelaySolve(gradient_workspace);
      }
      if (residual_fit) {
        const TecOffsetDelayValues fit = {.a = residual_fit->a + a,
                                          .b = residual_fit->b,
                                          .c = residual_fit->c + c};
        const double cost = TecOffsetDelayCost(x_data, y_data, weights, fit);
        if (cost < best_cost) {
          best_cost = cost;
          best_fit = fit;
        }
      }
    }
  }

  // Refinement: now that we're close to the optimized values, some values might
  // still wrap slightly differently then when using the grid values. Perform a
  // few more iteration on the best result.
  for (size_t repeat = 0; repeat != 5; ++repeat) {
    for (size_t i = 0; i != x_data.size(); ++i)
      residual_data[i] =
          Wrap(y_data[i] -
               UnwrappedModel(x_data[i], best_fit.a, best_fit.b, best_fit.c));
    std::optional<TecOffsetDelayValues> residual_fit;
    if (method == TecOffsetDelayFittingMethod::LeastSquares) {
      SetLinearSolverYValues(linear_workspace, residual_data);
      residual_fit = LinearSolve(linear_workspace, include_b);
    } else {
      std::swap(residual_data, gradient_workspace.y);
      residual_fit = GradientTecOffsetDelaySolve(gradient_workspace);
    }
    if (residual_fit) {
      const TecOffsetDelayValues fit = {residual_fit->a + best_fit.a,
                                        residual_fit->b + best_fit.b,
                                        residual_fit->c + best_fit.c};
      const double cost = TecOffsetDelayCost(x_data, y_data, weights, fit);
      if (cost < best_cost) {
        best_cost = cost;
        best_fit = fit;
      }
    }
  }

  return best_fit;
}

void PlotCostValues(const std::string& filename, std::span<const double> x_data,
                    std::span<const double> y_data,
                    std::span<const double> weights, size_t max_wraps) {
  assert(x_data.size() == y_data.size() && x_data.size() == weights.size());
  assert(x_data.size() >= 3);
  assert(x_data.front() != x_data.back());

  std::ofstream str(filename);
  const double x_range = x_data.back() - x_data.front();
  const double max_y_range = 2 * M_PI * max_wraps;
  const double max_a =
      max_y_range / (1.0 / x_data.front() - 1.0 / x_data.back());
  const double max_c = max_y_range / x_range;
  const size_t oversample = 20;
  // Perform the grid search
  TecOffsetDelayValues values{.a = 0.0, .b = 0.0, .c = 0.0};
  for (size_t a_index = 0; a_index != max_wraps * oversample * 2 + 2;
       ++a_index) {
    // Calculate position in range of [-1, 1]
    const double a_position =
        double(a_index) / double(max_wraps * oversample) - 1.0;
    values.a = a_position * max_a;

    for (size_t c_index = 0; c_index != max_wraps * oversample * 2 + 2;
         ++c_index) {
      const double c_position =
          double(c_index) / double(max_wraps * oversample) - 1.0;
      values.c = c_position * max_c;
      const double cost = TecOffsetDelayCost(x_data, y_data, weights, values);
      str << values.a << '\t' << values.c << '\t' << cost << '\n';
    }
  }
}

void EvaluateLinearTecOffsetValues(const TecOffsetDelayValues& fit,
                                   std::span<const double> x,
                                   std::span<double> y) {
  assert(x.size() == y.size());

  for (size_t i = 0; i != x.size(); ++i) {
    y[i] = Wrap(UnwrappedModel(x[i], fit.a, fit.b, fit.c));
  }
}
