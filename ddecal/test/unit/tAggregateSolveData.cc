// Copyright (C) 2026 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Equivalence test for ddecal.aggregate: aggregating the model data
// and visibilities over a solution interval should give the
// same results

#include "ddecal/CalibrationAggregation.h"
#include "ddecal/gain_solvers/AggregateSolveData.h"
#include "ddecal/gain_solvers/IterativeDiagonalSolver.h"
#include "ddecal/gain_solvers/IterativeScalarSolver.h"
#include "ddecal/gain_solvers/ScalarSolver.h"
#include "ddecal/gain_solvers/SolveData.h"
#include "ddecal/gain_solvers/SolverTools.h"

#include "base/DPBuffer.h"
#include "schaapcommon/threading/threadpool.h"

#include <boost/test/unit_test.hpp>

#include <aocommon/matrix2x2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <random>
#include <vector>

using dp3::base::DPBuffer;
using dp3::ddecal::DuoSolveData;
using dp3::ddecal::FullSolveData;
using dp3::ddecal::SolverDataUse;
using dp3::ddecal::UniSolveData;

namespace {

constexpr size_t kNPolarizations = 4;
constexpr size_t kNChannels = 4;
constexpr size_t kNChannelBlocks = 2;
constexpr size_t kNBaselines = 3;
constexpr size_t kNAntennas = 3;
// The interval length only has to be large enough for the random model data
// to span all lin = ndir * npol = 8 dimensions of the reduced Gram matrix.
constexpr size_t kNTimesteps = 8;  // the solution interval length
constexpr size_t kNDirections = 2;

// Baseline 2 is an autocorrelation (2 -> 2), which the solvers skip.
const std::vector<int> kAntennas1 = {0, 0, 0};
const std::vector<int> kAntennas2 = {1, 2, 0};

const std::array<size_t, 3> kShape{kNBaselines, kNChannels, kNPolarizations};

std::vector<std::string> DirectionNames() { return {"dir0", "dir1"}; }

/**
 * Create weighted data and model data for one solution interval, in the layout
 * ddecal::AssignAndWeight produces: one DPBuffer per timestep, each
 * [n_baselines, n_channels, n_correlations], with one model buffer per
 * direction.
 *
 * The visibilities are generated from known scalar gains applied to the model,
 * so the problem is well posed and the solver converges. Purely random data
 * makes the scalar solver iterate without reaching a stopping criterion.
 */
std::vector<DPBuffer> MakeWeightedBuffers(
    const std::vector<std::string>& names,
    const std::vector<size_t>& n_solutions, size_t n_timesteps = kNTimesteps) {
  std::mt19937 mt(42);
  std::uniform_real_distribution<float> uniform(-0.5, 0.5);
  const auto rnd = [&]() {
    return std::complex<float>(uniform(mt), uniform(mt));
  };

  // Known per-antenna scalar gains, one per solution, so the data can be
  // generated consistently.
  const size_t n_sol = n_solutions.empty() ? 1 : n_solutions.front();
  std::vector<std::vector<std::complex<double>>> gains(
      kNAntennas, std::vector<std::complex<double>>(n_sol));
  for (size_t ant = 0; ant < kNAntennas; ++ant) {
    for (size_t s = 0; s < n_sol; ++s) {
      // Keep the gains away from zero and of similar magnitude so the
      // normal equations stay well conditioned.
      gains[ant][s] =
          1.0 + 0.1 * std::complex<double>(static_cast<double>(ant + 1),
                                           static_cast<double>(s + 1));
    }
  }

  std::vector<std::unique_ptr<DPBuffer>> unweighted;
  std::vector<DPBuffer> weighted;
  for (size_t t = 0; t < n_timesteps; ++t) {
    unweighted.push_back(std::make_unique<DPBuffer>());
    weighted.emplace_back();
    unweighted.back()->GetData().resize(kShape);
    unweighted.back()->GetWeights().resize(kShape);
    unweighted.back()->GetFlags().resize(kShape);
    unweighted.back()->GetWeights().fill(1.0f);
    unweighted.back()->GetFlags().fill(false);
    for (size_t d = 0; d < names.size(); ++d) {
      unweighted.back()->AddData(names[d]);
      unweighted.back()->GetData(names[d]).resize(kShape);
    }

    for (size_t baseline = 0; baseline < kNBaselines; ++baseline) {
      const size_t ant1 = static_cast<size_t>(kAntennas1[baseline]);
      const size_t ant2 = static_cast<size_t>(kAntennas2[baseline]);
      for (size_t channel = 0; channel < kNChannels; ++channel) {
        // Sum the model of every direction, weighted by the gains of the two
        // antennas, to form the visibility.
        std::complex<float> vis[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        const std::complex<float> g1 =
            static_cast<std::complex<float>>(gains[ant1][0]);
        const std::complex<float> g2 =
            static_cast<std::complex<float>>(gains[ant2][0]);
        for (size_t d = 0; d < names.size(); ++d) {
          for (size_t p = 0; p < 4; ++p) {
            const std::complex<float> m = rnd();
            unweighted.back()->GetData(names[d])(baseline, channel, p) = m;
            vis[p] += g1 * m * std::conj(g2);
          }
        }
        for (size_t p = 0; p < 4; ++p) {
          unweighted.back()->GetData()(baseline, channel, p) = vis[p];
        }
      }
    }
  }
  dp3::ddecal::AssignAndWeight(unweighted, names, weighted, false, false);
  return weighted;
}

/**
 * The solver a given solve data type is used with in DP3, per SolverFactory.cc:
 * ScalarSolver for the full solve data and IterativeDiagonalSolver for the
 * diagonal (dual) one. Using the real solver is the point of this test, since
 * the aggregation only has to preserve what the solver sees.
 */
template <typename MatrixType>
std::unique_ptr<dp3::ddecal::SolverBase> MakeSolver() {
  // The solvers dispatch over schaapcommon's thread pool via RecursiveFor,
  // which blocks unless the pool has a size. SolverTester's constructor does
  // the same.
  schaapcommon::ThreadPool::GetInstance().SetNThreads(4);
  std::unique_ptr<dp3::ddecal::SolverBase> solver;
  if constexpr (std::is_same_v<MatrixType, aocommon::MC2x2FDiag>) {
    solver = std::make_unique<
        dp3::ddecal::IterativeDiagonalSolver<aocommon::MC2x2FDiag>>();
  } else if constexpr (std::is_same_v<MatrixType, std::complex<float>>) {
    // ScalarSolver only implements the full-solve-data overload; the scalar
    // solve data type needs the iterative scalar solver.
    solver = std::make_unique<
        dp3::ddecal::IterativeScalarSolver<std::complex<float>>>();
  } else {
    solver = std::make_unique<dp3::ddecal::ScalarSolver>();
  }
  solver->Initialize(kNAntennas, std::vector<size_t>(kNDirections, 1),
                     kNChannelBlocks);
  solver->SetMaxIterations(50);
  solver->SetAccuracy(1.0e-8);
  solver->SetStepSize(0.2);
  return solver;
}

/**
 * Reduce a whole interval of weighted buffers and build solve data from it.
 * This is the sequence DDECal itself performs: stream every timestep through
 * a SolveDataAggregator and construct the solve data from the reduced buffers.
 */
template <typename MatrixType>
dp3::ddecal::SolveData<MatrixType> AggregateSolveData(
    const std::vector<DPBuffer>& weighted,
    const std::vector<std::string>& names,
    const std::vector<size_t>& n_solutions,
    dp3::ddecal::SolverDataUse data_use) {
  dp3::ddecal::SolveDataAggregator aggregator(names, kNBaselines, n_solutions,
                                              data_use, kNChannels);
  for (const DPBuffer& buffer : weighted) aggregator.AddTimeSample(buffer);
  return dp3::ddecal::SolveData<MatrixType>(
      aggregator.Finish(), names, kNChannelBlocks, kNAntennas, n_solutions,
      kAntennas1, kAntennas2);
}

/**
 * Solve both the unreduced and the aggregated solve data, and require the
 * per-antenna solutions to agree. Returns the number of rows each system had.
 */
template <typename MatrixType>
std::pair<size_t, size_t> SolveAndCompare(
    const dp3::ddecal::SolveData<MatrixType>& unreduced,
    const dp3::ddecal::SolveData<MatrixType>& reduced) {
  // The solution tensor is indexed [channel_block][antenna * n_solutions *
  // n_solution_pols], so size it the way the solver expects. The number of
  // solution polarizations comes from the solver, not from the data: it is 1
  // for the scalar solvers and 2 for IterativeDiagonalSolver.
  const size_t n_solution_pols =
      MakeSolver<MatrixType>()->NSolutionPolarizations();
  const size_t solution_size = kNAntennas * kNDirections * n_solution_pols;
  std::vector<std::vector<std::complex<double>>> solutions_unreduced(
      kNChannelBlocks);
  std::vector<std::vector<std::complex<double>>> solutions_reduced(
      kNChannelBlocks);
  for (size_t cb = 0; cb < kNChannelBlocks; ++cb) {
    solutions_unreduced[cb].assign(solution_size, 1.0);
    solutions_reduced[cb].assign(solution_size, 1.0);
  }

  MakeSolver<MatrixType>()->Solve(unreduced, solutions_unreduced, 0.0);
  MakeSolver<MatrixType>()->Solve(reduced, solutions_reduced, 0.0);

  for (size_t cb = 0; cb < kNChannelBlocks; ++cb) {
    BOOST_REQUIRE_EQUAL(solutions_unreduced[cb].size(),
                        solutions_reduced[cb].size());
    for (size_t i = 0; i < solutions_unreduced[cb].size(); ++i) {
      const std::complex<double>& a = solutions_unreduced[cb][i];
      const std::complex<double>& b = solutions_reduced[cb][i];
      const bool a_nan = std::isnan(a.real()) || std::isnan(a.imag());
      const bool b_nan = std::isnan(b.real()) || std::isnan(b.imag());
      if (a_nan || b_nan) {
        BOOST_TEST_MESSAGE("cb " << cb << " idx " << i << ": unreduced=" << a
                                 << " reduced=" << b);
        BOOST_CHECK(a_nan && b_nan);
        continue;
      }
      BOOST_CHECK_SMALL(std::abs(a - b), 1.0e-3);
    }
  }
  return {unreduced.ChannelBlock(0).NVisibilities(),
          reduced.ChannelBlock(0).NVisibilities()};
}

}  // namespace

namespace {

/**
 * Run the suite with several threads so the batch-parallel accumulation is
 * actually exercised: RunConstrainedStaticFor only splits the batch axis when
 * more than one thread is available.
 */
struct UseMultipleThreads {
  UseMultipleThreads() {
    schaapcommon::ThreadPool::GetInstance().SetNThreads(4);
  }
};

}  // namespace

BOOST_FIXTURE_TEST_SUITE(aggregate_solve_data, UseMultipleThreads)

BOOST_AUTO_TEST_CASE(incremental_aggregation_matches_unreduced) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);
  const std::vector<DPBuffer> weighted =
      MakeWeightedBuffers(names, n_solutions);

  for (const auto data_use :
       {dp3::ddecal::SolverDataUse::kFull, dp3::ddecal::SolverDataUse::kDual,
        dp3::ddecal::SolverDataUse::kSingle}) {
    const size_t n_pols = data_use == dp3::ddecal::SolverDataUse::kFull   ? 4
                          : data_use == dp3::ddecal::SolverDataUse::kDual ? 2
                                                                          : 1;
    dp3::ddecal::SolveDataAggregator aggregator(names, kNBaselines, n_solutions,
                                                data_use, kNChannels);
    for (const DPBuffer& buffer : weighted) aggregator.AddTimeSample(buffer);
    const std::vector<DPBuffer> reduced = aggregator.Finish();
    BOOST_CHECK_EQUAL(reduced.size(), kNDirections * n_pols);
    BOOST_CHECK_THROW(aggregator.Finish(), std::runtime_error);
    BOOST_CHECK_THROW(aggregator.AddTimeSample(weighted.front()),
                      std::runtime_error);
  }

  // Every datause has to reproduce the unreduced solve, not just kFull.
  // kSingle is the case that a fold of XX alone would get wrong: the ordinary
  // scalar solve uses 0.5 * (XX + YY) for the data and the model alike.
  SolveAndCompare(
      FullSolveData(std::vector<DPBuffer>(weighted), names, kNChannelBlocks,
                    kNAntennas, n_solutions, kAntennas1, kAntennas2),
      AggregateSolveData<aocommon::MC2x2F>(weighted, names, n_solutions,
                                           SolverDataUse::kFull));
  SolveAndCompare(
      DuoSolveData(std::vector<DPBuffer>(weighted), names, kNChannelBlocks,
                   kNAntennas, n_solutions, kAntennas1, kAntennas2),
      AggregateSolveData<aocommon::MC2x2FDiag>(weighted, names, n_solutions,
                                               SolverDataUse::kDual));
  SolveAndCompare(
      UniSolveData(std::vector<DPBuffer>(weighted), names, kNChannelBlocks,
                   kNAntennas, n_solutions, kAntennas1, kAntennas2),
      AggregateSolveData<std::complex<float>>(weighted, names, n_solutions,
                                              SolverDataUse::kSingle));
}

BOOST_AUTO_TEST_CASE(incremental_aggregation_rejects_channels_at_finish) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);
  const std::vector<DPBuffer> weighted =
      MakeWeightedBuffers(names, n_solutions);
  const auto make_aggregator = [&]() {
    return std::make_unique<dp3::ddecal::SolveDataAggregator>(
        names, kNBaselines, n_solutions, dp3::ddecal::SolverDataUse::kSingle,
        kNChannels);
  };
  auto unmasked = make_aggregator();
  auto masked = make_aggregator();
  for (const DPBuffer& buffer : weighted) {
    unmasked->AddTimeSample(buffer);
    masked->AddTimeSample(buffer);
  }
  const std::vector<DPBuffer> expected = unmasked->Finish();
  BOOST_CHECK_THROW(masked->Finish({true}), std::runtime_error);
  const std::vector<DPBuffer> actual =
      masked->Finish({false, true, true, false});
  BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
  const std::vector<std::string> data_names{"", names[0], names[1]};
  for (size_t row = 0; row < actual.size(); ++row) {
    for (size_t baseline = 0; baseline < kNBaselines; ++baseline) {
      for (size_t channel = 0; channel < kNChannels; ++channel) {
        for (const std::string& name : data_names) {
          for (size_t pol = 0; pol < kNPolarizations; ++pol) {
            const auto value =
                actual[row].GetData(name)(baseline, channel, pol);
            if (channel == 1 || channel == 2) {
              BOOST_CHECK_EQUAL(value, std::complex<float>(0.0f));
            } else {
              BOOST_CHECK_EQUAL(
                  value, expected[row].GetData(name)(baseline, channel, pol));
            }
          }
        }
      }
    }
  }
}

/**
 * The reduced (aggregated) solve data must produce the same per-antenna
 * solutions as the unreduced solve data, for a solution interval of several
 * timesteps. Both are fed to the same solver; only the shape of the data
 * differs.
 *
 * The number of folded polarizations follows ddecal.datause, so this is checked
 * for both the full (4 correlations) and the dual (XX and YY) solve data. The
 * dual case additionally covers the correlation remapping: ToSpecificMatrix
 * reads data[0] and data[3], and the aggregation has to place the folded values
 * at exactly those indices.
 */
BOOST_AUTO_TEST_CASE(aggregated_matches_unreduced) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);
  const std::vector<DPBuffer> buffers = MakeWeightedBuffers(names, n_solutions);

  // Unreduced: the ordinary per-visibility solve data.
  const FullSolveData unreduced_full(std::vector<DPBuffer>(buffers), names,
                                     kNChannelBlocks, kNAntennas, n_solutions,
                                     kAntennas1, kAntennas2);
  const dp3::ddecal::DuoSolveData unreduced_dual(
      std::vector<DPBuffer>(buffers), names, kNChannelBlocks, kNAntennas,
      n_solutions, kAntennas1, kAntennas2);

  const FullSolveData reduced_full = AggregateSolveData<aocommon::MC2x2F>(
      buffers, names, n_solutions, dp3::ddecal::SolverDataUse::kFull);
  const dp3::ddecal::DuoSolveData reduced_dual =
      AggregateSolveData<aocommon::MC2x2FDiag>(
          buffers, names, n_solutions, dp3::ddecal::SolverDataUse::kDual);

  BOOST_TEST_MESSAGE("unreduced(full): n_directions = "
                     << unreduced_full.ChannelBlock(0).NDirections() << ", "
                     << "n_visibilities = "
                     << unreduced_full.ChannelBlock(0).NVisibilities());
  BOOST_TEST_MESSAGE("reduced(full): n_visibilities = "
                     << reduced_full.ChannelBlock(0).NVisibilities());
  BOOST_TEST_MESSAGE("reduced(dual): n_visibilities = "
                     << reduced_dual.ChannelBlock(0).NVisibilities());

  // Both systems must actually be solvable, i.e. have full column rank.
  BOOST_CHECK_GT(unreduced_full.ChannelBlock(0).NVisibilities(), 0u);
  BOOST_CHECK_GT(reduced_full.ChannelBlock(0).NVisibilities(), 0u);
  BOOST_CHECK_GT(reduced_dual.ChannelBlock(0).NVisibilities(), 0u);

  BOOST_TEST_MESSAGE("solving unreduced, reduced(full), reduced(dual)...");
  SolveAndCompare(unreduced_full, reduced_full);
  SolveAndCompare(unreduced_dual, reduced_dual);
}

BOOST_AUTO_TEST_CASE(aggregation_reduces_row_count) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);
  const std::vector<DPBuffer> buffers = MakeWeightedBuffers(names, n_solutions);

  const FullSolveData unreduced(std::vector<DPBuffer>(buffers), names,
                                kNChannelBlocks, kNAntennas, n_solutions,
                                kAntennas1, kAntennas2);
  const FullSolveData reduced = AggregateSolveData<aocommon::MC2x2F>(
      buffers, names, n_solutions, dp3::ddecal::SolverDataUse::kFull);

  // The interval axis went from kNTimesteps entries to lin = ndir * npol.
  // SolveData turns that into rows per channel block. Baseline 2 is an
  // autocorrelation, which SolveData skips, so only kNBaselines - 1 baselines
  // contribute rows.
  const size_t unreduced_rows = unreduced.ChannelBlock(0).NVisibilities();
  const size_t reduced_rows = reduced.ChannelBlock(0).NVisibilities();
  const size_t lin = kNDirections * kNPolarizations;
  const size_t channels_per_block = kNChannels / kNChannelBlocks;
  const size_t usable_baselines = kNBaselines - 1;  // minus the autocorrelation

  BOOST_TEST_MESSAGE("rows: unreduced = " << unreduced_rows << ", reduced = "
                                          << reduced_rows << ", lin = " << lin);
  BOOST_CHECK_EQUAL(unreduced_rows,
                    kNTimesteps * usable_baselines * channels_per_block);
  BOOST_CHECK_EQUAL(reduced_rows, lin * usable_baselines * channels_per_block);

  // The reduced data must be finite and non-degenerate: a singular or
  // non-finite reduced system cannot be solved at all.
  for (size_t cb = 0; cb < kNChannelBlocks; ++cb) {
    const auto& block = reduced.ChannelBlock(cb);
    for (size_t v = 0; v < block.NVisibilities(); ++v) {
      const aocommon::MC2x2F vis = block.Visibility(v);
      for (size_t p = 0; p < 4; ++p) {
        BOOST_CHECK(std::isfinite(vis.Get(p).real()));
        BOOST_CHECK(std::isfinite(vis.Get(p).imag()));
      }
      for (size_t d = 0; d < block.NDirections(); ++d) {
        const aocommon::MC2x2F model = block.ModelVisibility(d, v);
        for (size_t p = 0; p < 4; ++p) {
          BOOST_CHECK(std::isfinite(model.Get(p).real()));
          BOOST_CHECK(std::isfinite(model.Get(p).imag()));
        }
      }
    }
  }
}

/**
 * The point of the aggregation: peak model memory must not grow with the
 * solution interval. The reduced row count is lin * baselines * channels, so it
 * is the same for a short and for a long interval, while the unreduced one
 * scales with the number of timesteps. Folding twice as many timesteps must
 * therefore leave the reduced size untouched.
 *
 * The data has to stay solvable, so only long enough intervals are compared.
 */
BOOST_AUTO_TEST_CASE(aggregation_size_independent_of_interval) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);

  size_t reduced_rows = 0;
  size_t reduced_model_values = 0;
  size_t unreduced_rows = 0;
  // Fold intervals of 8 and 24 timesteps: three times as many samples must not
  // change the reduced size at all.
  for (const size_t n_timesteps : {kNTimesteps, 3 * kNTimesteps}) {
    const std::vector<DPBuffer> buffers =
        MakeWeightedBuffers(names, n_solutions, n_timesteps);

    const FullSolveData unreduced(std::vector<DPBuffer>(buffers), names,
                                  kNChannelBlocks, kNAntennas, n_solutions,
                                  kAntennas1, kAntennas2);
    const FullSolveData reduced = AggregateSolveData<aocommon::MC2x2F>(
        buffers, names, n_solutions, dp3::ddecal::SolverDataUse::kFull);

    BOOST_TEST_MESSAGE("n_timesteps = "
                       << n_timesteps << ": unreduced rows = "
                       << unreduced.ChannelBlock(0).NVisibilities()
                       << ", reduced rows = "
                       << reduced.ChannelBlock(0).NVisibilities());

    // The unreduced row count does grow with the interval, by construction.
    BOOST_CHECK_EQUAL(
        unreduced.ChannelBlock(0).NVisibilities(),
        n_timesteps * (kNBaselines - 1) * (kNChannels / kNChannelBlocks));

    // The model values held by the reduced data are what bounds the peak.
    const size_t model_values = reduced.ChannelBlock(0).NVisibilities() *
                                reduced.ChannelBlock(0).NDirections() *
                                kNPolarizations;

    if (reduced_rows == 0) {
      reduced_rows = reduced.ChannelBlock(0).NVisibilities();
      reduced_model_values = model_values;
      unreduced_rows = unreduced.ChannelBlock(0).NVisibilities();
    } else {
      // Tripling the interval left the reduced size untouched, while the
      // unreduced size tripled. That is the whole point of the aggregation.
      BOOST_CHECK_EQUAL(unreduced.ChannelBlock(0).NVisibilities(),
                        3 * unreduced_rows);
      BOOST_CHECK_EQUAL(reduced.ChannelBlock(0).NVisibilities(), reduced_rows);
      BOOST_CHECK_EQUAL(model_values, reduced_model_values);
    }
  }
}

BOOST_AUTO_TEST_CASE(aggregation_rejects_sub_solutions) {
  const std::vector<std::string> names = DirectionNames();
  const std::vector<size_t> n_solutions(kNDirections, 1);
  const std::vector<DPBuffer> buffers = MakeWeightedBuffers(names, n_solutions);

  const std::vector<size_t> split(kNDirections, 2);
  BOOST_CHECK_THROW(
      AggregateSolveData<aocommon::MC2x2F>(buffers, names, split,
                                           dp3::ddecal::SolverDataUse::kFull),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
