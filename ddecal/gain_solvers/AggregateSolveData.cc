// See the LICENSE file at the top-level directory of this distribution.

#include "AggregateSolveData.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <aocommon/logger.h>
#include <schaapcommon/threading/staticfor.h>
#include <xtensor/core/xnoalias.hpp>
#include <xtensor/views/xview.hpp>

#include "ddecal/CalibrationAggregation.h"

namespace dp3::ddecal {

namespace {

constexpr size_t kNCorrelations = 4;

/**
 * Fold one polarization column of a [rows, correlation] source block into the
 * destination column `dest`: the mean of the DP3 correlations it is built from.
 * kFull and kDual list a single correlation each, so this copies it; kSingle
 * lists the diagonal pair, which is what the ordinary scalar solve averages.
 * The layouts only ever list one or two correlations.
 *
 * Written as an xtensor expression so the copy and the sum operate on whole
 * columns at once; `dest` is a view into the staging buffer. xt::noalias keeps
 * the expression from materialising an intermediate buffer per folded column
 * (the build defines XTENSOR_FORCE_TEMPORARY_MEMORY_IN_ASSIGNMENTS).
 */
template <class Dest, class Source>
void FoldPolarization(Dest&& dest, const Source& source,
                      const std::vector<size_t>& correlations) {
  assert(!correlations.empty() && correlations.size() <= 2);
  if (correlations.size() == 1) {
    xt::noalias(dest) = xt::view(source, xt::all(), correlations[0]);
  } else {
    xt::noalias(dest) = (xt::view(source, xt::all(), correlations[0]) +
                         xt::view(source, xt::all(), correlations[1])) *
                        (1.0f / static_cast<float>(correlations.size()));
  }
}

}  // namespace

SolveDataAggregator::PolarizationLayout SolveDataAggregator::GetLayout(
    SolverDataUse data_use) {
  switch (data_use) {
    case SolverDataUse::kFull:
      return {{{0}, {1}, {2}, {3}}, {{0}, {1}, {2}, {3}}};
    case SolverDataUse::kDual:
      return {{{0}, {3}}, {{0}, {3}}};
    case SolverDataUse::kSingle:
      return {{{0, 3}}, {{0, 3}}};
  }
  throw std::runtime_error(
      "Unknown ddecal datause for calibration aggregation.");
}

SolveDataAggregator::SolveDataAggregator(
    const std::vector<std::string>& direction_names, size_t n_baselines,
    const std::vector<size_t>& n_solutions_per_direction,
    SolverDataUse data_use, size_t n_channels)
    : direction_names_(direction_names),
      layout_(GetLayout(data_use)),
      n_channels_(n_channels),
      n_directions_(direction_names.size()),
      n_baselines_(n_baselines),
      n_batches_(n_baselines * n_channels),
      n_pols_(layout_.source.size()),
      lin_(n_directions_ * n_pols_),
      staged_model_(xt::xtensor<std::complex<float>, 4>::from_shape(
          {n_baselines_, n_channels_, n_directions_, n_pols_})),
      staged_vis_(xt::xtensor<std::complex<float>, 3>::from_shape(
          {n_baselines_, n_channels_, n_pols_})),
      vis_shape_{n_batches_, 1, n_pols_},
      model_shape_{n_batches_, 1, n_directions_, n_pols_},
      aggregation_(std::make_unique<CalibrationAggregation>(
          n_batches_, n_directions_, n_pols_)) {
  if (n_directions_ == 0) {
    throw std::runtime_error(
        "Cannot aggregate model data without any direction.");
  }
  if (n_baselines_ == 0) {
    throw std::runtime_error("Cannot aggregate model data without baselines.");
  }
  for (const size_t n : n_solutions_per_direction) {
    if (n != 1) {
      throw std::runtime_error(
          "ddecal.aggregate=true requires all entries of "
          "solutions_per_direction to be 1, but found " +
          std::to_string(n) + ".");
    }
  }

  aggregation_->Reset();
}

SolveDataAggregator::~SolveDataAggregator() = default;

void SolveDataAggregator::AddTimeSample(const base::DPBuffer& buffer) {
  if (finished_) {
    throw std::runtime_error(
        "Cannot accumulate after aggregation is finished.");
  }
  const base::DPBuffer::DataType& reference = buffer.GetData();
  if (reference.shape().size() != 3) {
    throw std::runtime_error(
        "Cannot aggregate model data: expected buffers of shape "
        "[n_baselines, n_channels, n_correlations], got " +
        std::to_string(reference.shape().size()) + " dimensions.");
  }
  if (reference.shape(2) != kNCorrelations) {
    throw std::runtime_error("Cannot aggregate model data: expected " +
                             std::to_string(kNCorrelations) +
                             " correlations, got " +
                             std::to_string(reference.shape(2)) + ".");
  }
  if (reference.shape(0) != n_baselines_) {
    throw std::runtime_error("Cannot aggregate model data: the buffers hold " +
                             std::to_string(reference.shape(0)) +
                             " baselines but " + std::to_string(n_baselines_) +
                             " were expected.");
  }
  if (reference.shape(1) != n_channels_) {
    throw std::runtime_error(
        "Cannot aggregate model data: channel count changed within an "
        "interval.");
  }

  const base::DPBuffer::DataType& data = buffer.GetData();

  // Hoist the per-direction lookups out of the loops below: GetData resolves a
  // column name, so doing it inside the batch loop would repeat that lookup
  // once per folded value instead of once per timestep.
  std::vector<const base::DPBuffer::DataType*> model_data(n_directions_);
  for (size_t direction = 0; direction < n_directions_; ++direction) {
    model_data[direction] = &buffer.GetData(direction_names_[direction]);
  }

  // A model or visibility column is [n_baselines, n_channels, n_correlations].
  // Folding a polarization is a view of the correlation axis plus a sum, and
  // each direction is copied into its own slice of the staged model. The
  // baselines are staged in parallel: every baseline writes a disjoint slice
  // and only reads the input.
  schaapcommon::RunStaticFor<size_t>(
      0, n_baselines_, [&](size_t baseline_begin, size_t baseline_end) {
        for (size_t baseline = baseline_begin; baseline != baseline_end;
             ++baseline) {
          const auto vis = xt::view(data, baseline, xt::all(), xt::all());
          for (size_t p = 0; p != n_pols_; ++p) {
            FoldPolarization(xt::view(staged_vis_, baseline, xt::all(), p), vis,
                             layout_.source[p]);
          }
          for (size_t direction = 0; direction != n_directions_; ++direction) {
            const auto model = xt::view(*model_data[direction], baseline,
                                        xt::all(), xt::all());
            for (size_t p = 0; p != n_pols_; ++p) {
              FoldPolarization(
                  xt::view(staged_model_, baseline, xt::all(), direction, p),
                  model, layout_.source[p]);
            }
          }
        }
      });

  // The staging buffer of this single timestep is contiguous, so the
  // aggregation folds it with its outer-product kernel.
  aggregation_->Accumulate(staged_model_.data(), model_shape_,
                           {lin_, lin_, n_pols_, 1}, staged_vis_.data(),
                           vis_shape_, {n_pols_, n_pols_, 1});
  ++n_timesteps_;
}

std::vector<base::DPBuffer> SolveDataAggregator::Finish(
    const std::vector<bool>& rejected_channels) {
  if (finished_ || n_timesteps_ == 0) {
    throw std::runtime_error(
        "Cannot finish aggregation without new timesteps.");
  }
  if (!rejected_channels.empty() && rejected_channels.size() != n_channels_) {
    throw std::runtime_error(
        "Cannot finish aggregation: rejected channel count differs from "
        "input.");
  }
  finished_ = true;
  CalibrationAggregation& accumulator = *aggregation_;

  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
  };
  last_finish_timings_ = FinishTimings{};

  Clock::time_point phase_start = Clock::now();
  for (size_t batch = 0; batch < n_batches_; ++batch) {
    if (!rejected_channels.empty() && rejected_channels[batch % n_channels_]) {
      accumulator.ReplaceWithIdentity(batch);
    }
  }
  last_finish_timings_.identity = elapsed(phase_start);

  phase_start = Clock::now();
  accumulator.SquareRoot();
  last_finish_timings_.square_root = elapsed(phase_start);

  phase_start = Clock::now();
  accumulator.Apply();
  last_finish_timings_.apply = elapsed(phase_start);

  aocommon::Logger::Debug << "  Aggregated " << n_timesteps_
                          << " timestep(s) into " << n_batches_
                          << " reduced system(s) of size " << lin_ << "x"
                          << lin_ << ".\n";

  // sdp_cal_aggregate_square_root does not divide by a zero eigenvalue: it
  // thresholds them against max_eigenvalue * linear_size * epsilon and
  // clamps the reciprocal to zero, so a singular Gram matrix yields a
  // rank-deficient but entirely finite reduced system.
  // Such a batch simply contributes nothing to the per-antenna normal
  // equations, which is the right answer for a batch that carries no
  // information, so a flagged baseline or a fully flagged channel block
  // degrades gracefully instead of aborting the interval.
  //
  // Kept because it turns a silent all-NaN solution table into a diagnosable
  // message if that clamping is ever dropped upstream. Should that happen,
  // the message below becomes reachable, and its advice to switch datause or
  // disable aggregation would finally apply.
  {
    // Parallel over batches: the batches are independent, and this scan is
    // O(batches x lin^2) serial work on a box that is otherwise idle here.
    std::atomic<bool> non_finite{false};
    schaapcommon::StaticFor<size_t> check_loop;
    phase_start = Clock::now();
    check_loop.Run(0, n_batches_, [&](size_t begin, size_t end) {
      for (size_t batch = begin; batch < end && !non_finite; ++batch) {
        for (size_t r = 0; r < lin_ && !non_finite; ++r) {
          for (size_t p = 0; p < n_pols_; ++p) {
            const auto v = accumulator.ReducedVisibility()(batch, r, p);
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
              non_finite = true;
              break;
            }
          }
        }
      }
    });
    last_finish_timings_.check = elapsed(phase_start);
    if (non_finite) {
      throw std::runtime_error(
          "The aggregated model system is not finite, which means its Gram "
          "matrix M^H M could not be reduced. This should not happen: the "
          "Hermitian square root clamps eigenvalues below "
          "max_eigenvalue * n * epsilon to zero rather than taking their "
          "reciprocal, so a singular M^H M yields a finite but rank-deficient "
          "reduced system. M^H M is rank deficient when the model data "
          "folded over the solution interval does not span all "
          "n_directions * n_polarizations = " +
          std::to_string(lin_) +
          " dimensions, which typically happens when the model is nearly "
          "collinear across the interval, e.g. a predictable source whose "
          "model barely changes between timesteps, or linearly polarized data "
          "folded with all four correlations. Increase ddecal.solint, use "
          "ddecal.datause=dual, or disable ddecal.aggregate.");
    }
  }

  // Copy the reduced systems out into DPBuffers that SolveData can consume.
  //
  // The trick is the buffer vector itself. SolveData takes
  // std::vector<DPBuffer>, and its index is the solution-interval axis: one
  // buffer per timestep of the interval. The aggregation replaces those
  // timesteps with the lin rows of the reduced system, so the vector has lin
  // entries instead of solint, and buffer r carries reduced row r.
  //
  // Every buffer keeps the ordinary [n_baselines, n_channels, n_correlations]
  // shape, so SolveData walks its normal baseline x channel grid and gives each
  // row the antennas of its baseline. That is what makes the per-antenna
  // accumulation work on the reduced data: a reduced row still belongs to the
  // batch it came from, and the solver sums exactly the rows involving each
  // antenna. Because Q is Hermitian, Q^H Q = M^H M and Q^H b = M^H V, so those
  // sums reproduce the unreduced per-antenna normal equations.
  //
  // The model of direction d in buffer r is the 2x2 Q(batch, r, d, :), and the
  // visibility is b(batch, r, :).
  const xt::xtensor<std::complex<float>, 4>& q = accumulator.ReducedModel();
  const xt::xtensor<std::complex<float>, 3>& b =
      accumulator.ReducedVisibility();

  // Each buffer r is built completely independently (its own allocations,
  // fills and element writes into distinct arrays), so the loop is parallel
  // over r. For a wide solution this copy is tens of GB of zero-fills,
  // allocations and element writes, which would otherwise all run on one core
  // while the rest of the machine idles.
  std::vector<base::DPBuffer> buffers(lin_);
  schaapcommon::StaticFor<size_t> copy_loop;
  phase_start = Clock::now();
  copy_loop.Run(0, lin_, [&](size_t r_begin, size_t r_end) {
    for (size_t r = r_begin; r < r_end; ++r) {
      base::DPBuffer& buffer = buffers[r];
      buffer.GetData().resize({n_baselines_, n_channels_, kNCorrelations});
      // The reduced rows carry the weighting that AssignAndWeight already
      // folded into the data, so they are unit weight. Set them explicitly
      // rather than leaving the buffer without a weights array, because
      // SolveData treats a non-empty weights array as "this solve data has
      // weights" and other solvers and GetSolutionWeights() read it.
      buffer.GetWeights().resize({n_baselines_, n_channels_, kNCorrelations});
      buffer.GetWeights().fill(1.0f);

      for (size_t channel = 0; channel < n_channels_; ++channel) {
        if (rejected_channels.empty() || !rejected_channels[channel]) continue;
        xt::view(buffer.GetWeights(), xt::all(), channel, xt::all()).fill(0.0f);
      }
      buffer.GetData().fill(0.0f);
      for (const std::string& name : direction_names_) {
        buffer.AddData(name);
        buffer.GetData(name).fill(0.0f);
      }
      for (size_t baseline = 0; baseline < n_baselines_; ++baseline) {
        for (size_t channel = 0; channel < n_channels_; ++channel) {
          const size_t batch = baseline * n_channels_ + channel;
          if (!rejected_channels.empty() && rejected_channels[channel])
            continue;
          base::DPBuffer::DataType& data = buffer.GetData();
          for (size_t p = 0; p < n_pols_; ++p) {
            for (const size_t correlation : layout_.destination[p]) {
              data(baseline, channel, correlation) = b(batch, r, p);
            }
          }
          for (size_t d = 0; d < n_directions_; ++d) {
            base::DPBuffer::DataType& model_data =
                buffer.GetData(direction_names_[d]);
            for (size_t p = 0; p < n_pols_; ++p) {
              for (const size_t correlation : layout_.destination[p]) {
                model_data(baseline, channel, correlation) = q(batch, r, d, p);
              }
            }
          }
        }
      }
    }
  });
  last_finish_timings_.copy = elapsed(phase_start);

  return buffers;
}

}  // namespace dp3::ddecal
