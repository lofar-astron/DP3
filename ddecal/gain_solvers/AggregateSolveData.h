// See the LICENSE file at the top-level directory of this distribution.

#ifndef DDECAL_AGGREGATESOLVEDATA_H_
#define DDECAL_AGGREGATESOLVEDATA_H_

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "base/DPBuffer.h"
#include "ddecal/Settings.h"

namespace dp3::ddecal {

class CalibrationAggregation;

/**
 * Reduce the weighted model data of a solution interval with the
 * three-phase calibration aggregation, and hand the result to SolveData.
 *
 * This exists to reduce memory, not to increase speed: SolveData would
 * otherwise hold one row per (timestep, baseline, channel) of the interval,
 * so its peak model memory would grow with the solution interval length.
 * The reduction is streamed, one timestep per AddTimeSample() call, so the
 * timestep axis is never held in memory all at once.
 *
 * The aggregation is ska-sdp-func's. See CalibrationAggregation for how the
 * reduction works and what it assumes about the data.
 */
class SolveDataAggregator {
 public:
  /**
   * The baselines are only ever counted, never inspected: the aggregation
   * treats a baseline as an opaque batch id and does not need to know which
   * antennas it connects, nor how many antennas are in use. The per-antenna
   * bookkeeping is done by SolveData, which is given the antenna indices.
   *
   * @param direction_names Model data column names, one per direction.
   * @param n_baselines Number of baselines in the data to be accumulated.
   * @param n_solutions_per_direction Must be all ones; the aggregation folds
   *        the whole interval and cannot be partitioned per sub-solution.
   * @param data_use Number of polarizations to fold, see ddecal.datause.
   * @param n_channels Number of channels in the data to be accumulated.
   */
  SolveDataAggregator(const std::vector<std::string>& direction_names,
                      size_t n_baselines,
                      const std::vector<size_t>& n_solutions_per_direction,
                      SolverDataUse data_use, size_t n_channels);
  ~SolveDataAggregator();

  SolveDataAggregator(const SolveDataAggregator&) = delete;
  SolveDataAggregator& operator=(const SolveDataAggregator&) = delete;

  /**
   * Add one time sample: the weighted data and model data of a single timestep,
   * as produced by ddecal::AssignAndWeight. Must be called in timestep order,
   * and only before Finish().
   * @param weighted_buffer Weighted data and model data for one time sample.
   */
  void AddTimeSample(const base::DPBuffer& weighted_buffer);

  /**
   * Reduce the accumulated timesteps and return the reduced systems, one
   * DPBuffer per reduced row of the linear system.
   * @param rejected_channels Channels to leave out of the reduction. Empty
   *        means that no channel is rejected. Can only be given once, and not
   *        together with an earlier call to Finish().
   */
  std::vector<base::DPBuffer> Finish(
      const std::vector<bool>& rejected_channels = {});

  /**
   * Wall-clock breakdown of the last Finish() call, in seconds. Reported back
   * to DDECal so the timing printout can show where finalization time goes.
   * The parallel phases include their thread-barrier overhead, since that is
   * part of what they cost in wall time.
   */
  struct FinishTimings {
    double identity = 0.0;
    double square_root = 0.0;
    double apply = 0.0;
    double check = 0.0;
    double copy = 0.0;
  };

  const FinishTimings& last_finish_timings() const {
    return last_finish_timings_;
  }

  /// Number of reduced systems, i.e. baselines x channels.
  size_t n_batches() const { return n_batches_; }

  /// Size (lin_) of the reduced systems, i.e. n_directions x n_polarizations.
  size_t linear_size() const { return lin_; }

 private:
  /**
   * Which of DP3's four correlations each datause folds, and where the folded
   * value must be placed in the DPBuffer that SolveData later reads.
   *
   * There is one entry per folded polarization. A folded polarization is the
   * mean of the source correlations listed for it, and is copied to every
   * destination correlation listed for it.
   *
   * The source indices are DP3's native [XX, XY, YX, YY] order. The
   * destination indices differ from them: ToSpecificMatrix in SolveData.cc
   * reads
   *   - MC2x2F      -> all four correlations,
   *   - MC2x2FDiag  -> data[0] and data[3],
   *   - complex     -> 0.5 * (data[0] + data[3]),
   * so the diagonal pair must land at indices 0 and 3, and a single folded
   * polarization must land at both, so that the scalar reads its own value
   * back. Because the same conversion is applied to the model data, kSingle
   * has to fold the average of the diagonal pair too: averaging only XX would
   * make the aggregation solve a different system than the ordinary scalar
   * solve whenever XX and YY differ.
   */
  struct PolarizationLayout {
    std::vector<std::vector<size_t>> source;
    std::vector<std::vector<size_t>> destination;
  };

  static PolarizationLayout GetLayout(SolverDataUse data_use);

  std::vector<std::string> direction_names_;
  PolarizationLayout layout_;
  size_t n_channels_;
  size_t n_directions_;
  size_t n_baselines_;
  size_t n_batches_;
  size_t n_pols_;
  size_t lin_;
  /// Staging area for one timestep, laid out as Accumulate expects it:
  /// staged_model_ is [baseline, channel, direction, polarization] and
  /// staged_vis_ is [baseline, channel, polarization]. Their leading two axes
  /// flatten to the batch axis, so both drop straight into the rank-4 / rank-3
  /// shape the aggregation takes.
  xt::xtensor<std::complex<float>, 4> staged_model_;
  xt::xtensor<std::complex<float>, 3> staged_vis_;
  std::vector<size_t> vis_shape_;
  std::vector<size_t> model_shape_;
  /// Timesteps added over the whole interval, for the log and the empty check.
  size_t n_timesteps_ = 0;
  bool finished_ = false;
  FinishTimings last_finish_timings_;

  /**
   * The three-phase aggregation, which is the only part of this class that
   * needs the ska-sdp-func headers. Holding it behind a pointer is what keeps
   * those headers out of this file, and therefore out of DDECal.h.
   */
  std::unique_ptr<CalibrationAggregation> aggregation_;
};

}  // namespace dp3::ddecal

#endif  // DDECAL_AGGREGATESOLVEDATA_H_
