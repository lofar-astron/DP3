// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDE_REGULAR_SOLVER_BASE_H
#define DDE_REGULAR_SOLVER_BASE_H

#include "SolverBase.h"

namespace dp3 {
namespace ddecal {

class SolverBuffer;

class RegularSolverBase : public SolverBase {
 public:
  RegularSolverBase() : n_channels_(0) {}

  /**
   * Prepares the solver with the given dimensionality info
   * and antenna mapping.
   * The antenna arrays map the data provided in @Solve to the antennas.
   */
  virtual void Initialize(size_t nAntennas, size_t nDirections,
                          size_t nChannels, size_t nChannelBlocks,
                          const std::vector<int>& ant1,
                          const std::vector<int>& ant2);

  /**
   * Solves multi-directional Jones matrices. Takes the (single) measured data
   * and the (multi-directional) model data, and solves the optimization
   * problem that minimizes the norm of the differences.
   *
   * @param solver_buffer Buffer with unweighted data, weights and model data.
   * @param solutions The per-channel and per-antenna solutions.
   * solutions[ch] is a pointer for channelblock ch to antenna x directions x
   * pol solutions.
   * @param statStream Optional pointer to a stream for displaying statistics.
   */
  virtual SolveResult Solve(const SolverBuffer& solver_buffer,
                            std::vector<std::vector<DComplex>>& solutions,
                            double time, std::ostream* statStream) = 0;

 protected:
  int AntennaIndex1(size_t baseline) const { return ant1_[baseline]; }
  int AntennaIndex2(size_t baseline) const { return ant2_[baseline]; }

  /**
   * @param block A channel block index, less than NChannelBlocks().
   * @return The index of the first channel in the given channel block.
   */
  size_t FirstChannel(size_t block) const {
    return block * n_channels_ / n_channel_blocks_;
  }

  size_t NBaselines() const { return ant1_.size(); }
  size_t NChannels() const { return n_channels_; }

  std::vector<int> ant1_, ant2_;
  size_t n_channels_;
};

}  // namespace ddecal
}  // namespace dp3

#endif
