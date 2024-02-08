// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PREDICTBUFFER_H_
#define DP3_BASE_PREDICTBUFFER_H_

#include <cassert>
#include <complex>
#include <vector>

#include <EveryBeam/station.h>
#include <EveryBeam/telescope/telescope.h>

namespace dp3 {
namespace base {

class PredictBuffer {
 public:
  void resize(size_t n_threads, size_t n_correlations, size_t n_channels,
              size_t n_baselines, size_t n_stations, bool include_beam,
              bool full_beam) {
    if (include_beam) {
      // The full buffer is not used when Stokes I is used -- conditionally
      // allocating full/scalar will save some memory.
      if (full_beam) {
        full_beam_values_.resize(n_threads);
      } else {
        scalar_beam_values_.resize(n_threads);
      }

      for (size_t i = 0; i != n_threads; ++i) {
        if (full_beam) {
          full_beam_values_[i].resize(n_stations * n_channels);
        } else {
          scalar_beam_values_[i].resize(n_stations * n_channels);
        }
      }
    }

    full_beam_ = full_beam;
  }

  std::vector<aocommon::MC2x2>& GetFullBeamValues(size_t threadIndex) {
    assert(full_beam_);
    return full_beam_values_[threadIndex];
  }

  std::vector<everybeam::complex_t>& GetScalarBeamValues(size_t threadIndex) {
    assert(!full_beam_);
    return scalar_beam_values_[threadIndex];
  }

 private:
  std::vector<std::vector<aocommon::MC2x2>> full_beam_values_;
  std::vector<std::vector<everybeam::complex_t>> scalar_beam_values_;
  bool full_beam_{false};
};

}  // namespace base
}  // namespace dp3

#endif
