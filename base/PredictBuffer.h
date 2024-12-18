// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PREDICTBUFFER_H_
#define DP3_BASE_PREDICTBUFFER_H_

#include <cassert>
#include <complex>
#include <vector>

#include <aocommon/matrix2x2.h>

#include <EveryBeam/station.h>
#include <EveryBeam/telescope/telescope.h>

namespace dp3 {
namespace base {

class PredictBuffer {
 public:
  void Resize(size_t n_directions, size_t n_channels, size_t n_stations,
              bool full_beam) {
    // The full buffer is not used when Stokes I is used -- conditionally
    // allocating full/scalar will save some memory.
    if (full_beam) {
      full_beam_values_.resize(n_directions);
      for (std::vector<aocommon::MC2x2>& values : full_beam_values_)
        values.resize(n_stations * n_channels);
    } else {
      scalar_beam_values_.resize(n_directions);
      for (std::vector<std::complex<double>>& values : scalar_beam_values_)
        values.resize(n_stations * n_channels);
    }

    full_beam_ = full_beam;
    n_stations_ = n_stations;
  }

  aocommon::MC2x2* GetFullBeamValues(size_t direction_index) {
    assert(full_beam_);
    return full_beam_values_[direction_index].data();
  }

  std::complex<double>* GetScalarBeamValues(size_t direction_index) {
    assert(!full_beam_);
    return scalar_beam_values_[direction_index].data();
  }

  size_t NStations() const { return n_stations_; }

 private:
  std::vector<std::vector<aocommon::MC2x2>> full_beam_values_;
  std::vector<std::vector<std::complex<double>>> scalar_beam_values_;
  size_t n_stations_ = 0;
  bool full_beam_ = false;
};

}  // namespace base
}  // namespace dp3

#endif
