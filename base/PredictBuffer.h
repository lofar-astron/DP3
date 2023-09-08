// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PREDICTBUFFER_H_
#define DP3_BASE_PREDICTBUFFER_H_

#include <complex>
#include <vector>

#include <casacore/casa/Arrays/Cube.h>

#include <xtensor/xtensor.hpp>

#include <EveryBeam/station.h>
#include <EveryBeam/telescope/telescope.h>

#include "common/MatrixComplexDouble2x2.h"

namespace dp3 {
namespace base {

class PredictBuffer {
 public:
  typedef std::complex<double> dcomplex;

  void resize(size_t n_threads, size_t n_correlations, size_t n_channels,
              size_t n_baselines, size_t n_stations, bool include_beam) {
    model_visibilities_.resize(n_threads);
    for (size_t i = 0; i != n_threads; ++i) {
      model_visibilities_[i].resize({n_baselines, n_channels, n_correlations});
    }

    if (include_beam) {
      patch_model_visibilities_.resize(n_threads);
      // TODO the full buffer is not used when Stokes I is used -- conditionally
      // allocating full/scalar will save some memory.
      full_beam_values_.resize(n_threads);
      scalar_beam_values_.resize(n_threads);

      for (size_t i = 0; i != n_threads; ++i) {
        patch_model_visibilities_[i].resize(
            {n_baselines, n_channels, n_correlations});
        full_beam_values_[i].resize(n_stations * n_channels);
        scalar_beam_values_[i].resize(n_stations * n_channels);
      }
    }
  }

  xt::xtensor<std::complex<double>, 3>& GetModel(size_t threadIndex) {
    return model_visibilities_[threadIndex];
  }

  xt::xtensor<std::complex<double>, 3>& GetPatchModel(size_t threadIndex) {
    return patch_model_visibilities_[threadIndex];
  }

  std::vector<aocommon::MatrixComplexDouble2x2>& GetFullBeamValues(
      size_t threadIndex) {
    return full_beam_values_[threadIndex];
  }

  std::vector<everybeam::complex_t>& GetScalarBeamValues(size_t threadIndex) {
    return scalar_beam_values_[threadIndex];
  }

  std::vector<std::shared_ptr<everybeam::Station>>& GetStationList() {
    return station_list_;
  }

 private:
  std::vector<xt::xtensor<std::complex<double>, 3>> model_visibilities_;
  std::vector<xt::xtensor<std::complex<double>, 3>> patch_model_visibilities_;
  std::vector<std::vector<aocommon::MatrixComplexDouble2x2>> full_beam_values_;
  std::vector<std::vector<everybeam::complex_t>> scalar_beam_values_;
  std::vector<std::shared_ptr<everybeam::Station>> station_list_;
  std::shared_ptr<everybeam::telescope::Telescope> telescope_;
};

}  // namespace base
}  // namespace dp3

#endif
