// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PREDICT_MODEL_H_
#define DP3_BASE_PREDICT_MODEL_H_

#include <complex>
#include <vector>

#include <aocommon/xt/utensor.h>

namespace dp3::base {

class PredictModel {
 public:
  explicit PredictModel(size_t n_threads, size_t n_correlations,
                        size_t n_channels, size_t n_baselines,
                        bool include_beam)
      : model_visibilities_(n_threads) {
    for (size_t i = 0; i != n_threads; ++i) {
      model_visibilities_[i].resize({n_baselines, n_channels, n_correlations});
    }

    if (include_beam) {
      patch_model_visibilities_.resize(n_threads);
      for (size_t i = 0; i != n_threads; ++i) {
        patch_model_visibilities_[i].resize(
            {n_baselines, n_channels, n_correlations});
      }
    }
  }

  aocommon::xt::UTensor<std::complex<double>, 3>& GetModel(
      size_t thread_index) {
    return model_visibilities_[thread_index];
  }

  aocommon::xt::UTensor<std::complex<double>, 3>& GetPatchModel(
      size_t thread_index) {
    return patch_model_visibilities_[thread_index];
  }

 private:
  std::vector<aocommon::xt::UTensor<std::complex<double>, 3>>
      model_visibilities_;
  std::vector<aocommon::xt::UTensor<std::complex<double>, 3>>
      patch_model_visibilities_;
};

}  // namespace dp3::base

#endif
