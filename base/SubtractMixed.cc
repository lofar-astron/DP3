// SubtractMixed.cc: Subtract visibilities from a buffer after weighting by
// mixing coefficients.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "SubtractMixed.h"

namespace dp3 {
namespace base {

float subtract(size_t nBaseline, size_t nChannel,
               const_cursor<Baseline> baselines,
               cursor<std::complex<float>> data,
               const_cursor<std::complex<double>> model,
               const_cursor<std::complex<double>> weight) {
  float var_before = 0.0, var_after = 0.0;
  for (size_t bl = 0; bl < nBaseline; ++bl) {
    const size_t p = baselines->first;
    const size_t q = baselines->second;

    if (p != q) {
#pragma GCC ivdep
      for (size_t ch = 0; ch < nChannel; ++ch) {
        // Subtract weighted model from data.
        float data_abs = std::abs(*data);
        var_before += (!std::isnan(data_abs) ? data_abs * data_abs : 0.0f);
        *data -= std::complex<float>((*weight) * (*model));
        float res_abs = std::abs(*data);
        var_after += (!std::isnan(res_abs) ? res_abs * res_abs : 0.0f);
        ++weight;
        ++model;
        ++data;
        *data -= std::complex<float>((*weight) * (*model));
        ++weight;
        ++model;
        ++data;
        *data -= std::complex<float>((*weight) * (*model));
        ++weight;
        ++model;
        ++data;
        data_abs = std::abs(*data);
        var_before += (!std::isnan(data_abs) ? data_abs * data_abs : 0.0f);
        *data -= std::complex<float>((*weight) * (*model));
        res_abs = std::abs(*data);
        var_after += (!std::isnan(res_abs) ? res_abs * res_abs : 0.0f);
        ++weight;
        ++model;
        ++data;

        // Move to the next channel.
        weight -= 4;
        weight.forward(1);
        model -= 4;
        model.forward(1);
        data -= 4;
        data.forward(1);
      }  // Channels.

      weight.backward(1, nChannel);
      model.backward(1, nChannel);
      data.backward(1, nChannel);
    }

    // Move to the next baseline.
    weight.forward(2);
    model.forward(2);
    data.forward(2);
    ++baselines;
  }  // Baselines.
  return (var_before / (var_after + 1e-6f));
}

}  // namespace base
}  // namespace dp3
