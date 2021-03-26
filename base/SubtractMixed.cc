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

void subtract(size_t nBaseline, size_t nChannel,
              const_cursor<Baseline> baselines, cursor<fcomplex> data,
              const_cursor<dcomplex> model, const_cursor<dcomplex> weight) {
  for (size_t bl = 0; bl < nBaseline; ++bl) {
    const size_t p = baselines->first;
    const size_t q = baselines->second;

    if (p != q) {
      for (size_t ch = 0; ch < nChannel; ++ch) {
        // Subtract weighted model from data.
        *data -= static_cast<fcomplex>((*weight) * (*model));
        ++weight;
        ++model;
        ++data;
        *data -= static_cast<fcomplex>((*weight) * (*model));
        ++weight;
        ++model;
        ++data;
        *data -= static_cast<fcomplex>((*weight) * (*model));
        ++weight;
        ++model;
        ++data;
        *data -= static_cast<fcomplex>((*weight) * (*model));
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
}

}  // namespace base
}  // namespace dp3
