// SubtractNew.cc: Subtract visibilities from a buffer after weighting by
// mixing coefficients.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include "SubtractNew.h"

namespace dp3 {
namespace base {

void subtract(size_t nBaseline, size_t nChannel,
              const_cursor<Baseline> baselines, cursor<fcomplex> data,
              const_cursor<dcomplex> model, const_cursor<dcomplex> weight,
              std::vector<float>& ampl) {
  dcomplex vis[4];
  for (size_t bl = 0; bl < nBaseline; ++bl) {
    // Only for cross correlations.
    if (baselines->first != baselines->second) {
      for (size_t ch = 0; ch < nChannel; ++ch) {
        for (size_t cr = 0; cr < 4; ++cr) {
          vis[cr] = (*weight++) * (*model++);
          *data++ -= static_cast<fcomplex>(vis[cr]);
        }
        // Return subtracted amplitude for middle channel.
        if (ch == nChannel / 2) {
          ampl[bl] = (abs(vis[0]) + abs(vis[3])) * 0.5;
        }
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
    } else {
      ampl[bl] = 0.;
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
