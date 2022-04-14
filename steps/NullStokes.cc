// NullStokes.cc: DP3 step class for zeroing out Stokes parameters
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Matthijs van der Wild

#include "NullStokes.h"

#include "../base/DPInfo.h"

#include "../common/ParameterSet.h"

#include <iostream>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

NullStokes::NullStokes(InputStep& input, const common::ParameterSet& parset,
                       const std::string& prefix)
    : name_(prefix),
      modify_q_(parset.getBool(prefix + "modify_q", false)),
      modify_u_(parset.getBool(prefix + "modify_u", false)) {}

void NullStokes::show(std::ostream& os) const {
  os << "NullStokes " << name_ << '\n'
     << std::boolalpha << "modify_q " << modify_q_ << '\n'
     << "modify_u " << modify_u_ << '\n';
}

void NullStokes::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " NullStokes " << name_ << '\n';
}

bool NullStokes::process(const base::DPBuffer& buffer) {
  timer_.start();
  const casacore::IPosition& shape = buffer.getData().shape();
  const size_t n_channels = shape[1];
  const size_t n_baselines = shape[2];
  casacore::Array<casacore::Complex> data(shape);
  std::complex<float>* output_visibilities = data.data();
  const std::complex<float>* input_visibilities = buffer.getData().data();
  // The Stokes parameters are defined in terms of correlation of the electric
  // field. The visibilities corresponding to these correlations are denoted by
  // xx, xy, yx, yy. For example, Q = xx - yy, U = xy + yx. Here, the
  // visibilities are modified such that Q and/or U is zero, but I and V
  // unchanged. See e.g. eq. (4.28) and §4.5–4.7.2 of
  // https://link.springer.com/book/10.1007/978-3-319-44431-4
  for (size_t i = 0; i < n_baselines; ++i) {
    for (size_t j = 0; j < n_channels; ++j) {
      std::complex<float> xx = input_visibilities[0];
      std::complex<float> xy = input_visibilities[1];
      std::complex<float> yx = input_visibilities[2];
      std::complex<float> yy = input_visibilities[3];
      // Zero out appropriate Stokes parameter(s)
      if (modify_q_) {
        xx = (xx + yy) / 2.0f;
        yy = xx;
      }
      if (modify_u_) {
        xy = (xy - yx) / 2.0f;
        yx = -xy;
      }
      output_visibilities[0] = xx;
      output_visibilities[1] = xy;
      output_visibilities[2] = yx;
      output_visibilities[3] = yy;
      input_visibilities += 4;
      output_visibilities += 4;
    }
  }
  DPBuffer processedData(buffer);
  processedData.setData(data);
  timer_.stop();
  getNextStep()->process(processedData);
  return true;
}

void NullStokes::updateInfo(const base::DPInfo& info_in) {
  info() = info_in;
  if (modify_q_ || modify_u_) {
    info().setNeedVisData();
  }
}

void NullStokes::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

}  // namespace steps
}  // namespace dp3
