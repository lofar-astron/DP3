// NullStokes.cc: DP3 step class for zeroing out Stokes parameters
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Matthijs van der Wild

#include "NullStokes.h"

#include "base/DPInfo.h"

#include "common/ParameterSet.h"

#include <iostream>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

void MakeStokesIZero(std::complex<float>* visibilities) {
  std::complex<float>& xx = visibilities[0];
  std::complex<float>& yy = visibilities[3];
  xx = (xx - yy) * 0.5f;
  yy = -xx;
}

void MakeStokesQZero(std::complex<float>* visibilities) {
  std::complex<float>& xx = visibilities[0];
  std::complex<float>& yy = visibilities[3];
  xx = (xx + yy) * 0.5f;
  yy = xx;
}

void MakeStokesUZero(std::complex<float>* visibilities) {
  std::complex<float>& xy = visibilities[1];
  std::complex<float>& yx = visibilities[2];
  xy = (xy - yx) * 0.5f;
  yx = -xy;
}

void MakeStokesVZero(std::complex<float>* visibilities) {
  std::complex<float>& xy = visibilities[1];
  std::complex<float>& yx = visibilities[2];
  xy = (xy + yx) * 0.5f;
  yx = xy;
}

NullStokes::NullStokes(const common::ParameterSet& parset,
                       const std::string& prefix)
    : name_(prefix),
      modify_i_(parset.getBool(prefix + "modify_i", false)),
      modify_q_(parset.getBool(prefix + "modify_q", false)),
      modify_u_(parset.getBool(prefix + "modify_u", false)),
      modify_v_(parset.getBool(prefix + "modify_v", false)) {}

void NullStokes::show(std::ostream& os) const {
  os << "NullStokes " << name_ << '\n'
     << std::boolalpha << "modify_i " << modify_i_ << '\n'
     << "modify_q " << modify_q_ << '\n'
     << "modify_u " << modify_u_ << '\n'
     << "modify_v " << modify_v_ << '\n';
}

void NullStokes::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " NullStokes " << name_ << '\n';
}

bool NullStokes::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();
  std::complex<float>* visibilities = buffer->GetData().data();
  // The Stokes parameters are defined in terms of correlation of the electric
  // field. The visibilities corresponding to these correlations are denoted by
  // xx, xy, yx, yy. For example, Q = xx - yy, U = xy + yx. Here, the
  // visibilities are modified such that Q and/or U is zero, but I and V
  // unchanged. See e.g. eq. (4.28) and §4.5–4.7.2 of
  // https://link.springer.com/book/10.1007/978-3-319-44431-4
  for (std::size_t i = 0; i < buffer->GetData().size(); i += 4) {
    if (modify_i_) {
      MakeStokesIZero(&visibilities[i]);
    }
    if (modify_q_) {
      MakeStokesQZero(&visibilities[i]);
    }
    if (modify_u_) {
      MakeStokesUZero(&visibilities[i]);
    }
    if (modify_v_) {
      MakeStokesVZero(&visibilities[i]);
    }
  }
  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return true;
}

void NullStokes::updateInfo(const base::DPInfo& info_in) {
  Step::updateInfo(info_in);
}

void NullStokes::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

}  // namespace steps
}  // namespace dp3
