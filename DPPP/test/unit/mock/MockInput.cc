// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockInput.h"

namespace DP3 {
namespace DPPP {

MockInput::MockInput() {}
MockInput::~MockInput() {}

void MockInput::getUVW(const casacore::RefRows&, double, DPBuffer& buffer) {
  BOOST_TEST(!buffer.getUVW().empty());
}
void MockInput::getWeights(const casacore::RefRows&, DPBuffer& buffer) {
  BOOST_TEST(!buffer.getWeights().empty());
}
void MockInput::finish() { BOOST_ERROR("Unexpected finish() call"); }
void MockInput::show(std::ostream&) const {
  BOOST_ERROR("Unexpected show() call");
}

}  // namespace DPPP
}  // namespace DP3