// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockInput.h"

namespace dp3 {
namespace steps {

MockInput::MockInput() {}
MockInput::~MockInput() {}

void MockInput::getUVW(const casacore::RefRows&, double,
                       base::DPBuffer& buffer) {
  BOOST_TEST(!buffer.getUVW().empty());
}
void MockInput::getWeights(const casacore::RefRows&, base::DPBuffer& buffer) {
  BOOST_TEST(!buffer.getWeights().empty());
}
bool MockInput::getFullResFlags(const casacore::RefRows&,
                                base::DPBuffer& buffer) {
  BOOST_TEST(!buffer.getFullResFlags().empty());
  return true;
}
void MockInput::finish() { BOOST_ERROR("Unexpected finish() call"); }
void MockInput::show(std::ostream&) const {
  BOOST_ERROR("Unexpected show() call");
}

}  // namespace steps
}  // namespace dp3
