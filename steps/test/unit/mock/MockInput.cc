// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MockInput.h"

namespace dp3 {
namespace steps {

MockInput::MockInput() {}
MockInput::~MockInput() {}

common::Fields MockInput::getRequiredFields() const {
  throw std::runtime_error("Unexpected getRequiredFields call");
}

common::Fields MockInput::getProvidedFields() const {
  throw std::runtime_error("Unexpected getProvidedFields call");
}

void MockInput::finish() { BOOST_ERROR("Unexpected finish() call"); }
void MockInput::show(std::ostream&) const {
  BOOST_ERROR("Unexpected show() call");
}

}  // namespace steps
}  // namespace dp3
