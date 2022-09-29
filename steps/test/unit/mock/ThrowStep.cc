// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ThrowStep.h"

#include <boost/test/unit_test.hpp>

namespace dp3 {
namespace steps {
namespace test {

ThrowStep::ThrowStep() {}
ThrowStep::~ThrowStep() {}

common::Fields ThrowStep::getRequiredFields() const {
  throw std::runtime_error("Unexpected getRequiredFields call");
}

common::Fields ThrowStep::getProvidedFields() const {
  throw std::runtime_error("Unexpected getProvidedFields call");
}

void ThrowStep::updateInfo(const base::DPInfo&) {
  BOOST_ERROR("Unexpected updateInfo() call");
}

bool ThrowStep::process(const base::DPBuffer&) {
  throw std::runtime_error("Unexpected process call for regular data");
}

bool ThrowStep::process(std::unique_ptr<base::BDABuffer>) {
  throw std::runtime_error("Unexpected process call for BDA data");
}

void ThrowStep::finish() { BOOST_ERROR("Unexpected finish() call"); }

void ThrowStep::show(std::ostream&) const {
  BOOST_ERROR("Unexpected show() call");
}

}  // namespace test
}  // namespace steps
}  // namespace dp3
