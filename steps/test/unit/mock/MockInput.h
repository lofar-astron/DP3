// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_TEST_UNIT_MOCKINPUT_H_
#define DP3_STEPS_TEST_UNIT_MOCKINPUT_H_

#include <boost/test/unit_test.hpp>

#include "steps/InputStep.h"

namespace dp3 {
namespace steps {

class MockInput : public InputStep {
 public:
  MockInput();
  ~MockInput() override;

  common::Fields getRequiredFields() const override;
  common::Fields getProvidedFields() const override;

  void finish() override;
  void show(std::ostream&) const override;
};
}  // namespace steps
}  // namespace dp3

#endif
