// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PYDPSTEP_H
#define DPPP_PYDPSTEP_H

#include "../steps/Step.h"
#include "../common/ParameterSet.h"

namespace dp3 {
namespace pythondp3 {

class PyStep : public steps::Step {
 public:
  static Step::ShPtr create_instance(steps::InputStep* input,
                                     const common::ParameterSet& parset,
                                     const string& prefix);

 private:
  using Step::Step;
};

}  // namespace pythondp3
}  // namespace dp3

#endif  // DPPP_PYDPSTEP_H
