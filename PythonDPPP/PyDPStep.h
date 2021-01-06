// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_PYDPSTEP_H
#define DPPP_PYDPSTEP_H

#include "../DPPP/DPStep.h"
#include "../Common/ParameterSet.h"

namespace DP3 {
namespace DPPP {

class PyDPStep : public DPStep {
 public:
  static DPStep::ShPtr create_instance(DPInput* input,
                                       const ParameterSet& parset,
                                       const string& prefix);

 private:
  using DPStep::DPStep;
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DPPP_PYDPSTEP_H
