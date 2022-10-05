// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to set the beam keywords in a ms

#ifndef DP3_STEPS_SETBEAM_H_
#define DP3_STEPS_SETBEAM_H_

#include <casacore/measures/Measures/MDirection.h>

#include <EveryBeam/correctionmode.h>

#include "InputStep.h"
#include "../base/DPBuffer.h"
#include "../base/Direction.h"
#include "../common/ParameterSet.h"

namespace dp3 {
namespace steps {

/// @brief DPPP step class to set the beam keywords in a ms
class SetBeam final : public Step {
 public:
  /// Parameters are obtained from the parset using the given prefix.
  SetBeam(InputStep* input, const common::ParameterSet& parameters,
          const string& prefix);

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return {}; }

  bool process(const base::DPBuffer& buffer) override;

  void finish() override{};

  void updateInfo(const base::DPInfo& info) override;

  void show(std::ostream&) const override;

 private:
  string _name;
  std::vector<string> _directionStr;
  casacore::MDirection _direction;
  everybeam::CorrectionMode _mode;
};

}  // namespace steps
}  // namespace dp3

#endif
