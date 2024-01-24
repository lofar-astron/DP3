// Split.cc: DPPP step class to Split visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Split.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <dp3/base/DP3.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

Split::Split(const common::ParameterSet& parset, const string& prefix)
    : name_(prefix),
      replace_parameters_(parset.getStringVector(prefix + "replaceparms")),
      sub_steps_() {
  // For each of the parameters, the values for each substep
  std::vector<std::vector<string>> replace_values;
  replace_values.reserve(replace_parameters_.size());

  for (const std::string& replace_parameter : replace_parameters_) {
    replace_values.emplace_back(parset.getStringVector(replace_parameter));
    if (replace_values.back().size() != replace_values.front().size()) {
      throw std::runtime_error(
          "Each parameter in replaceparms should have the same number of "
          "items (expected " +
          std::to_string(replace_values.front().size()) + ", got " +
          std::to_string(replace_values.back().size()) + " for step " +
          replace_parameter + ")");
    }
  }

  // num_splits, the number of 'new streams' that the data are split into is
  // determined from the replaced parameters.
  const size_t num_splits =
      replace_values.empty() ? 0 : replace_values.front().size();

  // Make a shallow copy to work around constness of parset
  common::ParameterSet parsetCopy(parset);

  // Create the substeps
  const size_t num_parameters = replace_parameters_.size();
  for (size_t i = 0; i != num_splits; ++i) {
    for (size_t j = 0; j != num_parameters; ++j) {
      parsetCopy.replace(replace_parameters_[j], replace_values[j][i]);
    }
    std::shared_ptr<Step> first_step = base::MakeStepsFromParset(
        parsetCopy, prefix, "steps", "", true, steps::Step::MsType::kRegular);
    if (first_step) {
      first_step->setPrevStep(this);
      sub_steps_.push_back(std::move(first_step));
    }
  }
}

common::Fields Split::getRequiredFields() const {
  common::Fields fields;
  for (const std::shared_ptr<Step>& first_step : sub_steps_) {
    fields |= base::GetChainRequiredFields(first_step);
  }
  return fields;
}

void Split::SetFieldsToWrite(const common::Fields& fields) {
  // Forward the provided fields to each sub-pipeline.
  for (std::shared_ptr<Step>& first_step : sub_steps_) {
    base::SetChainProvidedFields(first_step, fields);
  }
}

void Split::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  for (std::shared_ptr<Step>& step : sub_steps_) {
    step->setInfo(infoIn);
  }
}

void Split::show(std::ostream& os) const {
  os << "Split " << name_ << '\n'
     << "  replace parameters:" << replace_parameters_ << '\n';
  // Show the steps.
  for (unsigned int i = 0; i < sub_steps_.size(); ++i) {
    os << "Split substep " << (i + 1) << " of " << sub_steps_.size() << '\n';
    std::shared_ptr<Step> step = sub_steps_[i];
    while (step) {
      step->show(os);
      step = step->getNextStep();
    }
  }
}

void Split::showTimings(std::ostream& os, double duration) const {
  for (unsigned int i = 0; i < sub_steps_.size(); ++i) {
    std::shared_ptr<Step> step = sub_steps_[i];
    while (step) {
      step->showTimings(os, duration);
      step = step->getNextStep();
    }
  }
}

bool Split::process(std::unique_ptr<DPBuffer> buffer) {
  for (std::shared_ptr<Step>& step : sub_steps_) {
    step->process(std::make_unique<DPBuffer>(*buffer));
  }
  return false;
}

void Split::finish() {
  // Let the next steps finish.
  for (std::shared_ptr<Step>& step : sub_steps_) {
    step->finish();
  }
}

void Split::setNextStep(std::shared_ptr<Step>) {
  throw std::runtime_error("Split must be the last step in a step chain");
}

}  // namespace steps
}  // namespace dp3
