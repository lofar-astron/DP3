// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Predict.h"

#include <ostream>
#include <string>

#include "Averager.h"
#include "BDAAverager.h"
#include "BDAExpander.h"
#include "OnePredict.h"
#include "Upsample.h"

#include <dp3/base/BDABuffer.h>

#include "../common/ParameterSet.h"

using dp3::base::BDABuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Predict::Predict(const common::ParameterSet& parset, const string& prefix,
                 MsType input_type)
    : ms_type_(input_type),
      predict_step_(std::make_shared<OnePredict>(parset, prefix,
                                                 std::vector<std::string>())) {
  Initialize(parset, prefix, input_type);
}

Predict::Predict(const common::ParameterSet& parset, const string& prefix,
                 const std::vector<std::string>& source_patterns,
                 MsType input_type)
    : ms_type_(input_type),
      predict_step_(
          std::make_shared<OnePredict>(parset, prefix, source_patterns)) {
  Initialize(parset, prefix, input_type);
}

void Predict::Initialize(const common::ParameterSet& parset,
                         const string& prefix, MsType input_type) {
  const unsigned int time_smearing_factor =
      parset.getUint(prefix + "correcttimesmearing", 1);

  // Create the steps that this step manages, in order of how they need to be
  // connected. These are called 'internal' steps here to differentiate them
  // from the other steps outside Predict, but they are not substeps.
  if (input_type == MsType::kBda) {
    internal_steps_.push_back(std::make_shared<BDAExpander>(prefix));
  }
  if (time_smearing_factor > 1) {
    internal_steps_.push_back(std::make_shared<Upsample>(
        prefix + "upsample", time_smearing_factor, true));
  }

  internal_steps_.push_back(predict_step_);

  if (time_smearing_factor > 1) {
    internal_steps_.push_back(std::make_shared<Averager>(prefix + "averager", 1,
                                                         time_smearing_factor));
  }
  if (input_type == MsType::kBda) {
    bda_averager_ = std::make_shared<BDAAverager>(parset, prefix, false);
    internal_steps_.push_back(bda_averager_);
  }

  // Connect the steps
  Step::setNextStep(internal_steps_.front());
  for (size_t i = 1; i < internal_steps_.size(); ++i) {
    internal_steps_[i - 1]->setNextStep(internal_steps_[i]);
  }
  // The last 'internal' step is connected in Predict::setNextStep()
}

void Predict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  if (bda_averager_) {
    bda_averager_->set_averaging_params(
        info().ntimeAvgs(), info().BdaChanFreqs(), info().BdaChanWidths());
  }
}

base::Direction Predict::GetFirstDirection() const {
  return predict_step_->GetFirstDirection();
}

void Predict::setNextStep(std::shared_ptr<Step> next_step) {
  internal_steps_.back()->setNextStep(next_step);
}

void Predict::SetOperation(const std::string& operation) {
  predict_step_->SetOperation(operation);
}

void Predict::SetThreadData(std::mutex* mutex) {
  predict_step_->SetThreadData(mutex);
}

void Predict::show(std::ostream& os) const { os << "Predict" << '\n'; }

bool Predict::process(std::unique_ptr<base::DPBuffer> buffer) {
  return getNextStep()->process(std::move(buffer));
}

bool Predict::process(std::unique_ptr<BDABuffer> bda_buffer) {
  bda_averager_->set_next_desired_buffersize(bda_buffer->GetNumberOfElements());
  return getNextStep()->process(std::move(bda_buffer));
}

void Predict::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
