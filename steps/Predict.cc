// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Predict.h"

#include "Averager.h"
#include "OnePredict.h"
#include "Upsample.h"

#include "../common/ParameterSet.h"

#include <ostream>
#include <string>

namespace dp3 {
namespace steps {

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix)
    : upsample_step_(),
      predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix)),
      averager_step_() {
  Initialize(input_step, parset, prefix);
}

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix,
                 const std::vector<std::string>& source_patterns)
    : upsample_step_(),
      predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix,
                                                 source_patterns)),
      averager_step_() {
  Initialize(input_step, parset, prefix);
}

void Predict::Initialize(InputStep& input_step,
                         const common::ParameterSet& parset,
                         const string& prefix) {
  const unsigned int time_smearing_factor =
      parset.getUint(prefix + "correcttimesmearing", 1);
  if (time_smearing_factor > 1) {
    upsample_step_ = std::make_shared<Upsample>(prefix + "upsample",
                                                time_smearing_factor, true);
    averager_step_ = std::make_shared<Averager>(input_step, prefix + "averager",
                                                1, time_smearing_factor);
    Step::setNextStep(upsample_step_);
    upsample_step_->setNextStep(predict_step_);
    predict_step_->setNextStep(averager_step_);
  } else {
    // Without time smearing, upsampling and averaging is not needed.
    Step::setNextStep(predict_step_);
  }
}

std::pair<double, double> Predict::GetFirstDirection() const {
  return predict_step_->GetFirstDirection();
}

void Predict::setNextStep(std::shared_ptr<Step> next_step) {
  if (averager_step_)
    averager_step_->setNextStep(next_step);
  else
    predict_step_->setNextStep(next_step);
}

void Predict::SetOperation(const std::string& operation) {
  predict_step_->SetOperation(operation);
}

void Predict::SetThreadData(aocommon::ThreadPool& pool, std::mutex* mutex) {
  predict_step_->SetThreadData(pool, mutex);
}

void Predict::SetPredictBuffer(
    std::shared_ptr<base::PredictBuffer> predict_buffer) {
  predict_step_->SetPredictBuffer(predict_buffer);
}

void Predict::show(std::ostream& os) const { os << "Predict" << '\n'; }

bool Predict::process(const base::DPBuffer& buffer) {
  return getNextStep()->process(buffer);
}

void Predict::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
