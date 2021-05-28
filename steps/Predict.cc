// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Predict.h"

#include "BdaPredict.h"
#include "OnePredict.h"

#include <ostream>

namespace dp3 {
namespace steps {

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix)
    : predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix)) {
  Step::setNextStep(predict_step_);
}

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix,
                 const std::vector<std::string>& source_patterns)
    : predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix,
                                                 source_patterns)) {
  Step::setNextStep(predict_step_);
}

std::pair<double, double> Predict::GetFirstDirection() const {
  return predict_step_->GetFirstDirection();
}

void Predict::setNextStep(std::shared_ptr<Step> next_step) {
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
  return predict_step_->process(buffer);
}

void Predict::finish() { predict_step_->finish(); }

}  // namespace steps
}  // namespace dp3
