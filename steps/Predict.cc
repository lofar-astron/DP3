// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Predict.h"
#include "Averager.h"
#include "BDAAverager.h"
#include "BDAExpander.h"
#include "OnePredict.h"
#include "Upsample.h"

#include "../base/BDABuffer.h"

#include "../common/ParameterSet.h"

#include <ostream>
#include <string>

using dp3::base::BDABuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix, MsType input_type)
    : ms_type_(input_type),
      predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix,
                                                 std::vector<std::string>())) {
  Initialize(input_step, parset, prefix, input_type);
}

Predict::Predict(InputStep& input_step, const common::ParameterSet& parset,
                 const string& prefix,
                 const std::vector<std::string>& source_patterns,
                 MsType input_type)
    : ms_type_(input_type),
      predict_step_(std::make_shared<OnePredict>(&input_step, parset, prefix,
                                                 source_patterns)) {
  Initialize(input_step, parset, prefix, input_type);
}

void Predict::Initialize(InputStep& input_step,
                         const common::ParameterSet& parset,
                         const string& prefix, MsType input_type) {
  if (input_type == MsType::kBda) {
    steps_before_predict_.push_back(std::make_shared<BDAExpander>(prefix));
  }

  const unsigned int time_smearing_factor =
      parset.getUint(prefix + "correcttimesmearing", 1);
  if (time_smearing_factor > 1) {
    steps_before_predict_.push_back(std::make_shared<Upsample>(
        prefix + "upsample", time_smearing_factor, true));
    steps_after_predict_.push_back(std::make_shared<Averager>(
        input_step, prefix + "averager", 1, time_smearing_factor));
  }

  if (input_type == MsType::kBda) {
    bda_averager_ =
        std::make_shared<BDAAverager>(input_step, parset, prefix, false);
    steps_after_predict_.push_back(bda_averager_);
  }

  if (!steps_before_predict_.empty()) {
    Step::setNextStep(steps_before_predict_.front());

    steps_before_predict_.push_back(predict_step_);

    for (size_t i = 1; i < steps_before_predict_.size(); ++i) {
      steps_before_predict_[i - 1]->setNextStep(steps_before_predict_[i]);
    }
    steps_before_predict_.back()->setNextStep(steps_after_predict_.front());

    for (size_t i = 1; i < steps_after_predict_.size(); ++i) {
      steps_after_predict_[i - 1]->setNextStep(steps_after_predict_[i]);
    }

  } else {
    // Without time smearing or bda, extra steps are not needed.
    Step::setNextStep(predict_step_);
  }
}

void Predict::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  if (bda_averager_) {
    bda_averager_->set_averaging_params(
        info().ntimeAvgs(), info().BdaChanFreqs(), info().BdaChanWidths());
  }
  info().setNeedVisData();
  info().setWriteData();
}

base::Direction Predict::GetFirstDirection() const {
  return predict_step_->GetFirstDirection();
}

void Predict::setNextStep(std::shared_ptr<Step> next_step) {
  if (!steps_after_predict_.empty()) {
    steps_after_predict_.back()->setNextStep(next_step);
  } else {
    predict_step_->setNextStep(next_step);
  }
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

bool Predict::process(std::unique_ptr<BDABuffer> bda_buffer) {
  bda_averager_->set_next_desired_buffersize(bda_buffer->GetNumberOfElements());
  return getNextStep()->process(std::move(bda_buffer));
}

void Predict::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
