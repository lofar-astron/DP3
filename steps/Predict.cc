// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Predict.h"

#include <ostream>
#include <string>

#include "Averager.h"
#include "BDAAverager.h"
#include "BdaExpander.h"
#include "OnePredict.h"
#include "Upsample.h"

#include <dp3/base/BdaBuffer.h>

#include "../common/ParameterSet.h"

using dp3::base::BdaBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

/**
 * Internal Step class which clears the meta data changed flag if that flag is
 * also clear in the Predict step.
 * If Predict creates a chain consisting of Upsample, Averager, BdaExpander,
 * and/or BdaAverager steps, the meta data changed flag will be set at the end
 * of that chain, since those steps indeed change the data structure.
 * However, since the steps in the chain cancel each other, Predict should
 * restore the flag to its value at the beginning of the chain.
 */
class RestoreMetaDataChangedStep : public Step {
 public:
  RestoreMetaDataChangedStep(Step* reference_step)
      : reference_step_(reference_step) {}

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return {}; }

  void updateInfo(const DPInfo& info_in) override {
    Step::updateInfo(info_in);
    if (!reference_step_->getInfoIn().metaChanged())
      GetWritableInfoOut().clearMetaChanged();
  }

  bool accepts(MsType dt) const override {
    return dt == reference_step_->outputs();
  }

  MsType outputs() const override { return reference_step_->outputs(); }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override {
    return getNextStep()->process(std::move(buffer));
  }

  bool process(std::unique_ptr<base::BdaBuffer> bda_buffer) override {
    return getNextStep()->process(std::move(bda_buffer));
  }

  void finish() override { getNextStep()->finish(); }

  void show(std::ostream& os) const override {}

 private:
  /// Non-owning pointer to the Predict step that contains this step, which
  /// is the first step in the step chain created by Predict.
  Step* reference_step_;
};

Predict::Predict(const common::ParameterSet& parset, const std::string& prefix,
                 MsType input_type)
    : ms_type_(input_type),
      predict_step_(std::make_shared<OnePredict>(parset, prefix,
                                                 std::vector<std::string>())) {
  Initialize(parset, prefix, input_type);
}

Predict::Predict(const common::ParameterSet& parset, const std::string& prefix,
                 const std::vector<std::string>& source_patterns,
                 MsType input_type)
    : ms_type_(input_type),
      predict_step_(
          std::make_shared<OnePredict>(parset, prefix, source_patterns)) {
  Initialize(parset, prefix, input_type);
}

void Predict::Initialize(const common::ParameterSet& parset,
                         const std::string& prefix, MsType input_type) {
  const unsigned int time_smearing_factor =
      parset.getUint(prefix + "correcttimesmearing", 1);

  // Create the steps that this step manages, in order of how they need to be
  // connected. These are called 'internal' steps here to differentiate them
  // from the other steps outside Predict, but they are not substeps.
  if (input_type == MsType::kBda) {
    internal_steps_.push_back(std::make_shared<BdaExpander>(prefix));
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
    bda_averager_ = std::make_shared<BdaAverager>(parset, prefix, false);
    internal_steps_.push_back(bda_averager_);
  }

  internal_steps_.push_back(std::make_shared<RestoreMetaDataChangedStep>(this));

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
    bda_averager_->set_averaging_params(getInfoOut().ntimeAvgs(),
                                        getInfoOut().BdaChanFreqs(),
                                        getInfoOut().BdaChanWidths());
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

bool Predict::process(std::unique_ptr<BdaBuffer> bda_buffer) {
  bda_averager_->set_next_desired_buffersize(bda_buffer->GetNumberOfElements());
  return getNextStep()->process(std::move(bda_buffer));
}

void Predict::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
