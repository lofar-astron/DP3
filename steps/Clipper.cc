// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Clipper.h"
#include "NullStep.h"
#include "MsColumnReader.h"

#include <iostream>
#include <xtensor/xcomplex.hpp>

#include "../base/FlagCounter.h"
#include <dp3/base/DP3.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Clipper::Clipper(const common::ParameterSet& parset, const std::string& prefix)
    : name_(prefix),
      counter_(0),
      time_step_(parset.getInt(prefix + "timestep", 5)),
      ampl_max_(parset.getFloat(prefix + "amplmax", 0.0)) {
  predict_step_ =
      std::make_shared<OnePredict>(parset, prefix, std::vector<std::string>());
  result_step_ = std::make_shared<ResultStep>();
  predict_step_->setNextStep(result_step_);
}

void Clipper::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);
  if (ampl_max_ == 0.0) {
    std::string antennaSet(info_in.antennaSet());
    if (antennaSet.substr(0, 3) == "LBA") {
      ampl_max_ = 50.0;
    } else {
      ampl_max_ = 5.0;
    }
  }
  predict_step_->updateInfo(info_in);
}

void Clipper::show(std::ostream& os) const {
  os << "Clipper " << name_ << '\n';
  os << "  time step:     " << time_step_ << '\n';
  os << "  max amplitude: " << ampl_max_ << '\n';
  predict_step_->show(os);
}

void Clipper::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " Clipper " << name_ << '\n';
}

bool Clipper::process(std::unique_ptr<DPBuffer> buffer) {
  timer_.start();

  if (counter_ % time_step_ == 0) {
    std::unique_ptr<DPBuffer> substep_buffer =
        std::make_unique<DPBuffer>(*buffer);
    predict_step_->process(std::move(substep_buffer));
    std::unique_ptr<DPBuffer> result_buffer = result_step_->take();

    last_flags_ = xt::eval(xt::abs(result_buffer->GetData())) > ampl_max_;
    counter_ = 0;
  }
  counter_++;

  buffer->GetFlags() = buffer->GetFlags() || last_flags_;

  timer_.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

void Clipper::finish() { getNextStep()->finish(); }

}  // namespace steps
}  // namespace dp3
