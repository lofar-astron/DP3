// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Clipper.h"
#include "NullStep.h"
#include "MsColumnReader.h"

#include <iostream>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xview.hpp>

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
      frequency_step_(parset.getInt(prefix + "freqstep", 1)),
      flag_all_correlations_(
          parset.getBool(prefix + "flagallcorrelations", true)),
      max_amplitude_(parset.getFloat(prefix + "amplmax", 0.0)) {
  SetPredict(
      std::make_shared<OnePredict>(parset, prefix, std::vector<std::string>()));
}

void Clipper::SetPredict(std::shared_ptr<Step> substep) {
  predict_step_ = substep;
  result_step_ = std::make_shared<ResultStep>();
  predict_step_->setNextStep(result_step_);
}

void Clipper::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);
  // Create a metadata copy for the internal Predict step.
  DPInfo predict_info = info_in;

  if (max_amplitude_ == 0.0) {
    std::string antennaSet(info_in.antennaSet());
    if (antennaSet.substr(0, 3) == "LBA") {
      max_amplitude_ = 50.0;
    } else {
      max_amplitude_ = 5.0;
    }
  }

  // Adapt the channel frequencies and widths for the prediction
  const std::vector<double>& input_frequencies = infoOut().chanFreqs();
  const std::vector<double>& input_widths = infoOut().chanWidths();
  const size_t n_channels_in = infoOut().nchan();
  const size_t n_channels_out =
      (n_channels_in + frequency_step_ - 1) / frequency_step_;
  std::vector<double> predict_frequencies(n_channels_out);
  std::vector<double> predict_widths(n_channels_out);
  for (size_t channel = 0; channel < n_channels_out; ++channel) {
    predict_frequencies[channel] = input_frequencies[channel * frequency_step_];
    predict_widths[channel] = input_widths[channel * frequency_step_];
  }
  predict_info.setChannels(std::move(predict_frequencies),
                           std::move(predict_widths));
  predict_step_->setInfo(predict_info);
}

void Clipper::show(std::ostream& os) const {
  os << "Clipper " << name_ << '\n';
  os << "  time step:             " << time_step_ << '\n';
  os << "  frequency step:        " << frequency_step_ << '\n';
  os << "  maximum amplitude:     " << max_amplitude_ << '\n';
  os << "  flag all correlations: " << std::boolalpha << flag_all_correlations_
     << '\n';

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
    // Copy the buffer, and reduce the number of channels in the buffer
    common::Fields predict_fields = predict_step_->getRequiredFields();
    std::unique_ptr<DPBuffer> substep_buffer =
        std::make_unique<DPBuffer>(*buffer, predict_fields);
    predict_step_->process(std::move(substep_buffer));
    std::unique_ptr<DPBuffer> result_buffer = result_step_->take();

    // After the prediction, flag the predicted values that exceed
    // max_amplitude_. Copy each resulting flag to the buffer
    // frequency_step_ times.
    last_flags_.resize(buffer->GetFlags().shape());
    const size_t n_channels_in = info().nchan();
    const size_t n_channels_out =
        (n_channels_in + frequency_step_ - 1) / frequency_step_;
    for (size_t channel = 0; channel < n_channels_out; ++channel) {
      const size_t start = channel * frequency_step_;
      const size_t steps = std::min(frequency_step_, n_channels_in - start);
      const size_t stop = start + steps;
      xt::xtensor<bool, 2> result_flags =
          xt::view(xt::abs(result_buffer->GetData()) > max_amplitude_,
                   xt::all(), channel, xt::all());
      // flag all correlations corresponding to a single baseline if
      // a single correlation for this baseline exceeds amplmax_
      if (flag_all_correlations_) {
        xt::xarray<int> baseline_flag_count = xt::sum(result_flags, 1);
        auto selected_baselines =
            xt::flatten_indices(xt::argwhere(baseline_flag_count));
        auto correlations_to_flag =
            xt::view(result_flags, xt::keep(selected_baselines), xt::all());
        correlations_to_flag = true;
      }
      for (size_t i = start; i < stop; i++) {
        xt::view(last_flags_, xt::all(), i, xt::all()) = result_flags;
      }
    }
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
