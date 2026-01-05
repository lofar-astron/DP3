// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Interpolate.h"

#include "base/DPBuffer.h"
#include "base/DPInfo.h"
#include "base/DP3.h"

#include "common/buffered_lane.h"
#include "common/ParameterSet.h"
#include "common/StringTools.h"

#include <iostream>
#include <iomanip>
#include <thread>

#include <aocommon/threadpool.h>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

Interpolate::Interpolate(const common::ParameterSet& parset,
                         const std::string& prefix)
    : name_(prefix),
      interpolated_pos_(0),
      window_size_(parset.getUint(prefix + "windowsize", 15)) {
  if (window_size_ % 2 != 1)
    throw std::runtime_error(
        "Window size of Interpolate action should be an odd number");

  kernel_lookup_.reserve(window_size_ * window_size_);
  for (int t = 0; t != int(window_size_); ++t) {
    int y = t - int(window_size_ / 2);
    for (int ch = 0; ch != int(window_size_); ++ch) {
      int x = ch - int(window_size_ / 2);
      double window_dist = double(x * x + y * y);
      // Gaussian function with sigma = 1
      // (evaluated with double prec, then converted to floats)
      double w = std::exp(window_dist * -0.5);
      kernel_lookup_.emplace_back(w);
    }
  }
}

void Interpolate::show(std::ostream& os) const {
  os << "Interpolate " << name_ << '\n';
  os << "  windowsize:     " << window_size_ << '\n';
}

void Interpolate::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " Interpolate " << name_ << '\n';
}

bool Interpolate::process(std::unique_ptr<base::DPBuffer> input_buffer) {
  timer_.start();
  // Collect the data in buffers.
  buffers_.emplace_back(std::move(input_buffer));
  // If we have a full window of data, interpolate everything
  // up to the middle of the window
  if (buffers_.size() >= window_size_) {
    size_t mid = window_size_ / 2;
    while (interpolated_pos_ <= mid) {
      interpolateTimestep(interpolated_pos_);
      ++interpolated_pos_;
    }
    // Buffers are only pushed to the next step when they are completely
    // out of the window. This is because flags need to be set to false,
    // however the flag information of the entire window is needed during
    // interpolation, so these can only be set to false after processing.
    sendFrontBufferToNextStep();
  }
  timer_.stop();
  return true;
}

void Interpolate::sendFrontBufferToNextStep() {
  auto front_buffer = std::move(buffers_.front());
  const auto& data_shape = front_buffer->GetData().shape();
  size_t num_pol = data_shape[2];
  size_t num_channels = data_shape[1];
  size_t num_baselines = data_shape[0];
  size_t n = num_pol * num_channels * num_baselines;
  // Set all flags to false
  bool* flags = front_buffer->GetFlags().data();
  std::complex<float>* data = front_buffer->GetData().data();
  std::fill(flags, flags + n, false);
  // Flag NaN values (values for which the entire window was flagged on input)
  for (size_t i = 0; i != n; ++i) {
    if (!std::isfinite(data[i].real()) || !std::isfinite(data[i].imag())) {
      // The datum value is also set to 0, because NaNs sometimes give problems
      // in certain software, even when they are flagged (e.g. in Sagecal).
      data[i] = 0.0;
      flags[i] = true;
    }
  }

  timer_.stop();
  getNextStep()->process(std::move(front_buffer));
  timer_.start();

  buffers_.pop_front();
  --interpolated_pos_;
}

void Interpolate::finish() {
  timer_.start();

  // Interpolate everything up to the end of the window
  while (interpolated_pos_ < buffers_.size()) {
    interpolateTimestep(interpolated_pos_);
    ++interpolated_pos_;
  }
  while (!buffers_.empty()) {
    sendFrontBufferToNextStep();
  }

  timer_.stop();

  getNextStep()->finish();
}

#define BUFFER_SIZE 1024

void Interpolate::interpolateTimestep(size_t index) {
  const auto& data_shape = buffers_.front()->GetData().shape();
  const size_t num_pol = data_shape[2];
  const size_t num_channels = data_shape[1];
  const size_t num_per_baseline = num_pol * num_channels;
  const size_t num_baselines = data_shape[0];

  lane_.resize(aocommon::ThreadPool::GetInstance().NThreads() * BUFFER_SIZE);
  common::lane_write_buffer<Sample> buflane(&lane_, BUFFER_SIZE);
  aocommon::ThreadPool::GetInstance().StartParallelExecution(
      [&](size_t) { interpolationThread(); });

  for (size_t bl = 0; bl < num_baselines; ++bl) {
    bool* flags = buffers_[index]->GetFlags().data() + bl * num_per_baseline;
    for (size_t ch = 0; ch != num_channels; ++ch) {
      for (size_t p = 0; p != num_pol; ++p) {
        if (*flags) {
          buflane.emplace(index, bl, ch, p);
        }
        ++flags;
      }
    }
  }
  buflane.write_end();

  aocommon::ThreadPool::GetInstance().FinishParallelExecution();
}

void Interpolate::interpolationThread() {
  common::lane_read_buffer<Sample> lane_buffer(&lane_, BUFFER_SIZE);
  Sample sample;
  while (lane_buffer.read(sample)) {
    interpolateSample(sample.timestep, sample.baseline, sample.channel,
                      sample.pol);
  }
}

void Interpolate::interpolateSample(size_t timestep, size_t baseline,
                                    size_t channel, size_t pol) {
  const auto& data_shape = buffers_.front()->GetData().shape();
  const size_t num_pol = data_shape[2];
  const size_t num_channels = data_shape[1];
  const size_t timestep_begin =
      (timestep > window_size_ / 2) ? (timestep - window_size_ / 2) : 0;
  const size_t timestep_end =
      std::min(timestep + window_size_ / 2 + 1, buffers_.size());
  const size_t channel_begin =
      (channel > window_size_ / 2) ? (channel - window_size_ / 2) : 0;
  const size_t channel_end =
      std::min(channel + window_size_ / 2 + 1, num_channels);

  std::complex<float> value_sum = 0.0;
  float window_sum = 0.0;

  for (size_t t = timestep_begin; t != timestep_end; ++t) {
    std::complex<float>* data =
        buffers_[t]->GetData().data() +
        (baseline * num_channels + channel_begin) * num_pol + pol;
    const bool* flags = buffers_[t]->GetFlags().data() +
                        (baseline * num_channels + channel_begin) * num_pol +
                        pol;
    const float* row =
        &kernel_lookup_[window_size_ * (t + int(window_size_ / 2) - timestep)];
    for (size_t ch = channel_begin; ch != channel_end; ++ch) {
      if (!*flags) {
        int x = ch + int(window_size_ / 2) - channel;
        float w = row[x];
        value_sum += *data * w;
        window_sum += w;
      }

      data += num_pol;
      flags += num_pol;
    }
  }
  // This write is multithreaded, but is allowed because this value is never
  // read from in the loops above (because flagged values are skipped).
  std::complex<float>& value =
      buffers_[timestep]
          ->GetData()
          .data()[(baseline * num_channels + channel) * num_pol + pol];
  if (window_sum != 0.0)
    value = value_sum / window_sum;
  else
    value = std::complex<float>(std::numeric_limits<float>::quiet_NaN(),
                                std::numeric_limits<float>::quiet_NaN());
}

}  // namespace steps
}  // namespace dp3
