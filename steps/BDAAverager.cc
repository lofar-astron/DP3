// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "BDAAverager.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

#include <dp3/base/BdaBuffer.h>
#include "../common/Epsilon.h"
#include "../common/ParameterSet.h"

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

BdaAverager::BaselineBuffer::BaselineBuffer(std::size_t _time_factor,
                                            std::size_t n_input_channels,
                                            std::size_t n_output_channels,
                                            std::size_t n_correlations)
    : times_added(0),
      time_factor(_time_factor),
      input_channel_indices(),
      starttime(0.0),
      interval(0.0),
      exposure(0.0),
      data(),
      weights(n_output_channels * n_correlations, 0.0f),
      uvw{0.0f, 0.0f, 0.0f} {
  // Determine the start and end input channel for each output channel.
  // This list always ends with the number of output channels, which allows
  // loops over input_channel_indices[o] until input_channel_indices[o+1],
  // where o is an output channel index.
  assert(n_output_channels <= n_input_channels);
  input_channel_indices.reserve(n_output_channels + 1);
  for (std::size_t i = 0; i <= n_output_channels; ++i) {
    input_channel_indices.push_back((i * n_input_channels) / n_output_channels);
  }
  assert(input_channel_indices.size() == n_output_channels + 1);
  assert(n_input_channels == input_channel_indices.back());
}

void BdaAverager::BaselineBuffer::Clear() {
  times_added = 0;
  starttime = 0.0;
  interval = 0.0;
  exposure = 0.0;
  for (auto& name_vector : data) {
    std::fill(name_vector.second.begin(), name_vector.second.end(),
              std::complex<float>{0.0f, 0.0f});
  }
  std::fill(weights.begin(), weights.end(), 0.0f);
  std::fill(uvw, uvw + 3, 0.0);
}

BdaAverager::BdaAverager(const common::ParameterSet& parset,
                         const std::string& prefix,
                         const bool use_weights_and_flags)
    : timer_("BDA Averager"),
      bl_threshold_time_(parset.getDouble(prefix + "timebase", 0.0)),
      bl_threshold_channel_(parset.getDouble(prefix + "frequencybase", 0.0)),
      max_interval_(parset.getDouble(prefix + "maxinterval", 0.0)),
      min_channels_(parset.getUint(prefix + "minchannels", 1)),
      name_(prefix),
      maxfreqfactor_(1),
      maxtimefactor_(1),
      next_rownr_(0),
      bda_pool_size_(0),
      bda_buffer_(),
      baseline_buffers_(),
      expected_input_shape_(),
      use_weights_and_flags_(use_weights_and_flags) {}

BdaAverager::~BdaAverager() = default;

void BdaAverager::set_averaging_params(
    std::vector<unsigned int> baseline_factors,
    std::vector<std::vector<double>> freqs,
    std::vector<std::vector<double>> widths) {
  if (baseline_factors.empty() || freqs.empty() || widths.empty()) {
    throw std::invalid_argument(
        "One or more empty arguments while setting bda averaging parameters");
  }
  baseline_factors_ = std::move(baseline_factors);
  freqs_ = std::move(freqs);
  widths_ = std::move(widths);
}

void BdaAverager::updateInfo(const DPInfo& _info) {
  if (_info.nchan() != _info.chanFreqs().size() ||
      !_info.channelsAreRegular()) {
    throw std::invalid_argument("Invalid info in BDA averager");
  }
  if (!baseline_factors_.empty()) {
    if ((baseline_factors_.size() != _info.nbaselines()) ||
        (freqs_.size() != _info.nbaselines()) ||
        (widths_.size() != _info.nbaselines())) {
      throw std::invalid_argument(
          "Invalid averaging parameters in BDA averager");
    }
  }

  Step::updateInfo(_info);
  infoOut().setIsBDAIntervalFactorInteger(true);

  expected_input_shape_ = {_info.nbaselines(), _info.nchan(), _info.ncorr()};

  std::vector<std::vector<double>> freqs(_info.nbaselines());
  std::vector<std::vector<double>> widths(_info.nbaselines());

  const std::vector<double>& lengths = _info.getBaselineLengths();
  std::vector<unsigned int> baseline_factors;
  baseline_factors.reserve(_info.nbaselines());

  // Sum the relative number of channels of each baseline, for
  // determining the BDA output buffer size.
  float relative_channels = 0.0;

  // Track the maximum number of channels, as that should always fit.
  std::size_t max_channels = 0;

  // Apply the length thresholds to all baselines.
  baseline_buffers_.clear();
  baseline_buffers_.reserve(_info.nbaselines());

  for (std::size_t i = 0; i < _info.nbaselines(); ++i) {
    std::size_t factor_time;
    std::size_t nchan;

    if (baseline_factors_.empty()) {
      // Determine the time averaging factor. Ignore max_interval_ if it is
      // 0.0.
      factor_time = std::floor(bl_threshold_time_ / std::max(lengths[i], 0.1));
      if (max_interval_ > 0.0 &&
          factor_time * _info.timeInterval() > max_interval_) {
        factor_time = std::floor(max_interval_ / _info.timeInterval());
      }
      factor_time = std::max(factor_time, std::size_t{1});

      maxtimefactor_ = std::max(maxtimefactor_, factor_time);

      // Determine the number of channels in the output.
      nchan = _info.nchan();
      if (bl_threshold_channel_ > 0.0) {
        nchan = std::ceil(lengths[i] / bl_threshold_channel_ * double(nchan));
        if (nchan > _info.nchan()) {
          nchan = _info.nchan();
        } else if (nchan < min_channels_) {
          nchan = min_channels_;
        }
      }
    } else {
      factor_time = baseline_factors_[i];
      nchan = widths_[i].size();
    }
    baseline_factors.emplace_back(factor_time);
    baseline_buffers_.emplace_back(factor_time, _info.nchan(), nchan,
                                   _info.ncorr());
    relative_channels += float(nchan) / float(factor_time);
    max_channels = std::max(max_channels, nchan);

    // Calculate the center frequency and width of the averaged channels.
    const std::vector<std::size_t>& indices =
        baseline_buffers_.back().input_channel_indices;
    const auto input_widths_begin = _info.chanWidths().begin();
    freqs[i].reserve(nchan);
    widths[i].reserve(nchan);
    for (std::size_t ch = 0; ch < nchan; ++ch) {
      size_t freqfactor = indices[ch + 1] - indices[ch];
      if (freqfactor > maxfreqfactor_) {
        maxfreqfactor_ = freqfactor;
      }
      freqs[i].push_back(0.5 * (_info.chanFreqs()[indices[ch]] +
                                _info.chanFreqs()[indices[ch + 1] - 1]));
      widths[i].push_back(std::accumulate(input_widths_begin + indices[ch],
                                          input_widths_begin + indices[ch + 1],
                                          0.0));
    }
  }

  if (!baseline_factors_.empty()) {
    for (size_t i = 0; i < freqs_.size(); i++) {
      if (!common::EpsilonEqual(freqs_[i], freqs[i], 1.0e-3) ||
          !common::EpsilonEqual(widths_[i], widths[i], 1.0e-3)) {
        throw std::runtime_error(
            "Frequency averaging specified is not supported");
      }
    }
    freqs_.clear();
    widths_.clear();
  }

  std::size_t bda_channels = std::ceil(relative_channels);
  bda_channels = std::max(bda_channels, max_channels);

  bda_pool_size_ = _info.ncorr() * bda_channels;

  infoOut().update(std::move(baseline_factors));
  infoOut().setChannels(std::move(freqs), std::move(widths));
}

bool BdaAverager::process(std::unique_ptr<base::DPBuffer> buffer) {
  common::NSTimer::StartStop sstime(timer_);

  if (!bda_buffer_) {
    const std::vector<std::string> data_names = buffer->GetDataNames();

    if (baseline_buffers_.front().data.empty() && !data_names.empty()) {
      // In the first process() call, add data buffers to the baseline buffers.
      const std::complex<float> kZeroData(0.0f, 0.0f);
      for (BaselineBuffer& bb : baseline_buffers_) {
        const std::size_t data_size = bb.weights.size();
        for (const std::string& name : data_names) {
          // Use emplace_hint, since data_names is sorted.
          bb.data.emplace_hint(
              bb.data.end(), name,
              std::vector<std::complex<float>>(data_size, kZeroData));
        }
      }
    }

    SetBdaBuffer(data_names);
  }

  const DPBuffer::UvwType& uvw = buffer->GetUvw();
  const std::size_t n_correlations = getInfoOut().ncorr();

  if (buffer->GetData().shape() != expected_input_shape_ ||
      (use_weights_and_flags_ &&
       (buffer->GetFlags().shape() != expected_input_shape_ ||
        buffer->GetWeights().shape() != expected_input_shape_))) {
    throw std::runtime_error("BdaAverager: Invalid buffer shape");
  }

  for (std::size_t b = 0; b < baseline_buffers_.size(); ++b) {
    BaselineBuffer& bb = baseline_buffers_[b];
    ++bb.times_added;

    if (1 == bb.times_added) {
      bb.starttime = buffer->GetTime() - getInfoOut().timeInterval() / 2;
    }
    bb.interval += getInfoOut().timeInterval();
    bb.exposure += buffer->GetExposure();

    std::size_t offset = 0;

    if (!use_weights_and_flags_) {
      for (std::size_t och = 0; och < bb.input_channel_indices.size() - 1;
           ++och) {
        for (std::size_t ich = bb.input_channel_indices[och];
             ich < bb.input_channel_indices[och + 1]; ++ich) {
          for (std::size_t corr = 0; corr < n_correlations; ++corr) {
            for (auto& [data_name, data_vector] : bb.data) {
              data_vector[offset + corr] +=
                  buffer->GetData(data_name)(b, ich, corr);
            }
            bb.weights[offset + corr] += 1.0;
          }
        }
        offset += n_correlations;
      }
    } else {
      for (std::size_t och = 0; och < bb.input_channel_indices.size() - 1;
           ++och) {
        for (std::size_t ich = bb.input_channel_indices[och];
             ich < bb.input_channel_indices[och + 1]; ++ich) {
          for (std::size_t corr = 0; corr < n_correlations; ++corr) {
            if (!buffer->GetFlags()(b, ich, corr)) {
              const float weight = buffer->GetWeights()(b, ich, corr);
              for (auto& [data_name, data_vector] : bb.data) {
                data_vector[offset + corr] +=
                    buffer->GetData(data_name)(b, ich, corr) * weight;
              }
              bb.weights[offset + corr] += weight;
            }
          }
        }
        offset += n_correlations;
      }
    }
    bb.uvw[0] += uvw(b, 0);
    bb.uvw[1] += uvw(b, 1);
    bb.uvw[2] += uvw(b, 2);

    if (bb.times_added == bb.time_factor) {
      AddBaseline(b);  // BaselineBuffer is complete: Add it.
    }
  }

  // Send out full buffers immediately / don't wait for the next iteration
  // with detecting that there's no space left.
  if (0 == bda_buffer_->GetRemainingCapacity()) {
    getNextStep()->process(std::move(bda_buffer_));
    SetBdaBuffer(buffer->GetDataNames());
  }

  return true;
}

void BdaAverager::finish() {
  for (std::size_t b = 0; b < baseline_buffers_.size(); ++b) {
    if (baseline_buffers_[b].times_added > 0) {
      AddBaseline(b);
    }
  }

  if (bda_buffer_ && bda_buffer_->GetNumberOfElements() > 0) {
    getNextStep()->process(std::move(bda_buffer_));
  }
  bda_buffer_.reset();

  getNextStep()->finish();
}

void BdaAverager::AddBaseline(std::size_t baseline_nr) {
  BdaAverager::BaselineBuffer& bb = baseline_buffers_[baseline_nr];
  assert(bb.times_added > 0);
  const std::size_t n_channels = getInfoOut().chanFreqs(baseline_nr).size();

  // Divide data values by their total weight.
  float total_weight = 0.0f;
  for (std::size_t index = 0; index < bb.weights.size(); ++index) {
    const float weight = bb.weights[index];
    if (weight > 0) {
      for (auto& data_name_vector : bb.data) {
        data_name_vector.second[index] /= weight;
      }
    }
    total_weight += weight;
  }

  if (total_weight > 0) {
    const double factor = 1.0 / bb.times_added;
    bb.uvw[0] *= factor;
    bb.uvw[1] *= factor;
    bb.uvw[2] *= factor;
  }

  if (bda_buffer_->GetRemainingCapacity() < n_channels * getInfoOut().ncorr()) {
    const std::vector<std::string> data_names = bda_buffer_->GetDataNames();
    getNextStep()->process(std::move(bda_buffer_));
    SetBdaBuffer(data_names);
  }

  const size_t row = bda_buffer_->GetRows().size();
  bda_buffer_->AddRow(bb.starttime + bb.interval / 2, bb.interval, bb.exposure,
                      baseline_nr, n_channels, getInfoOut().ncorr(), nullptr,
                      nullptr, bb.weights.data(), bb.uvw);
  for (const auto& [data_name, data_vector] : bb.data) {
    std::copy(data_vector.begin(), data_vector.end(),
              bda_buffer_->GetData(row, data_name));
  }

  bb.Clear();  // Prepare baseline for a next iteration.
}

void BdaAverager::set_next_desired_buffersize(unsigned int buffer_size) {
  fixed_size_bda_buffers_.push(
      std::make_unique<BdaBuffer>(buffer_size, getProvidedFields()));
}

void BdaAverager::SetBdaBuffer(const std::vector<std::string>& data_names) {
  if (fixed_size_bda_buffers_.empty()) {
    bda_buffer_ =
        std::make_unique<BdaBuffer>(bda_pool_size_, getProvidedFields());
  } else {
    bda_buffer_ = std::move(fixed_size_bda_buffers_.front());
    fixed_size_bda_buffers_.pop();
  }
  for (const std::string& name : data_names) {
    bda_buffer_->AddData(name);
  }
}

void BdaAverager::show(std::ostream& os) const {
  os << "BdaAverager " << name_ << '\n';
  os << "  timebase:        " << bl_threshold_time_ << "s\n";
  os << "  max interval:    " << max_interval_ << "s\n";
  os << "  frequencybase:   " << bl_threshold_channel_ << '\n';
  os << "  min channels:    " << min_channels_ << "\n";
  os << "  max time factor: " << maxtimefactor_ << '\n';
  os << "  max freq factor: " << maxfreqfactor_ << '\n';
}

}  // namespace steps
}  // namespace dp3
