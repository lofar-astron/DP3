// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include "BDAAverager.h"

#include "BDABuffer.h"
#include "DPInput.h"
#include "../Common/ParameterSet.h"

#include <boost/make_unique.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace DP3 {
namespace DPPP {

BDAAverager::BaselineBuffer::BaselineBuffer(std::size_t _time_factor,
                                            std::size_t n_input_channels,
                                            std::size_t n_output_channels,
                                            std::size_t n_correlations)
    : times_added(0),
      time_factor(_time_factor),
      input_channel_indices(),
      time(0.0),
      interval(0.0),
      data(n_output_channels * n_correlations, {0.0f, 0.0f}),
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

void BDAAverager::BaselineBuffer::Clear() {
  times_added = 0;
  time = 0.0;
  interval = 0.0;
  std::fill(data.begin(), data.end(), std::complex<float>{0.0f, 0.0f});
  std::fill(weights.begin(), weights.end(), 0.0f);
  std::fill(uvw, uvw + 3, 0.0);
}

BDAAverager::BDAAverager(DPInput& input, const DP3::ParameterSet& parset,
                         const std::string& prefix)
    : input_(input),
      timer_("BDA Averager"),
      bl_threshold_time_(parset.getDouble(prefix + "timethresholdlength", 0.0)),
      bl_threshold_channel_(
          parset.getDouble(prefix + "freqthresholdlength", 0.0)),
      max_interval_(parset.getDouble(prefix + "maxinterval", 0.0)),
      min_channels_(parset.getUint(prefix + "minchannels", 1)),
      next_rownr_(0),
      bda_pool_size_(0),
      bda_buffer_(),
      baseline_buffers_() {}

BDAAverager::~BDAAverager() {}

void BDAAverager::updateInfo(const DPInfo& _info) {
  if (_info.nchan() != _info.chanFreqs().size() ||
      !_info.channelsAreRegular()) {
    throw std::invalid_argument("Invalid info in BDA averager");
  }

  DPStep::updateInfo(_info);
  info().setNeedVisData();

  std::vector<std::vector<double>> freqs(_info.nbaselines());
  std::vector<std::vector<double>> widths(_info.nbaselines());

  const std::vector<double>& lengths = _info.getBaselineLengths();

  // Sum the relative number of channels of each baseline, for
  // determining the BDA output buffer size.
  float relative_channels = 0.0;

  // Track the maximum number of channels, as that should always fit.
  std::size_t max_channels = 0;

  // Apply the length thresholds to all baselines.
  baseline_buffers_.clear();
  baseline_buffers_.reserve(_info.nbaselines());
  for (std::size_t i = 0; i < _info.nbaselines(); ++i) {
    // Determine the time averaging factor. Ignore max_interval_ if it is 0.0.
    std::size_t factor_time = std::floor(bl_threshold_time_ / lengths[i]);
    if (max_interval_ > 0.0 &&
        factor_time * _info.timeInterval() > max_interval_) {
      factor_time = std::floor(max_interval_ / _info.timeInterval());
    }
    factor_time = std::max(factor_time, std::size_t{1});

    // Determine the number of channels in the output.
    std::size_t nchan = _info.nchan();
    if (bl_threshold_channel_ > 0.0) {
      nchan = std::ceil(lengths[i] / bl_threshold_channel_ * double(nchan));
      if (nchan > _info.nchan()) {
        nchan = _info.nchan();
      } else if (nchan < min_channels_) {
        nchan = min_channels_;
      }
    }

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
      freqs[i].push_back(0.5 * (_info.chanFreqs()[indices[ch]] +
                                _info.chanFreqs()[indices[ch + 1] - 1]));
      widths[i].push_back(std::accumulate(input_widths_begin + indices[ch],
                                          input_widths_begin + indices[ch + 1],
                                          0.0));
    }
  }

  std::size_t bda_channels = std::ceil(relative_channels);
  bda_channels = std::max(bda_channels, max_channels);

  bda_pool_size_ = _info.ncorr() * bda_channels;
  bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);

  info().set(std::move(freqs), std::move(widths));
}

bool BDAAverager::process(const DPBuffer& buffer) {
  NSTimer::StartStop sstime(timer_);

  DPBuffer dummy_buffer;
  const casacore::Cube<float>& weights =
      input_.fetchWeights(buffer, dummy_buffer, timer_);
  const casacore::Matrix<double>& uvw =
      input_.fetchUVW(buffer, dummy_buffer, timer_);

  const casacore::IPosition expected_shape(3, info().ncorr(), info().nchan(),
                                           baseline_buffers_.size());
  if (buffer.getData().shape() != expected_shape ||
      buffer.getFlags().shape() != expected_shape ||
      weights.shape() != expected_shape) {
    throw std::runtime_error("BDAAverager: Invalid buffer shape");
  }

  assert(bda_buffer_);

  for (std::size_t b = 0; b < baseline_buffers_.size(); ++b) {
    BaselineBuffer& bb = baseline_buffers_[b];
    ++bb.times_added;

    if (1 == bb.times_added) {
      bb.time = buffer.getTime();
    }
    bb.interval += buffer.getExposure();

    std::complex<float>* bb_data = bb.data.data();
    float* bb_weights = bb.weights.data();
    float total_weight = 0.0f;

    for (std::size_t och = 0; och < bb.input_channel_indices.size() - 1;
         ++och) {
      for (std::size_t ich = bb.input_channel_indices[och];
           ich < bb.input_channel_indices[och + 1]; ++ich) {
        for (std::size_t corr = 0; corr < info().ncorr(); ++corr) {
          if (!buffer.getFlags()(corr, ich, b)) {
            const float weight = weights(corr, ich, b);

            bb_data[corr] += buffer.getData()(corr, ich, b) * weight;
            bb_weights[corr] += weight;
            total_weight += weight;
          }
        }
      }
      bb_data += info().ncorr();
      bb_weights += info().ncorr();
    }

    bb.uvw[0] += uvw(0, b) * total_weight;
    bb.uvw[1] += uvw(1, b) * total_weight;
    bb.uvw[2] += uvw(2, b) * total_weight;

    if (bb.times_added == bb.time_factor) {
      AddBaseline(b);  // BaselineBuffer is complete: Add it.
      bb.Clear();      // Prepare baseline for the next iteration.
    }
  }

  // Send out full buffers immediately / don't wait for the next iteration
  // with detecting that there's no space left.
  if (0 == bda_buffer_->GetRemainingCapacity()) {
    getNextStep()->process(std::move(bda_buffer_));
    bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
  }

  return true;
}  // namespace DPPP

void BDAAverager::finish() {
  assert(bda_buffer_);

  for (std::size_t b = 0; b < baseline_buffers_.size(); ++b) {
    if (baseline_buffers_[b].times_added > 0) {
      AddBaseline(b);
      baseline_buffers_[b].Clear();
    }
  }

  if (bda_buffer_->GetNumberOfElements() > 0) {
    getNextStep()->process(std::move(bda_buffer_));
  }
  bda_buffer_.reset();

  getNextStep()->finish();
}

void BDAAverager::AddBaseline(std::size_t baseline_nr) {
  BDAAverager::BaselineBuffer& bb = baseline_buffers_[baseline_nr];
  assert(bb.times_added > 0);
  const std::size_t nchan = info().chanFreqs(baseline_nr).size();

  // Divide data values by their total weight.
  float* weights = bb.weights.data();
  float total_weight = 0.0f;
  for (std::complex<float>& d : bb.data) {
    d /= *weights;
    total_weight += *weights;
    ++weights;
  }

  const double factor = 1.0 / total_weight;
  bb.uvw[0] *= factor;
  bb.uvw[1] *= factor;
  bb.uvw[2] *= factor;

  if (bda_buffer_->GetRemainingCapacity() < nchan * info().ncorr()) {
    getNextStep()->process(std::move(bda_buffer_));
    bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
  }

  bda_buffer_->AddRow(bb.time, bb.interval, next_rownr_, baseline_nr,
                      nchan * info().ncorr(), bb.data.data(), nullptr,
                      bb.weights.data(), nullptr, bb.uvw);
}

void BDAAverager::show(std::ostream& stream) const { stream << "BDAAverager"; }

}  // namespace DPPP
}  // namespace DP3
