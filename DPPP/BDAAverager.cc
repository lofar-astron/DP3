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
#include "../Common/ParameterSet.h"

#include <boost/make_unique.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>

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

BDAAverager::BDAAverager(double bl_threshold_time, double bl_threshold_channel)
    : bl_threshold_time_(bl_threshold_time),
      bl_threshold_channel_(bl_threshold_channel),
      next_rownr_(0),
      bda_pool_size_(0),
      bda_buffer_(),
      baseline_buffers_() {}

BDAAverager::~BDAAverager() {}

void BDAAverager::updateInfo(const DPInfo& info_in) {
  DPStep::updateInfo(info_in);

  const std::vector<double>& lengths = info_in.getBaselineLengths();
  std::vector<size_t> baseline_factors(info_in.nbaselines());

  // Sum the relative number of channels of each baseline, for
  // determining the BDA output buffer size.
  float relative_channels = 0.0;

  // Track the maximum number of channels, as that should always fit.
  std::size_t max_channels = 0;

  // Apply the length thresholds to all baselines.
  baseline_buffers_.clear();
  baseline_buffers_.reserve(info_in.nbaselines());
  for (std::size_t i = 0; i < info_in.nbaselines(); ++i) {
    std::size_t factor_time = std::floor(bl_threshold_time_ / lengths[i]);
    factor_time = std::max(factor_time, std::size_t(1));

    // Determine the number of channels in the output.
    std::size_t nchan =
        std::ceil(lengths[i] / bl_threshold_channel_ * info_in.nchan());
    if (nchan > info_in.nchan()) {
      nchan = info_in.nchan();
    } else if (nchan < 1) {
      nchan = 1;
    }

    baseline_factors.emplace_back(factor_time);
    baseline_buffers_.emplace_back(factor_time, info_in.nchan(), nchan,
                                   info_in.ncorr());
    relative_channels += float(nchan) / float(factor_time);
    max_channels = std::max(max_channels, nchan);
  }

  std::size_t bda_channels = std::ceil(relative_channels);
  bda_channels = std::max(bda_channels, max_channels);

  bda_pool_size_ = info_in.ncorr() * bda_channels;
  bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);

  info().setBDAFactors(baseline_factors);
}

bool BDAAverager::process(const DPBuffer& buffer) {
  const casacore::IPosition& shape = buffer.getData().shape();
  if (shape.size() != 3u || shape[0] != info().ncorr() ||
      shape[1] != info().nchan() ||
      shape[2] != ssize_t(baseline_buffers_.size())) {
    throw std::runtime_error("Invalid buffer shape");
  }

  assert(bda_buffer_);

  for (std::size_t b = 0; b < baseline_buffers_.size(); ++b) {
    BaselineBuffer& bb = baseline_buffers_[b];
    ++bb.times_added;

    if (1 == bb.times_added) {
      bb.time = buffer.getTime();
    }
    bb.interval += buffer.getExposure();

    std::complex<float>* data = bb.data.data();
    float* weights = bb.weights.data();
    float total_weight = 0.0f;

    for (std::size_t och = 0; och < bb.input_channel_indices.size() - 1;
         ++och) {
      for (std::size_t ich = bb.input_channel_indices[och];
           ich < bb.input_channel_indices[och + 1]; ++ich) {
        for (std::size_t corr = 0; corr < info().ncorr(); ++corr) {
          if (!buffer.getFlags()(corr, ich, b)) {
            const float weight = buffer.getWeights()(corr, ich, b);

            data[corr] += buffer.getData()(corr, ich, b) * weight;
            weights[corr] += weight;
            total_weight += weight;
          }
        }
      }
      data += info().ncorr();
      weights += info().ncorr();
    }

    bb.uvw[0] += buffer.getUVW()(0, b) * total_weight;
    bb.uvw[1] += buffer.getUVW()(1, b) * total_weight;
    bb.uvw[2] += buffer.getUVW()(2, b) * total_weight;

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
  const std::size_t nchan = bb.input_channel_indices.size() - 1;

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

  bda_buffer_->AddRow(bb.time, bb.interval, next_rownr_, baseline_nr, nchan,
                      info().ncorr(), bb.data.data(), nullptr,
                      bb.weights.data(), nullptr, bb.uvw);
}

void BDAAverager::show(std::ostream& stream) const { stream << "BDAAverager"; }

}  // namespace DPPP
}  // namespace DP3
