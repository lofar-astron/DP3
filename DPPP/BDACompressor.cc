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

#include "BDACompressor.h"

#include "BDABuffer.h"
#include "../Common/ParameterSet.h"

#include <boost/make_unique.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>

namespace DP3 {
namespace DPPP {

BDACompressor::Baseline::Baseline(std::size_t _factor, std::size_t _n_channels,
                                  std::size_t _n_correlations)
    : added(0),
      factor(_factor),
      data(_n_channels * _n_correlations, {0.0f, 0.0f}),
      weights(_n_channels * _n_correlations, 0.0f),
      summed_weight(0.0f),
      uvw{0.0f, 0.0f, 0.0f} {}

BDACompressor::BDACompressor()
    : next_rownr_(0), bda_pool_size_(0), bda_buffer_(), baselines_() {}

BDACompressor::~BDACompressor() {}

void BDACompressor::updateInfo(const DPInfo& info) {
  DPStep::updateInfo(info);

  const std::vector<double>& lengths = info.getBaselineLengths();

  // Determine the averaging threshold.
  // For baselines longer than the threshold, the averaging factor is 1, which
  // means no averaging.
  // For baselines shorter than the threshold, the averaging factor is
  // threshold / length, rounded down, with a minimum of 1.

  // For now, the threshold is the length of the longest baseline divided by 4.
  // In the future, the parset should configure this threshold factor.
  const double kThresholdFactor = 4.0;
  const double threshold =
      *std::max_element(lengths.begin(), lengths.end()) / kThresholdFactor;

  // Apply the threshold to all baselines.
  baselines_.clear();
  baselines_.reserve(info.nbaselines());
  for (std::size_t i = 0; i < info.nbaselines(); ++i) {
    std::size_t factor = std::floor(threshold / lengths[i]);
    factor = std::min(factor, std::size_t(1));
    baselines_.emplace_back(factor, info.nchan(), info.ncorr());
  }

  // BDA buffers will hold kBdaRatio times less elements than 'buffer'.
  const std::size_t kBdaRatio = 8;
  bda_pool_size_ = info.ncorr() * info.nchan() * baselines_.size() / kBdaRatio;
  bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
}

bool BDACompressor::process(const DPBuffer& buffer) {
  const std::size_t buffer_size = buffer.getData().size();
  if (buffer_size != info().ncorr() * info().nchan() * baselines_.size()) {
    throw std::runtime_error("Invalid buffer size");
  }

  assert(bda_buffer_);

  for (std::size_t b = 0; b < baselines_.size(); ++b) {
    Baseline& baseline = baselines_[b];
    std::complex<float>* data = baseline.data.data();
    float* weights = baseline.weights.data();
    float total_weight = 0.0;
    for (std::size_t ch = 0; ch < info().nchan(); ++ch) {
      for (std::size_t corr = 0; corr < info().ncorr(); ++corr) {
        if (!buffer.getFlags()(corr, ch, b)) {
          const float weight = buffer.getWeights()(corr, ch, b);

          *data += buffer.getData()(corr, ch, b) * weight;
          *weights += weight;
          total_weight += weight;
        }
      }
    }
    baseline.uvw[0] += buffer.getUVW()(0, b) * total_weight;
    baseline.uvw[1] += buffer.getUVW()(1, b) * total_weight;
    baseline.uvw[2] += buffer.getUVW()(2, b) * total_weight;
    baseline.summed_weight += total_weight;

    // Check if the BDA baseline is complete.
    ++baseline.added;
    if (baseline.added == baseline.factor) {
      const float factor = 1.0f / baseline.factor;
      for (std::complex<float>& d : baseline.data) {
        d *= factor;
      }
      for (float& w : baseline.weights) {
        w *= factor;
      }

      baseline.uvw[0] /= baseline.summed_weight;
      baseline.uvw[1] /= baseline.summed_weight;
      baseline.uvw[2] /= baseline.summed_weight;

      const double time =
          buffer.getTime() - ((baseline.factor - 1) * buffer.getExposure());
      const double interval = buffer.getExposure() * baseline.factor;

      if (!bda_buffer_->AddRow(time, interval, next_rownr_, b, info().nchan(),
                               info().ncorr(), baseline.data.data(), nullptr,
                               baseline.weights.data(), nullptr,
                               baseline.uvw)) {
        // BDA buffer is full. Send it away and create a new one.
        getNextStep()->process(std::move(bda_buffer_));
        bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);

        if (!bda_buffer_->AddRow(time, interval, next_rownr_, b, info().nchan(),
                                 info().ncorr(), baseline.data.data(), nullptr,
                                 baseline.weights.data(), nullptr,
                                 baseline.uvw)) {
          throw std::runtime_error("Empty BDA buffer has no space");
        }
      }

      baseline.added = 0;
    }
  }

  return true;
}

void BDACompressor::finish() {
  getNextStep()->process(std::move(bda_buffer_));
  getNextStep()->finish();
}

void BDACompressor::show(std::ostream& stream) const {
  stream << "BDACompressor";
}

}  // namespace DPPP
}  // namespace DP3
