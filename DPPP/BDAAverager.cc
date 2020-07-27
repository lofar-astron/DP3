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

BDAAverager::Baseline::Baseline(std::size_t _time_factor,
                                std::size_t _n_channels,
                                std::size_t _n_correlations)
    : added(0),
      time_factor(_time_factor),
      time(0.0),
      interval(0.0),
      data(_n_channels * _n_correlations, {0.0f, 0.0f}),
      weights(_n_channels * _n_correlations, 0.0f),
      summed_weight(0.0f),
      uvw{0.0f, 0.0f, 0.0f} {}

void BDAAverager::Baseline::Clear() {
  added = 0;
  time = 0.0;
  interval = 0.0;
  std::fill(data.begin(), data.end(), std::complex<float>{0.0f, 0.0f});
  std::fill(weights.begin(), weights.end(), 0.0f);
  summed_weight = 0.0f;
  std::fill(uvw, uvw + 3, 0.0);
}

BDAAverager::BDAAverager()
    : next_rownr_(0), bda_pool_size_(0), bda_buffer_(), baselines_() {}

BDAAverager::~BDAAverager() {}

void BDAAverager::updateInfo(const DPInfo& info) {
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

  // Sum the relative contribution of each baseline to the output, for
  // determining the number of rows in the BDA output buffers.
  double bda_baselines = 0.0;

  // Apply the threshold to all baselines.
  baselines_.clear();
  baselines_.reserve(info.nbaselines());
  for (std::size_t i = 0; i < info.nbaselines(); ++i) {
    std::size_t factor = std::floor(threshold / lengths[i]);
    factor = std::max(factor, std::size_t(1));
    baselines_.emplace_back(factor, info.nchan(), info.ncorr());
    bda_baselines += 1.0 / factor;
  }

  const std::size_t bda_rows = std::ceil(bda_baselines);
  bda_pool_size_ = info.ncorr() * info.nchan() * bda_rows;
  bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
}

bool BDAAverager::process(const DPBuffer& buffer) {
  const std::size_t buffer_size = buffer.getData().size();
  if (buffer_size != info().ncorr() * info().nchan() * baselines_.size()) {
    throw std::runtime_error("Invalid buffer size");
  }

  assert(bda_buffer_);

  for (std::size_t b = 0; b < baselines_.size(); ++b) {
    Baseline& baseline = baselines_[b];
    ++baseline.added;

    if (1 == baseline.added) {
      baseline.time = buffer.getTime();
    }
    baseline.interval += buffer.getExposure();

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

          ++data;
          ++weights;
        }
      }
    }
    baseline.uvw[0] += buffer.getUVW()(0, b) * total_weight;
    baseline.uvw[1] += buffer.getUVW()(1, b) * total_weight;
    baseline.uvw[2] += buffer.getUVW()(2, b) * total_weight;
    baseline.summed_weight += total_weight;

    if (baseline.added == baseline.time_factor) {
      AddBaseline(b);    // Baseline is complete: Add it.
      baseline.Clear();  // Prepare baseline for the next iteration.
    }
  }

  // Send out full buffers immediately / don't wait for the next iteration
  // with detecting that there's no space left.
  if (0 == bda_buffer_->GetRemainingCapacity()) {
    getNextStep()->process(std::move(bda_buffer_));
    bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
  }

  return true;
}

void BDAAverager::finish() {
  assert(bda_buffer_);

  for (std::size_t b = 0; b < baselines_.size(); ++b) {
    if (baselines_[b].added > 0) {
      AddBaseline(b);
      baselines_[b].Clear();
    }
  }

  if (bda_buffer_->GetNumberOfElements() > 0) {
    getNextStep()->process(std::move(bda_buffer_));
  }
  bda_buffer_.reset();

  getNextStep()->finish();
}

void BDAAverager::AddBaseline(std::size_t baseline_nr) {
  BDAAverager::Baseline& baseline = baselines_[baseline_nr];
  assert(baseline.added > 0);

  // Divide data and weight by the number of added values.
  if (baseline.added != 1u) {
    const float factor = 1.0f / baseline.added;
    for (std::complex<float>& d : baseline.data) {
      d *= factor;
    }
    for (float& w : baseline.weights) {
      w *= factor;
    }
  }

  baseline.uvw[0] /= baseline.summed_weight;
  baseline.uvw[1] /= baseline.summed_weight;
  baseline.uvw[2] /= baseline.summed_weight;

  if (bda_buffer_->GetRemainingCapacity() < info().nchan() * info().ncorr()) {
    getNextStep()->process(std::move(bda_buffer_));
    bda_buffer_ = boost::make_unique<BDABuffer>(bda_pool_size_);
  }

  bda_buffer_->AddRow(baseline.time, baseline.interval, next_rownr_,
                      baseline_nr, info().nchan(), info().ncorr(),
                      baseline.data.data(), nullptr, baseline.weights.data(),
                      nullptr, baseline.uvw);
}

void BDAAverager::show(std::ostream& stream) const { stream << "BDAAverager"; }

}  // namespace DPPP
}  // namespace DP3
