// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverBuffer.h"

#include <cassert>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <dp3/base/DPBuffer.h>

using dp3::base::DPBuffer;

namespace {
const size_t kNCorrelations = 4;
}  // namespace

namespace dp3 {
namespace ddecal {

SolverBuffer::SolverBuffer() : data_(), model_buffers_() {}

void SolverBuffer::AssignAndWeight(
    const std::vector<DPBuffer>& unweighted_data_buffers,
    std::vector<std::vector<DPBuffer>>&& model_buffers) {
  const std::size_t n_times = model_buffers.size();
  assert(unweighted_data_buffers.size() >= n_times);
  data_.resize(n_times);
  model_buffers_ = std::move(model_buffers);

  for (std::size_t timestep = 0; timestep != n_times; ++timestep) {
    const aocommon::xt::Span<std::complex<float>, 3>& unweighted_data =
        unweighted_data_buffers[timestep].GetData();
    const aocommon::xt::Span<float, 3>& weights =
        unweighted_data_buffers[timestep].GetWeights();
    assert(unweighted_data.shape() == weights.shape());
    assert(timestep == 0 ||
           unweighted_data.shape() ==
               unweighted_data_buffers.front().GetData().shape());
    assert(kNCorrelations == unweighted_data.shape(2));

    // Flag all non-finite values in the data and model data buffers.
    // There is one flag for all correlations, so if the data for one
    // correlation has a non-finite value, the other correlations are also
    // flagged. Although the initialization to false can be avoided by directly
    // assigning flags for the first correlation, a benchmark showed that the
    // (simpler) code below gives slightly better performance.
    xt::xtensor<bool, 2> flags(
        {unweighted_data.shape(0), unweighted_data.shape(1)}, false);
    for (std::size_t correlation = 0; correlation < kNCorrelations;
         ++correlation) {
      flags |= !xt::isfinite(
          xt::view(unweighted_data, xt::all(), xt::all(), correlation));
      for (DPBuffer& model_buffer : model_buffers_[timestep]) {
        flags |= !xt::isfinite(xt::view(model_buffer.GetData(), xt::all(),
                                        xt::all(), correlation));
      }
    }

    // Copy and weigh the data. Weigh the model data.
    // If the flag is set, set both the data and model data to zero.
    // Storing the result in an xtensor (and not in an expression) ensures
    // that the expression, including the sqrt() calls, is evaluated once.
    const xt::xtensor<float, 3> weights_sqrt =
        xt::where(xt::view(flags, xt::all(), xt::all(), xt::newaxis()), 0.0f,
                  xt::sqrt(weights));

    data_[timestep] = unweighted_data * weights_sqrt;
    for (DPBuffer& model_buffer : model_buffers_[timestep]) {
      model_buffer.GetData() *= weights_sqrt;
    }
  }
}

const std::complex<float>* SolverBuffer::DataPointer(size_t time_index,
                                                     size_t baseline,
                                                     size_t channel) const {
  return &data_[time_index](baseline, channel, 0);
}

const std::complex<float>* SolverBuffer::ModelDataPointer(
    size_t time_index, size_t direction, size_t baseline,
    size_t channel) const {
  const DPBuffer& buffer = model_buffers_[time_index][direction];
  return &buffer.GetData()(baseline, channel, 0);
}

}  // namespace ddecal
}  // namespace dp3
