// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SolverTools.h"

#include <cassert>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace {
const size_t kNCorrelations = 4;
}  // namespace

namespace dp3::ddecal {

void AssignAndWeight(
    std::vector<std::unique_ptr<base::DPBuffer>>& unweighted_buffers,
    const std::vector<std::string>& direction_names,
    std::vector<base::DPBuffer>& weighted_buffers,
    bool keep_unweighted_model_data) {
  const std::size_t n_times = unweighted_buffers.size();
  assert(weighted_buffers.size() >= n_times);

  for (std::size_t timestep = 0; timestep != n_times; ++timestep) {
    base::DPBuffer& unweighted_buffer = *unweighted_buffers[timestep];
    base::DPBuffer& weighted_buffer = weighted_buffers[timestep];
    const aocommon::xt::Span<std::complex<float>, 3>& unweighted_data =
        unweighted_buffer.GetData();
    const aocommon::xt::Span<float, 3>& weights =
        unweighted_buffer.GetWeights();
    assert(unweighted_data.shape() == weights.shape());
    assert(timestep == 0 || unweighted_data.shape() ==
                                unweighted_buffers.front()->GetData().shape());
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
      for (const std::string& name : direction_names) {
        flags |= !xt::isfinite(xt::view(unweighted_buffer.GetData(name),
                                        xt::all(), xt::all(), correlation));
      }
    }

    // Copy and weigh the data. Weigh the model data.
    // If the flag is set, set both the data and model data to zero.
    // Storing the result in an xtensor (and not in an expression) ensures
    // that the expression, including the sqrt() calls, is evaluated once.
    const xt::xtensor<float, 3> weights_sqrt =
        xt::where(xt::view(flags, xt::all(), xt::all(), xt::newaxis()), 0.0f,
                  xt::sqrt(weights));

    weighted_buffer.ResizeData(unweighted_data.shape());
    weighted_buffer.GetData() = unweighted_data * weights_sqrt;
    for (const std::string& name : direction_names) {
      if (keep_unweighted_model_data) {
        if (!weighted_buffer.HasData(name)) weighted_buffer.AddData(name);
        weighted_buffer.GetData(name) =
            unweighted_buffer.GetData(name) * weights_sqrt;
      } else {
        weighted_buffer.MoveData(unweighted_buffer, name, name);
        weighted_buffer.GetData(name) *= weights_sqrt;
      }
    }
  }
}

}  // namespace dp3::ddecal
