// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_TOOLS_H_
#define DDECAL_SOLVER_TOOLS_H_

#include <memory>
#include <vector>

#include <dp3/base/DPBuffer.h>

namespace dp3::ddecal {

/**
 * Weighs input data by multiplying it with weight values.
 * All arguments must have the same shape.
 *
 * This function is a manually simdized version of
 * xt::noalias(out) = in * weights;
 * (xt::noalias avoids an intermediate buffer for the result.)
 *
 * On an Intel Xeon E5-2660, this function is ~3 x faster than
 * xt::noalias(out) = in * weights;
 * when keep_unweighted_model_data == false in AssignAndWeight.
 *
 * @param in Input data buffer.
 * @param out Output data buffer. May be equal to 'in'.
 * @param weights Weights buffer.
 */
void Weigh(const base::DPBuffer::DataType& in, base::DPBuffer::DataType& out,
           const base::DPBuffer::WeightsType& weights);

/**
 * This function takes (unweighted) data and model data, as well as a
 * weights array, weights these and writes the result.
 *
 * @param unweighted_buffers A vector with one buffer for each timestep. Each
 * buffer should contain unweighted data and the corresponding weight values.
 * @param direction_names A list with the names of the model data buffers
 * in the DPBuffers.
 * @param weighted_buffers A vector with one buffer for each timestep.
 * This function writes the weighted data to these buffers, using the same
 * direction name. Existing data with that name is overwritten.
 * @param keep_unweighted_model_data If true, keep the model data in
 * unweighted_buffers and create new model data buffers in weighted_buffers.
 * If false, avoid creating new model data buffers by moving the model data
 * from unweighted_buffers to weighted_buffers and weighing the data in-place.
 */
void AssignAndWeight(
    std::vector<std::unique_ptr<base::DPBuffer>>& unweighted_buffers,
    const std::vector<std::string>& direction_names,
    std::vector<base::DPBuffer>& weighted_buffers,
    bool keep_unweighted_model_data);

}  // namespace dp3::ddecal

#endif  // DDECAL_SOLVER_TOOLS_H
