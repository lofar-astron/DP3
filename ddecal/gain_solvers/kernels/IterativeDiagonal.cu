// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeDiagonal.h"

#include <cuComplex.h>
#include <math_constants.h>

#include "Common.h"
#include "Complex.h"
#include "MatrixComplex2x2.h"

#define BLOCK_SIZE 128

template <bool Add>
__device__ void AddOrSubtract(size_t vis_index, size_t n_solutions,
                              const unsigned int* antenna_pairs,
                              const unsigned int* solution_map,
                              const cuDoubleComplex* solutions,
                              const cuM2x2FloatComplex* model,
                              const cuM2x2FloatComplex* residual_in,
                              cuM2x2FloatComplex* residual_out) {
  const uint32_t antenna_1 = antenna_pairs[vis_index * 2 + 0];
  const uint32_t antenna_2 = antenna_pairs[vis_index * 2 + 1];
  const size_t solution_index = solution_map[vis_index];
  const cuDoubleComplex* solution_1 =
      &solutions[(antenna_1 * n_solutions + solution_index) * 2];
  const cuDoubleComplex* solution_2 =
      &solutions[(antenna_2 * n_solutions + solution_index) * 2];

  const cuFloatComplex solution_1_0 = cuComplexDoubleToFloat(solution_1[0]);
  const cuFloatComplex solution_1_1 = cuComplexDoubleToFloat(solution_1[1]);
  const cuFloatComplex solution_2_0_conj =
      cuComplexDoubleToFloat(cuConj(solution_2[0]));
  const cuFloatComplex solution_2_1_conj =
      cuComplexDoubleToFloat(cuConj(solution_2[1]));

  const cuM2x2FloatComplex contribution(
      cuCmulf(cuCmulf(solution_1_0, model[vis_index][0]), solution_2_0_conj),
      cuCmulf(cuCmulf(solution_1_0, model[vis_index][1]), solution_2_1_conj),
      cuCmulf(cuCmulf(solution_1_1, model[vis_index][2]), solution_2_0_conj),
      cuCmulf(cuCmulf(solution_1_1, model[vis_index][3]), solution_2_1_conj));

  if (Add) {
    residual_out[vis_index] = residual_in[vis_index] + contribution;
  } else {
    residual_out[vis_index] = residual_in[vis_index] - contribution;
  }
}

__device__ void SolveDirection(size_t vis_index, size_t n_visibilities,
                               size_t n_direction_solutions, size_t n_solutions,
                               const unsigned int* antenna_pairs,
                               const unsigned int* solution_map,
                               const cuDoubleComplex* solutions,
                               const cuM2x2FloatComplex* model,
                               const cuM2x2FloatComplex* residual,
                               cuFloatComplex* numerator, float* denominator) {
  // Load correct variables to compute on.
  const size_t antenna_1 = antenna_pairs[vis_index * 2];
  const size_t antenna_2 = antenna_pairs[vis_index * 2 + 1];
  const size_t solution_index = solution_map[vis_index];

  const cuDoubleComplex* solution_antenna_1 =
      &solutions[(antenna_1 * n_solutions + solution_index) * 2];
  const cuDoubleComplex* solution_antenna_2 =
      &solutions[(antenna_2 * n_solutions + solution_index) * 2];

  const size_t rel_solution_index = solution_index - solution_map[0];

  // Calculate the contribution of this baseline for both antennas
  // For antenna2,
  // data_ba = data_ab^H, etc., therefore, numerator and denominator
  // become:
  // - num = data_ab^H * solutions_a * model_ab
  // - den = norm(model_ab^H * solutions_a)
  for (size_t i = 0; i < 2; i++) {
    const size_t antenna = antenna_pairs[vis_index * 2 + i];

    cuM2x2FloatComplex result;
    cuM2x2FloatComplex changed_model;

    if (i == 0) {
      const cuM2x2FloatComplexDiagonal solution(
          make_cuFloatComplex(solution_antenna_2[0].x, solution_antenna_2[0].y),
          make_cuFloatComplex(solution_antenna_2[1].x,
                              solution_antenna_2[1].y));
      changed_model = solution * cuConj(model[vis_index]);
      result = residual[vis_index] * changed_model;
    } else {
      const cuM2x2FloatComplexDiagonal solution(
          make_cuFloatComplex(solution_antenna_1[0].x, solution_antenna_1[0].y),
          make_cuFloatComplex(solution_antenna_1[1].x,
                              solution_antenna_1[1].y));
      changed_model = solution * model[vis_index];
      result = cuConj(residual[vis_index]) * changed_model;
    }

    const size_t full_solution_index =
        antenna * n_direction_solutions + rel_solution_index;

    // Atomic reduction into global memory
    atomicAdd(&numerator[full_solution_index * 2 + 0].x, result[0].x);
    atomicAdd(&numerator[full_solution_index * 2 + 0].y, result[0].y);
    atomicAdd(&numerator[full_solution_index * 2 + 1].x, result[3].x);
    atomicAdd(&numerator[full_solution_index * 2 + 1].y, result[3].y);

    atomicAdd(&denominator[full_solution_index * 2],
              cuNorm(changed_model[0]) + cuNorm(changed_model[2]));
    atomicAdd(&denominator[full_solution_index * 2 + 1],
              cuNorm(changed_model[1]) + cuNorm(changed_model[3]));
  }
}

__global__ void SolveDirectionKernel(
    size_t n_visibilities, size_t n_direction_solutions, size_t n_solutions,
    const unsigned int* antenna_pairs, const unsigned int* solution_map,
    const cuDoubleComplex* solutions, const cuM2x2FloatComplex* model,
    const cuM2x2FloatComplex* residual_in, cuM2x2FloatComplex* residual_temp,
    cuFloatComplex* numerator, float* denominator) {
  const size_t vis_index = blockIdx.x * blockDim.x + threadIdx.x;

  if (vis_index >= n_visibilities) {
    return;
  }

  AddOrSubtract<true>(vis_index, n_solutions, antenna_pairs, solution_map,
                      solutions, model, residual_in, residual_temp);

  SolveDirection(vis_index, n_visibilities, n_direction_solutions, n_solutions,
                 antenna_pairs, solution_map, solutions, model, residual_temp,
                 numerator, denominator);
}

void LaunchSolveDirectionKernel(
    cudaStream_t stream, size_t n_visibilities, size_t n_direction_solutions,
    size_t n_solutions, size_t direction, cu::DeviceMemory& antenna_pairs,
    cu::DeviceMemory& solution_map, cu::DeviceMemory& solutions,
    cu::DeviceMemory& model, cu::DeviceMemory& residual_in,
    cu::DeviceMemory& residual_temp, cu::DeviceMemory& numerator,
    cu::DeviceMemory& denominator) {
  const size_t block_dim = BLOCK_SIZE;
  const size_t grid_dim = (n_visibilities + block_dim) / block_dim;

  const size_t direction_offset = direction * n_visibilities;
  const unsigned int* solution_map_direction =
      Cast<const unsigned int>(solution_map) + direction_offset;
  const cuM2x2FloatComplex* model_direction =
      Cast<const cuM2x2FloatComplex>(model) + direction_offset;
  SolveDirectionKernel<<<grid_dim, block_dim, 0, stream>>>(
      n_visibilities, n_direction_solutions, n_solutions,
      Cast<const unsigned int>(antenna_pairs), solution_map_direction,
      Cast<const cuDoubleComplex>(solutions), model_direction,
      Cast<const cuM2x2FloatComplex>(residual_in),
      Cast<cuM2x2FloatComplex>(residual_temp), Cast<cuFloatComplex>(numerator),
      Cast<float>(denominator));
}

__global__ void SubtractKernel(size_t n_directions, size_t n_visibilities,
                               size_t n_solutions,
                               const unsigned int* antenna_pairs,
                               const unsigned int* solution_map,
                               const cuDoubleComplex* solutions,
                               const cuFloatComplex* model,
                               cuM2x2FloatComplex* residual) {
  const size_t vis_index = blockIdx.x * blockDim.x + threadIdx.x;

  if (vis_index >= n_visibilities) {
    return;
  }

  for (size_t direction = 0; direction < n_directions; direction++) {
    const size_t direction_offset = direction * n_visibilities;
    const unsigned int* solution_map_direction =
        solution_map + direction_offset;
    const cuFloatComplex* model_direction = model + (4 * direction_offset);
    AddOrSubtract<false>(
        vis_index, n_solutions, antenna_pairs, solution_map_direction,
        solutions, reinterpret_cast<const cuM2x2FloatComplex*>(model_direction),
        residual, residual);  // in-place
  }
}

void LaunchSubtractKernel(cudaStream_t stream, size_t n_directions,
                          size_t n_visibilities, size_t n_solutions,
                          cu::DeviceMemory& antenna_pairs,
                          cu::DeviceMemory& solution_map,
                          cu::DeviceMemory& solutions, cu::DeviceMemory& model,
                          cu::DeviceMemory& residual) {
  const size_t block_dim = BLOCK_SIZE;
  const size_t grid_dim = (n_visibilities + block_dim) / block_dim;

  SubtractKernel<<<grid_dim, block_dim, 0, stream>>>(
      n_directions, n_visibilities, n_solutions,
      Cast<const unsigned int>(antenna_pairs),
      Cast<const unsigned int>(solution_map),
      Cast<const cuDoubleComplex>(solutions), Cast<const cuFloatComplex>(model),
      Cast<cuM2x2FloatComplex>(residual));
}

__global__ void SolveNextSolutionKernel(unsigned int n_antennas,
                                        unsigned int n_direction_solutions,
                                        const unsigned int n_solutions,
                                        const unsigned int* solution_map,
                                        const cuFloatComplex* numerator,
                                        const float* denominator,
                                        cuDoubleComplex* next_solutions) {
  const size_t antenna = blockIdx.x * blockDim.x + threadIdx.x;

  if (antenna >= n_antennas) {
    return;
  }

  for (size_t relative_solution = 0; relative_solution < n_direction_solutions;
       relative_solution++) {
    const size_t solution_index = relative_solution + solution_map[0];
    cuDoubleComplex* destination =
        &next_solutions[(antenna * n_solutions + solution_index) * 2];
    const size_t index = antenna * n_direction_solutions + relative_solution;

    for (size_t pol = 0; pol < 2; pol++) {
      if (denominator[index * 2 + pol] == 0.0) {
        destination[pol] = {CUDART_NAN, CUDART_NAN};
      } else {
        // The CPU code performs this compuation in double-precision,
        // however single-precision also seems sufficiently accurate.
        destination[pol] = {
            numerator[index * 2 + pol].x / denominator[index * 2 + pol],
            numerator[index * 2 + pol].y / denominator[index * 2 + pol]};
      }
    }
  }
}

void LaunchSolveNextSolutionKernel(
    cudaStream_t stream, size_t n_antennas, size_t n_visibilities,
    size_t n_direction_solutions, size_t n_solutions, size_t direction,
    cu::DeviceMemory& antenna_pairs, cu::DeviceMemory& solution_map,
    cu::DeviceMemory& next_solutions, cu::DeviceMemory& numerator,
    cu::DeviceMemory& denominator) {
  const size_t block_dim = BLOCK_SIZE;
  const size_t grid_dim = (n_antennas + block_dim) / block_dim;

  const size_t direction_offset = direction * n_visibilities;
  const unsigned int* solution_map_direction =
      Cast<const unsigned int>(solution_map) + direction_offset;
  SolveNextSolutionKernel<<<grid_dim, block_dim, 0, stream>>>(
      n_antennas, n_direction_solutions, n_solutions, solution_map_direction,
      Cast<const cuFloatComplex>(numerator), Cast<const float>(denominator),
      Cast<cuDoubleComplex>(next_solutions));
}

__global__ void StepKernel(const size_t n_visibilities,
                           const cuDoubleComplex* solutions,
                           cuDoubleComplex* next_solutions, bool phase_only,
                           double step_size) {
  const size_t vis_index = blockIdx.x * blockDim.x + threadIdx.x;

  if (vis_index >= n_visibilities) {
    return;
  }

  if (phase_only) {
    // In phase only mode, a step is made along the complex circle,
    // towards the shortest direction.
    double phase_from = cuCarg(solutions[vis_index]);
    double distance = cuCarg(next_solutions[vis_index]) - phase_from;
    if (distance > CUDART_PI)
      distance = distance - 2.0 * CUDART_PI;
    else if (distance < -CUDART_PI)
      distance = distance + 2.0 * CUDART_PI;

    next_solutions[vis_index] =
        cuCpolar(1.0, phase_from + step_size * distance);
  } else {
    next_solutions[vis_index] =
        cuCadd(cuCmul(solutions[vis_index], (1.0 - step_size)),
               cuCmul(next_solutions[vis_index], step_size));
  }
}

void LaunchStepKernel(cudaStream_t stream, size_t n_visibilities,
                      cu::DeviceMemory& solutions,
                      cu::DeviceMemory& next_solutions, bool phase_only,
                      double step_size) {
  const size_t block_dim = BLOCK_SIZE;
  const size_t grid_dim = (n_visibilities + block_dim) / block_dim;

  StepKernel<<<grid_dim, block_dim, 0, stream>>>(
      n_visibilities, Cast<const cuDoubleComplex>(solutions),
      Cast<cuDoubleComplex>(next_solutions), phase_only, step_size);
}
