// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_GAIN_SOLVERS_KERNELS_ITERATIVEDIAGONAL_H_
#define DP3_DDECAL_GAIN_SOLVERS_KERNELS_ITERATIVEDIAGONAL_H_

#include <complex>
#include <cuda_runtime.h>

#include <cudawrappers/cu.hpp>

void LaunchSubtractKernel(cudaStream_t stream, size_t n_directions,
                          size_t n_visibilities, size_t n_solutions,
                          cu::DeviceMemory& antenna_pairs,
                          cu::DeviceMemory& solution_map,
                          cu::DeviceMemory& solutions, cu::DeviceMemory& model,
                          cu::DeviceMemory& residual);

void LaunchSolveNextSolutionKernel(
    cudaStream_t stream, size_t n_antennas, size_t n_visibilities,
    size_t n_direction_solutions, size_t n_solutions, size_t direction,
    cu::DeviceMemory& antenna_pairs, cu::DeviceMemory& solution_map,
    cu::DeviceMemory& next_solutions, cu::DeviceMemory& numerator,
    cu::DeviceMemory& denominator);

void LaunchSolveDirectionKernel(
    cudaStream_t stream, size_t n_visibilities, size_t n_direction_solutions,
    size_t n_solutions, size_t direction, cu::DeviceMemory& antenna_pairs,
    cu::DeviceMemory& solution_map, cu::DeviceMemory& solutions,
    cu::DeviceMemory& model, cu::DeviceMemory& residual_in,
    cu::DeviceMemory& residual_temp, cu::DeviceMemory& numerator,
    cu::DeviceMemory& denominator);

void LaunchStepKernel(cudaStream_t stream, size_t n_visibilities,
                      cu::DeviceMemory& solutions,
                      cu::DeviceMemory& next_solutions, bool phase_only,
                      double step_size);

#endif  // DP3_DDECAL_GAIN_SOLVERS_KERNELS_ITERATIVEDIAGONAL_H_