// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IterativeDiagonalSolverCuda.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include <cuda_runtime.h>
#include <nvToolsExt.h>

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

#include "kernels/IterativeDiagonal.h"

using aocommon::MC2x2F;
using aocommon::MC2x2FDiag;

namespace {

size_t SizeOfModel(size_t n_directions, size_t n_visibilities) {
  return n_directions * n_visibilities * sizeof(MC2x2F);
}

size_t SizeOfResidual(size_t n_visibilities) {
  return n_visibilities * sizeof(MC2x2F);
}

size_t SizeOfSolutions(size_t n_visibilities) {
  return n_visibilities * sizeof(std::complex<double>);
}

size_t SizeOfAntennaPairs(size_t n_visibilities) {
  return n_visibilities * 2 * sizeof(uint32_t);
}

size_t SizeOfSolutionMap(size_t n_directions, size_t n_visibilities) {
  return n_directions * n_visibilities * sizeof(uint32_t);
}

size_t SizeOfNextSolutions(size_t n_visibilities) {
  return n_visibilities * sizeof(std::complex<double>);
}

size_t SizeOfNumerator(size_t n_antennas, size_t n_direction_solutions) {
  return n_antennas * n_direction_solutions * sizeof(MC2x2FDiag);
}

size_t SizeOfDenominator(size_t n_antennas, size_t n_direction_solutions) {
  return n_antennas * n_direction_solutions * 2 * sizeof(float);
}

void SolveDirection(
    const dp3::ddecal::SolveData::ChannelBlockData& channel_block_data,
    cu::Stream& stream, size_t n_antennas, size_t n_solutions, size_t direction,
    cu::DeviceMemory& device_residual_in,
    cu::DeviceMemory& device_residual_temp,
    cu::DeviceMemory& device_solution_map, cu::DeviceMemory& device_solutions,
    cu::DeviceMemory& device_model, cu::DeviceMemory& device_next_solutions,
    cu::DeviceMemory& device_antenna_pairs, cu::DeviceMemory& device_numerator,
    cu::DeviceMemory& device_denominator) {
  // Calculate this equation, given ant a:
  //
  //          sum_b data_ab * solutions_b * model_ab^*
  // sol_a =  ----------------------------------------
  //             sum_b norm(model_ab * solutions_b)
  const size_t n_direction_solutions =
      channel_block_data.NSolutionsForDirection(direction);
  const size_t n_visibilities = channel_block_data.NVisibilities();

  // Initialize values to 0
  device_numerator.zero(SizeOfNumerator(n_antennas, n_direction_solutions),
                        stream);
  device_denominator.zero(SizeOfDenominator(n_antennas, n_direction_solutions),
                          stream);

  stream.memcpyDtoDAsync(device_residual_temp, device_residual_in,
                         SizeOfResidual(n_visibilities));

  LaunchSolveDirectionKernel(
      stream, n_visibilities, n_direction_solutions, n_solutions, direction,
      device_antenna_pairs, device_solution_map, device_solutions, device_model,
      device_residual_in, device_residual_temp, device_numerator,
      device_denominator);

  LaunchSolveNextSolutionKernel(
      stream, n_antennas, n_visibilities, n_direction_solutions, n_solutions,
      direction, device_antenna_pairs, device_solution_map,
      device_next_solutions, device_numerator, device_denominator);
}

void PerformIteration(
    bool phase_only, double step_size,
    const dp3::ddecal::SolveData::ChannelBlockData& channel_block_data,
    cu::Stream& stream, size_t n_antennas, size_t n_solutions,
    size_t n_directions, cu::DeviceMemory& device_solution_map,
    cu::DeviceMemory& device_solutions, cu::DeviceMemory& device_next_solutions,
    cu::DeviceMemory& device_residual, cu::DeviceMemory& device_residual_temp,
    cu::DeviceMemory& device_model, cu::DeviceMemory& device_antenna_pairs,
    cu::DeviceMemory& device_numerator, cu::DeviceMemory& device_denominator) {
  const size_t n_visibilities = channel_block_data.NVisibilities();

  // Subtract all directions with their current solutions
  // In-place: residual -> residual
  LaunchSubtractKernel(stream, n_directions, n_visibilities, n_solutions,
                       device_antenna_pairs, device_solution_map,
                       device_solutions, device_model, device_residual);

  for (size_t direction = 0; direction != n_directions; direction++) {
    // Be aware that we purposely still use the subtraction with 'old'
    // solutions, because the new solutions have not been constrained yet. Add
    // this direction back before solving

    // Out-of-place: residual -> residual_temp
    SolveDirection(channel_block_data, stream, n_antennas, n_solutions,
                   direction, device_residual, device_residual_temp,
                   device_solution_map, device_solutions, device_model,
                   device_next_solutions, device_antenna_pairs,
                   device_numerator, device_denominator);
  }

  LaunchStepKernel(stream, n_visibilities, device_solutions,
                   device_next_solutions, phase_only, step_size);
}

std::tuple<size_t, size_t, size_t> ComputeArrayDimensions(
    const dp3::ddecal::SolveData& data) {
  size_t max_n_direction_solutions = 0;
  size_t max_n_visibilities = 0;
  size_t max_n_directions = 0;

  for (size_t ch_block = 0; ch_block < data.NChannelBlocks(); ch_block++) {
    const dp3::ddecal::SolveData::ChannelBlockData& channel_block_data =
        data.ChannelBlock(ch_block);
    max_n_visibilities =
        std::max(max_n_visibilities, channel_block_data.NVisibilities());
    max_n_directions =
        std::max(max_n_directions, channel_block_data.NDirections());
    for (size_t direction = 0; direction < channel_block_data.NDirections();
         direction++) {
      max_n_direction_solutions =
          std::max(max_n_direction_solutions,
                   static_cast<size_t>(
                       channel_block_data.NSolutionsForDirection(direction)));
    }
  }

  return std::make_tuple(max_n_direction_solutions, max_n_visibilities,
                         max_n_directions);
}
}  // namespace

namespace dp3 {
namespace ddecal {

IterativeDiagonalSolverCuda::IterativeDiagonalSolverCuda(bool keep_buffers)
    : keep_buffers_{keep_buffers} {
  cu::init();
  device_ = std::make_unique<cu::Device>(0);
  context_ = std::make_unique<cu::Context>(0, *device_);
  context_->setCurrent();
  execute_stream_ = std::make_unique<cu::Stream>();
  host_to_device_stream_ = std::make_unique<cu::Stream>();
  device_to_host_stream_ = std::make_unique<cu::Stream>();
}

void IterativeDiagonalSolverCuda::AllocateGPUBuffers(const SolveData& data) {
  size_t max_n_direction_solutions = 0;
  size_t max_n_visibilities = 0;
  size_t max_n_directions = 0;
  std::tie(max_n_direction_solutions, max_n_visibilities, max_n_directions) =
      ComputeArrayDimensions(data);

  gpu_buffers_.numerator = std::make_unique<cu::DeviceMemory>(
      SizeOfNumerator(NAntennas(), max_n_direction_solutions));
  gpu_buffers_.denominator = std::make_unique<cu::DeviceMemory>(
      SizeOfDenominator(NAntennas(), max_n_direction_solutions));
  // Allocating two buffers allows double buffering.
  for (size_t i = 0; i < 2; i++) {
    gpu_buffers_.antenna_pairs.emplace_back(
        SizeOfAntennaPairs(max_n_visibilities));
    gpu_buffers_.solution_map.emplace_back(
        SizeOfSolutionMap(max_n_directions, max_n_visibilities));
    gpu_buffers_.solutions.emplace_back(SizeOfSolutions(max_n_visibilities));
    gpu_buffers_.next_solutions.emplace_back(
        SizeOfNextSolutions(max_n_visibilities));
    gpu_buffers_.model.emplace_back(
        SizeOfModel(max_n_directions, max_n_visibilities));
  }

  // We need two buffers for residual like above to facilitate double-buffering,
  // the third buffer is used for the per-direction add/subtract.
  for (size_t i = 0; i < 3; i++) {
    gpu_buffers_.residual.emplace_back(SizeOfResidual(max_n_visibilities));
  }
}

void IterativeDiagonalSolverCuda::AllocateHostBuffers(const SolveData& data) {
  host_buffers_.next_solutions =
      std::make_unique<cu::HostMemory>(SizeOfNextSolutions(NVisibilities()));
  for (size_t ch_block = 0; ch_block < NChannelBlocks(); ch_block++) {
    const SolveData::ChannelBlockData& channel_block_data =
        data.ChannelBlock(ch_block);
    const size_t n_directions = channel_block_data.NDirections();
    const size_t n_visibilities = channel_block_data.NVisibilities();
    host_buffers_.model.emplace_back(SizeOfModel(n_directions, n_visibilities));
    host_buffers_.residual.emplace_back(SizeOfResidual(n_visibilities));
    host_buffers_.solutions.emplace_back(SizeOfSolutions(n_visibilities));
    host_buffers_.antenna_pairs.emplace_back(
        SizeOfAntennaPairs(n_visibilities));
    host_buffers_.solution_map.emplace_back(
        SizeOfSolutionMap(n_directions, n_visibilities));
    uint32_t* antenna_pairs =
        static_cast<uint32_t*>(host_buffers_.antenna_pairs[ch_block]);
    for (size_t visibility_index = 0; visibility_index < n_visibilities;
         visibility_index++) {
      antenna_pairs[visibility_index * 2 + 0] =
          channel_block_data.Antenna1Index(visibility_index);
      antenna_pairs[visibility_index * 2 + 1] =
          channel_block_data.Antenna2Index(visibility_index);
    }
  }
}
void IterativeDiagonalSolverCuda::DeallocateHostBuffers() {
  host_buffers_.next_solutions.reset();
  host_buffers_.model.clear();
  host_buffers_.residual.clear();
  host_buffers_.solutions.clear();
  host_buffers_.antenna_pairs.clear();
  host_buffers_.solution_map.clear();
  host_buffers_initialized_ = false;
}

void IterativeDiagonalSolverCuda::CopyHostToHost(
    size_t ch_block, bool first_iteration, const SolveData& data,
    const std::vector<std::complex<double>>& solutions, cu::Stream& stream) {
  const SolveData::ChannelBlockData& channel_block_data =
      data.ChannelBlock(ch_block);
  const size_t n_directions = channel_block_data.NDirections();
  const size_t n_visibilities = channel_block_data.NVisibilities();
  cu::HostMemory& host_model = host_buffers_.model[ch_block];
  cu::HostMemory& host_solutions = host_buffers_.solutions[ch_block];
  stream.memcpyHtoHAsync(host_model, &channel_block_data.ModelVisibility(0, 0),
                         SizeOfModel(n_directions, n_visibilities));
  stream.memcpyHtoHAsync(host_solutions, solutions.data(),
                         SizeOfSolutions(n_visibilities));
  if (first_iteration) {
    cu::HostMemory& host_residual = host_buffers_.residual[ch_block];
    cu::HostMemory& host_solution_map = host_buffers_.solution_map[ch_block];
    stream.memcpyHtoHAsync(host_residual, &channel_block_data.Visibility(0),
                           SizeOfResidual(n_visibilities));
    stream.memcpyHtoHAsync(host_solution_map,
                           channel_block_data.SolutionMapData(),
                           SizeOfSolutionMap(n_directions, n_visibilities));
  }
}

void IterativeDiagonalSolverCuda::CopyHostToDevice(size_t ch_block,
                                                   size_t buffer_id,
                                                   cu::Stream& stream,
                                                   cu::Event& event,
                                                   const SolveData& data) {
  const dp3::ddecal::SolveData::ChannelBlockData& channel_block_data =
      data.ChannelBlock(ch_block);

  const size_t n_directions = channel_block_data.NDirections();
  const size_t n_visibilities = channel_block_data.NVisibilities();

  cu::HostMemory& host_solution_map = host_buffers_.solution_map[ch_block];
  cu::HostMemory& host_antenna_pairs = host_buffers_.antenna_pairs[ch_block];
  cu::HostMemory& host_model = host_buffers_.model[ch_block];
  cu::HostMemory& host_residual = host_buffers_.residual[ch_block];
  cu::HostMemory& host_solutions = host_buffers_.solutions[ch_block];
  cu::DeviceMemory& device_solution_map = gpu_buffers_.solution_map[buffer_id];
  cu::DeviceMemory& device_antenna_pairs =
      gpu_buffers_.antenna_pairs[buffer_id];
  cu::DeviceMemory& device_model = gpu_buffers_.model[buffer_id];
  cu::DeviceMemory& device_residual = gpu_buffers_.residual[buffer_id];
  cu::DeviceMemory& device_solutions = gpu_buffers_.solutions[buffer_id];

  stream.memcpyHtoDAsync(device_solution_map, host_solution_map,
                         SizeOfSolutionMap(n_directions, n_visibilities));
  stream.memcpyHtoDAsync(device_model, host_model,
                         SizeOfModel(n_directions, n_visibilities));
  stream.memcpyHtoDAsync(device_residual, host_residual,
                         SizeOfResidual(n_visibilities));
  stream.memcpyHtoDAsync(device_antenna_pairs, host_antenna_pairs,
                         SizeOfAntennaPairs(n_visibilities));
  stream.memcpyHtoDAsync(device_solutions, host_solutions,
                         SizeOfSolutions(n_visibilities));

  stream.record(event);
}

void IterativeDiagonalSolverCuda::PostProcessing(
    size_t& iteration, double time, bool has_previously_converged,
    bool& has_converged, bool& constraints_satisfied, bool& done,
    SolverBase::SolveResult& result,
    std::vector<std::vector<std::complex<double>>>& solutions,
    SolutionSpan& next_solutions, std::vector<double>& step_magnitudes,
    std::ostream* stat_stream) {
  constraints_satisfied =
      ApplyConstraints(iteration, time, has_previously_converged, result,
                       next_solutions, stat_stream);

  double avg_squared_diff;
  has_converged =
      AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                      avg_squared_diff, step_magnitudes);
  iteration++;

  has_previously_converged = has_converged || has_previously_converged;

  done = ReachedStoppingCriterion(iteration, has_converged,
                                  constraints_satisfied, step_magnitudes);
}

IterativeDiagonalSolver::SolveResult IterativeDiagonalSolverCuda::Solve(
    const SolveData& data,
    std::vector<std::vector<std::complex<double>>>& solutions, double time,
    std::ostream* stat_stream) {
  PrepareConstraints();
  context_->setCurrent();

  const bool phase_only = GetPhaseOnly();
  const double step_size = GetStepSize();

  SolveResult result;

  /*
   * Allocate buffers
   */
  if (!host_buffers_initialized_) {
    AllocateHostBuffers(data);
    host_buffers_initialized_ = true;
  }
  if (!gpu_buffers_initialized_) {
    AllocateGPUBuffers(data);
    gpu_buffers_initialized_ = true;
  }

  const std::array<size_t, 4> next_solutions_shape = {
      NChannelBlocks(), NAntennas(), NSolutions(), NSolutionPolarizations()};
  std::complex<double>* next_solutions_ptr = *(host_buffers_.next_solutions);
  SolutionSpan next_solutions =
      aocommon::xt::CreateSpan(next_solutions_ptr, next_solutions_shape);

  /*
   * Allocate events for each channel block
   */
  std::vector<cu::Event> input_copied_events(NChannelBlocks());
  std::vector<cu::Event> compute_finished_events(NChannelBlocks());
  std::vector<cu::Event> output_copied_events(NChannelBlocks());

  /*
   * Start iterating
   */
  size_t iteration = 0;
  bool has_converged = false;
  bool has_previously_converged = false;
  bool constraints_satisfied = false;
  bool done = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(GetMaxIterations());

  do {
    MakeSolutionsFinite2Pol(solutions);

    nvtxRangeId_t nvts_range_gpu = nvtxRangeStart("GPU");

    for (size_t ch_block = 0; ch_block < NChannelBlocks(); ch_block++) {
      const SolveData::ChannelBlockData& channel_block_data =
          data.ChannelBlock(ch_block);
      const int buffer_id = ch_block % 2;
      // Copy input data for first channel block
      if (ch_block == 0) {
        CopyHostToHost(ch_block, iteration == 0, data, solutions[ch_block],
                       *host_to_device_stream_);

        CopyHostToDevice(ch_block, buffer_id, *host_to_device_stream_,
                         input_copied_events[0], data);
      }

      // As soon as input_copied_events[0] is triggered, the input data is
      // copied to the GPU and the host buffers could theoretically be reused.
      // However, since the size of these buffers may differ, every channel
      // block has its own set of host buffers anyway.
      // Before starting kernel execution for the current channel block (on a
      // different stream), the copy of data for the next channel block (if any)
      // is scheduled using a second set of GPU buffers.
      if (ch_block < NChannelBlocks() - 1) {
        CopyHostToHost(ch_block + 1, iteration == 0, data,
                       solutions[ch_block + 1], *host_to_device_stream_);

        // Since the computation of channel block <n> and <n + 2> share the same
        // set of GPU buffers, wait for the compute_finished event to be
        // triggered before overwriting their contents.
        if (ch_block > 1) {
          host_to_device_stream_->wait(compute_finished_events[ch_block - 2]);
        }

        CopyHostToDevice(ch_block + 1, (ch_block + 1) % 2,
                         *host_to_device_stream_,
                         input_copied_events[ch_block + 1], data);
      }

      // Wait for input of the current channel block to be copied
      execute_stream_->wait(input_copied_events[ch_block]);

      // Wait for output buffer to be free
      if (ch_block > 1) {
        execute_stream_->wait(output_copied_events[ch_block - 2]);
      }

      // Start iteration (dtod copies and kernel execution only)
      PerformIteration(phase_only, step_size, channel_block_data,
                       *execute_stream_, NAntennas(), NSolutions(),
                       NDirections(), gpu_buffers_.solution_map[buffer_id],
                       gpu_buffers_.solutions[buffer_id],
                       gpu_buffers_.next_solutions[buffer_id],
                       gpu_buffers_.residual[buffer_id],
                       gpu_buffers_.residual[2], gpu_buffers_.model[buffer_id],
                       gpu_buffers_.antenna_pairs[buffer_id],
                       *gpu_buffers_.numerator, *gpu_buffers_.denominator);

      execute_stream_->record(compute_finished_events[ch_block]);

      // Wait for the computation to finish
      device_to_host_stream_->wait(compute_finished_events[ch_block]);

      // Copy next solutions back to host
      const size_t n_visibilities = next_solutions.shape(1) *
                                    next_solutions.shape(2) *
                                    next_solutions.shape(3);
      device_to_host_stream_->memcpyDtoHAsync(
          &next_solutions(ch_block, 0, 0, 0),
          gpu_buffers_.next_solutions[buffer_id],
          SizeOfNextSolutions(n_visibilities));

      // Record that the output is copied
      device_to_host_stream_->record(output_copied_events[ch_block]);
    }  // end for ch_block

    // Wait for next solutions to be copied
    device_to_host_stream_->synchronize();

    nvtxRangeEnd(nvts_range_gpu);

    // CPU-only postprocessing
    nvtxRangeId_t nvtx_range_cpu = nvtxRangeStart("CPU");
    PostProcessing(iteration, time, has_previously_converged, has_converged,
                   constraints_satisfied, done, result, solutions,
                   next_solutions, step_magnitudes, stat_stream);
    nvtxRangeEnd(nvtx_range_cpu);
  } while (!done);

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged iterations can be distinguished from converged ones.
  if (has_converged && constraints_satisfied) {
    result.iterations = iteration;
  } else {
    result.iterations = iteration + 1;
  }

  if (!keep_buffers_) DeallocateHostBuffers();
  return result;
}

}  // namespace ddecal
}  // namespace dp3
