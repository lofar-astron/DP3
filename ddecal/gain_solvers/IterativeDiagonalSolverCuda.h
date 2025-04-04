// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_GAIN_SOLVERS_ITERATIVE_DIAGONAL_SOLVER_CUDA_H_
#define DDECAL_GAIN_SOLVERS_ITERATIVE_DIAGONAL_SOLVER_CUDA_H_

#include <vector>

#include <cudawrappers/cu.hpp>

#include "IterativeDiagonalSolver.h"
#include "SolverBase.h"
#include "SolveData.h"
#include "../../common/Timer.h"

namespace dp3 {
namespace ddecal {

template <typename VisMatrix>
class IterativeDiagonalSolverCuda final : public SolverBase {
 public:
  IterativeDiagonalSolverCuda(bool keep_buffers = false);
  SolveResult Solve(const SolveData<VisMatrix>& data,
                    std::vector<std::vector<DComplex>>& solutions, double time,
                    std::ostream* stat_stream) override;

  size_t NSolutionPolarizations() const override { return 2; }

  bool SupportsDdSolutionIntervals() const override { return true; }

 private:
  void AllocateGPUBuffers(const SolveData<VisMatrix>& data);
  void DeallocateHostBuffers();
  void AllocateHostBuffers(const SolveData<VisMatrix>& data);

  void CopyHostToHost(size_t ch_block, bool first_iteration,
                      const SolveData<VisMatrix>& data,
                      const std::vector<DComplex>& solutions,
                      cu::Stream& stream);

  void CopyHostToDevice(size_t ch_block, size_t buffer_id, cu::Stream& stream,
                        cu::Event& event, const SolveData<VisMatrix>& data);

  void PostProcessing(size_t& iteration, double time,
                      bool has_previously_converged, bool& has_converged,
                      bool& constraints_satisfied, bool& done,
                      SolverBase::SolveResult& result,
                      std::vector<std::vector<DComplex>>& solutions,
                      SolutionSpan& next_solutions,
                      std::vector<double>& step_magnitudes,
                      std::ostream* stat_stream);

  /// If this variable if false gpu buffers are not initialized
  bool gpu_buffers_initialized_ = false;
  /// if this variable is false host buffers are not initialized
  bool host_buffers_initialized_ = false;

  /// If this variable is false the host buffers are deallocated
  /// at the end of each Solve call
  bool keep_buffers_ = false;

  std::unique_ptr<cu::Device> device_;
  std::unique_ptr<cu::Context> context_;
  std::unique_ptr<cu::Stream> execute_stream_;
  std::unique_ptr<cu::Stream> host_to_device_stream_;
  std::unique_ptr<cu::Stream> device_to_host_stream_;

  /**
   * GPUBuffers hold the GPU memory used in ::Solve()
   *
   * The GPU memory is of type cu::DeviceMemory. This is a wrapper around a
   * plain CUdeviceptr, provided by the cudawrappers library.
   *
   * To facilitate double-buffering, most of the memory is allocated twice (in
   * ::AllocateGPUBuffers) and stored in a vector. There are three exceptions:
   *  - residual: three buffers are used
   *  - numerator: a single buffer is used
   *  - denominator: a single buffer is used
   *
   *
   * For each element in the struct, the comment
   * "<x>[a][b], y" denotes:
   *   - x: the number of elements in the vector (if used)
   *   - a: length of the first dimension
   *   - b: length of the second dimension (if any)
   *   - y: data type
   *
   * For instance for antenna_pairs, <2> denotes
   * that the vector has two elements. [n_antennas][2] denotes that every
   * element is a 2D array with dimensions (n_antennas, 2). Finally, uint32_t
   * denotes the data type.
   */
  struct GPUBuffers {
    // <2>[n_antennas][2], uint32_t
    std::vector<cu::DeviceMemory> antenna_pairs;
    // <2>[n_directions][n_visibilities], uint32_t
    std::vector<cu::DeviceMemory> solution_map;
    // <2>[n_visibilities], DComplex
    std::vector<cu::DeviceMemory> solutions;
    // <2>[n_visibilities], DComplex
    std::vector<cu::DeviceMemory> next_solutions;
    // <2>[n_directions][n_visibilities], MC2x2F
    std::vector<cu::DeviceMemory> model;
    // <3>[n_visibilities], MC2x2F
    std::vector<cu::DeviceMemory> residual;
    // [n_antennas][n_directions], MC2x2FDiag
    std::unique_ptr<cu::DeviceMemory> numerator;
    // [n_antennas][n_directions_solutions], float
    std::unique_ptr<cu::DeviceMemory> denominator;
  } gpu_buffers_;

  /**
   * HostBuffers hold the host memory used in ::Solve()
   *
   * The host memory is of type cu::HostMemory. This is a wrapper around a
   * plain void*, provided by the cudawrappers library.
   *
   * These buffers contain a copy of data that is elsewhere in host memory.
   * Using an extra host-to-cuda-host-memory copy and then a
   * cuda-host-memory-to-gpu copy is faster than a direct host-to-gpu copy.
   */
  struct HostBuffers {
    // <n_channelblocks>[n_directions][n_visibilities], MC2x2F
    std::vector<cu::HostMemory> model;
    // <n_channelblocks>[n_visibilities], MC2x2F
    std::vector<cu::HostMemory> residual;
    // <n_channelblocks>[n_visibilities], DComplex
    std::vector<cu::HostMemory> solutions;
    // [n_channelblocks][n_antennas][n_polarizations], DComplex
    std::unique_ptr<cu::HostMemory> next_solutions;
    // <n_channelblocks>[n_visibilities], std::pair<uin32_t, uint32_t>
    std::vector<cu::HostMemory> antenna_pairs;
    // <n_channelblocks>[n_directions][n_visibilities], uint32_t
    std::vector<cu::HostMemory> solution_map;
  } host_buffers_;
};

}  // namespace ddecal
}  // namespace dp3

#endif  // DDECAL_GAIN_SOLVERS_ITERATIVE_DIAGONAL_SOLVER_CUDA_H_
