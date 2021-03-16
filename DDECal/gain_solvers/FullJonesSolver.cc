// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "FullJonesSolver.h"

#include "../linear_solvers/QRSolver.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>

using aocommon::ParallelFor;

#include <iomanip>
#include <iostream>

namespace DP3 {
namespace DPPP {

SolverBase::SolveResult FullJonesSolver::Solve(
    const std::vector<Complex*>& unweighted_data,
    const std::vector<float*>& weights,
    std::vector<std::vector<Complex*> >&& unweighted_model_data,
    std::vector<std::vector<DComplex> >& solutions, double time,
    std::ostream* stat_stream) {
  // This algorithm is basically the same as the scalar algorithm,
  // but visibility values are extended to 2x2 matrices and concatenated
  // in the matrix equations as block matrices. One difference is that
  // order of operations are important because of the non-commutativity of
  // matrix multiplication, as well as that A^H B^H = [BA]^H.
  //
  // The approach:
  // First we pre-apply the left-hand solutions to the model to make JM. Each
  // 2x2 coherence matrix Ji is matrix-multied by the lh solutions, for all
  // directions, and visibilities (times x channels).
  //   JMi = Ji Mi
  // These are stacked in matrix JM :
  //        JM0_d0 JM1_d0 ...
  //   JM = JM0_d1 JM1_d1
  //        ...
  // such that JM is a (2D) rows x (2N) col matrix, N=nvis, D=ndir.
  // The solved 2D x 2 solution matrix is similarly formed with the solution
  // values:
  //       ( J0 )
  //   J = ( J1 )
  //       ( .. )
  // And the 2N x 2 visibility matrix as well:
  //       ( V0 )
  //   V = ( V1 )
  //       ( .. )
  // And we solve the equation:
  //   'JM' J^H = V
  // With dimensions:
  //   [ 2N x 2D ] [ 2D x 2 ] = [ 2N x 2 ]

  const size_t n_times = unweighted_data.size();

  buffer_.AssignAndWeight(unweighted_data, weights,
                          std::move(unweighted_model_data));

  for (size_t i = 0; i != constraints_.size(); ++i)
    constraints_[i]->PrepareIteration(false, 0, false);

  std::vector<std::vector<DComplex> > next_solutions(n_channel_blocks_);

  SolveResult result;
#ifndef NDEBUG
  if (solutions.size() != n_channel_blocks_) {
    std::cout << "Error: 'solutions' parameter does not have the right shape\n";
    result.iterations = 0;
    return result;
  }
#endif

  result.results.resize(constraints_.size());

  // Dimensions for each channelblock:
  // Model matrix ant x [2N x 2D] and visibility matrix ant x [2N x 2],
  // The following loop allocates all structures
  std::vector<std::vector<Matrix> > g_times_cs(n_channel_blocks_);
  std::vector<std::vector<Matrix> > vs(n_channel_blocks_);
  for (size_t ch_block = 0; ch_block != n_channel_blocks_; ++ch_block) {
    next_solutions[ch_block].resize(n_directions_ * n_antennas_ * 4);
    const size_t channelIndexStart = ch_block * n_channels_ / n_channel_blocks_,
                 channelIndexEnd =
                     (ch_block + 1) * n_channels_ / n_channel_blocks_,
                 cur_channel_block_size = channelIndexEnd - channelIndexStart;
    g_times_cs[ch_block].resize(n_antennas_);
    vs[ch_block].resize(n_antennas_);

    for (size_t ant = 0; ant != n_antennas_; ++ant) {
      // Model matrix [2N x 2D] and visibility matrix [2N x 2]
      // Space for the auto correlation is also reserved, but they will be set
      // to 0.
      size_t m = n_antennas_ * n_times * cur_channel_block_size * 2;
      size_t n = n_directions_ * 2, n_rhs = 2;
      g_times_cs[ch_block][ant] = Matrix(m, n);
      vs[ch_block][ant] = Matrix(std::max(n, m), n_rhs);
    }
  }

  ///
  /// Start iterating
  ///
  size_t iteration = 0, constrained_iterations = 0;
  bool has_converged = false, has_previously_converged = false,
       constraints_satisfied = false, has_stalled = false;

  std::vector<double> step_magnitudes;
  step_magnitudes.reserve(max_iterations_);

  do {
    MakeSolutionsFinite4Pol(solutions);

    ParallelFor<size_t> loop(n_threads_);
    loop.Run(0, n_channel_blocks_, [&](size_t chBlock, size_t /*thread*/) {
      PerformIteration(chBlock, g_times_cs[chBlock], vs[chBlock],
                       solutions[chBlock], next_solutions[chBlock]);
    });

    Step(solutions, next_solutions);

    if (stat_stream) {
      (*stat_stream) << iteration << '\t';
    }

    constraints_satisfied = true;
    for (size_t i = 0; i != constraints_.size(); ++i) {
      constraints_satisfied =
          constraints_[i]->Satisfied() && constraints_satisfied;
      constraints_[i]->PrepareIteration(has_previously_converged, iteration,
                                        iteration + 1 >= max_iterations_);
      result.results[i] =
          constraints_[i]->Apply(next_solutions, time, stat_stream);
    }

    if (!constraints_satisfied) constrained_iterations = iteration + 1;

    double avgSquaredDiff;
    has_converged =
        AssignSolutions(solutions, next_solutions, !constraints_satisfied,
                        avgSquaredDiff, step_magnitudes, 4);
    if (stat_stream) {
      (*stat_stream) << step_magnitudes.back() << '\t' << avgSquaredDiff
                     << '\n';
    }
    iteration++;

    has_previously_converged = has_converged || has_previously_converged;

    if (detect_stalling_ && constraints_satisfied)
      has_stalled = DetectStall(iteration, step_magnitudes);

  } while (iteration < max_iterations_ &&
           (!has_converged || !constraints_satisfied) && !has_stalled);

  // When we have not converged yet, we set the nr of iterations to the max+1,
  // so that non-converged solves can be distinguished from converged ones.
  if ((!has_converged || !constraints_satisfied) && !has_stalled)
    result.iterations = iteration + 1;
  else
    result.iterations = iteration;
  result.constraint_iterations = constrained_iterations;
  return result;
}

void FullJonesSolver::PerformIteration(size_t channelBlockIndex,
                                       std::vector<Matrix>& g_times_cs,
                                       std::vector<Matrix>& vs,
                                       const std::vector<DComplex>& solutions,
                                       std::vector<DComplex>& next_solutions) {
  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    g_times_cs[ant].SetZero();
    vs[ant].SetZero();
  }

  const size_t channel_index_start =
                   channelBlockIndex * n_channels_ / n_channel_blocks_,
               channel_index_end =
                   (channelBlockIndex + 1) * n_channels_ / n_channel_blocks_,
               cur_channel_block_size = channel_index_end - channel_index_start,
               n_times = buffer_.Data().size();

  // The following loop fills the matrices for all antennas
  for (size_t timeIndex = 0; timeIndex != n_times; ++timeIndex) {
    std::vector<const Complex*> model_ptrs(n_directions_);
    for (size_t baseline = 0; baseline != ant1_.size(); ++baseline) {
      size_t antenna1 = ant1_[baseline];
      size_t antenna2 = ant2_[baseline];
      if (antenna1 != antenna2) {
        // This equation is solved:
        //   J_1 M J_2^H = V
        // for visibilities of the 'normal' correlation ant1 x ant2^H.
        // Since in this equation antenna2 is solved, the solve matrices are
        // called gTimesC2 and v2. The index into these matrices is depending
        // on antenna1, hence the index is called dataIndex1.
        //
        // To use visibilities of correlation ant2 x ant1^H to solve ant1, an
        // extra Herm transpose on M and V is required. The equation is:
        //   J_2 M^H J_1^H = V^H,
        // and the relevant matrices/index are called gTimesC1, v1 and
        // dataIndex2.
        Matrix& g_times_c1 = g_times_cs[antenna1];
        Matrix& g_times_c2 = g_times_cs[antenna2];
        Matrix& v1 = vs[antenna1];
        Matrix& v2 = vs[antenna2];
        for (size_t d = 0; d != n_directions_; ++d)
          model_ptrs[d] =
              &buffer_.ModelData()[timeIndex][d][(channel_index_start +
                                                  baseline * n_channels_) *
                                                 4];
        const Complex* dataPtr =
            &buffer_.Data()[timeIndex]
                           [(channel_index_start + baseline * n_channels_) * 4];
        for (size_t ch = channel_index_start; ch != channel_index_end; ++ch) {
          const size_t data_index1 = 2 * (ch - channel_index_start +
                                          (timeIndex + antenna1 * n_times) *
                                              cur_channel_block_size),
                       data_index2 = 2 * (ch - channel_index_start +
                                          (timeIndex + antenna2 * n_times) *
                                              cur_channel_block_size);

          for (size_t d = 0; d != n_directions_; ++d) {
            aocommon::MC2x2 modelMat(model_ptrs[d]), gTimesC1Mat, gTimesC2Mat;
            size_t solIndex1 = (antenna1 * n_directions_ + d) * 4;
            size_t solIndex2 = (antenna2 * n_directions_ + d) * 4;
            aocommon::Matrix2x2::ATimesB(
                gTimesC2Mat.Data(), &solutions[solIndex1], modelMat.Data());
            aocommon::Matrix2x2::ATimesHermB(
                gTimesC1Mat.Data(), &solutions[solIndex2], modelMat.Data());
            for (size_t p = 0; p != 4; ++p) {
              g_times_c2(data_index1 + (p / 2), d * 2 + p % 2) = gTimesC2Mat[p];
              g_times_c1(data_index2 + (p / 2), d * 2 + p % 2) = gTimesC1Mat[p];
            }

            model_ptrs[d] += 4;  // Goto the next 2x2 matrix.
          }
          for (size_t p = 0; p != 4; ++p) {
            v1(data_index2 + (p % 2), p / 2) = std::conj(*dataPtr);
            v2(data_index1 + (p / 2), p % 2) =
                *dataPtr;  // note that this also performs the Herm transpose
            ++dataPtr;     // Goto the next element of the 2x2 matrix.
          }
        }
      }
    }
  }

  // The matrices have been filled; compute the linear solution
  // for each antenna.

  size_t m = n_antennas_ * n_times * cur_channel_block_size * 2;
  size_t n = n_directions_ * 2, nrhs = 2;
  QRSolver solver(m, n, nrhs);

  for (size_t ant = 0; ant != n_antennas_; ++ant) {
    // solve x^H in [g C] x^H  = v
    bool success = solver.Solve(g_times_cs[ant].data(), vs[ant].data());
    Matrix& x = vs[ant];
    if (success && x(0, 0) != Complex(0.0, 0.0)) {
      for (size_t d = 0; d != n_directions_; ++d) {
        for (size_t p = 0; p != 4; ++p) {
          // The conj transpose is also performed at this point (note swap of %
          // and /)
          next_solutions[(ant * n_directions_ + d) * 4 + p] =
              std::conj(x(d * 2 + p % 2, p / 2));
        }
      }
    } else {
      for (size_t i = 0; i != n_directions_ * 4; ++i) {
        next_solutions[ant * n_directions_ * 4 + i] =
            std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

}  // namespace DPPP
}  // namespace DP3
