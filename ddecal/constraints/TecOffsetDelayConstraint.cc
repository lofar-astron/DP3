#include "TecOffsetDelayConstraint.h"

#include <schaapcommon/threading/dynamicfor.h>
#include <schaapcommon/threading/staticfor.h>

#include <xtensor/core/xmath.hpp>
#include <xtensor/views/xview.hpp>

#include "TecOffsetDelayFitting.h"

namespace dp3::ddecal {

TecOffsetDelayConstraint::TecOffsetDelayConstraint(bool include_offset,
                                                   size_t max_wraps,
                                                   bool do_phase_referencing)
    : include_offset_(include_offset),
      do_phase_referencing_(do_phase_referencing),
      max_wraps_(max_wraps) {}

void TecOffsetDelayConstraint::SetWeights(const std::vector<double>& weights) {
  weights_ = weights;
}

void TecOffsetDelayConstraint::Initialize(
    size_t n_antennas, const std::vector<uint32_t>& solutions_per_direction,
    const std::vector<double>& frequencies) {
  Constraint::Initialize(n_antennas, solutions_per_direction, frequencies);

  assert(!frequencies.empty());
  // We reference the frequency to the first channel to prevent very high values
  reference_frequency_ = frequencies.front();
  relative_frequencies_.resize(frequencies.size());
  for (size_t i = 0; i != frequencies.size(); ++i) {
    relative_frequencies_[i] = frequencies[i] / reference_frequency_;
  }

  const size_t n_threads = schaapcommon::ThreadPool::GetInstance().NThreads();
  thread_data_.resize(n_threads);
  for (ThreadData& local_data : thread_data_) {
    local_data.data.resize(NChannelBlocks());
    local_data.weights.resize(NChannelBlocks());
  }

  weights_.assign(NChannelBlocks() * NAntennas(), 1.0);

  results_.clear();
  ConstraintResult& tec_result = results_.emplace_back();
  tec_result.vals.resize(NAntennas() * NSubSolutions());
  tec_result.weights.resize(NAntennas() * NSubSolutions());
  tec_result.axes = "ant,dir,freq";
  tec_result.name = "tec";
  tec_result.dims.resize(3);
  tec_result.dims[0] = NAntennas();
  tec_result.dims[1] = NSubSolutions();
  tec_result.dims[2] = 1;
  if (include_offset_) {
    ConstraintResult& offset_result = results_.emplace_back(results_.back());
    offset_result.name = "phase";
  }
  ConstraintResult& delay_result = results_.emplace_back(results_.back());
  delay_result.name = "delay";
}

void TecOffsetDelayConstraint::Apply(SolutionSpan& solutions, double time) {
  assert(solutions.shape(0) == NChannelBlocks());
  assert(solutions.shape(1) == NAntennas());
  assert(solutions.shape(2) == NSubSolutions());
  assert(solutions.shape(3) ==
         1);  // single polarization (phase only) solutions

  // Divide out the reference antenna
  if (do_phase_referencing_) ApplyReferenceAntenna(solutions);

  // Use a DynamicFor, since the ternarySearch functions in PhaseFitter.cc run
  // a variable number of iterations.
  schaapcommon::DynamicFor<size_t> loop;
  loop.Run(0, NAntennas() * NSubSolutions(),
           [&](size_t antenna_and_solution_index, size_t thread) {
             const size_t antenna_index =
                 antenna_and_solution_index / NSubSolutions();
             const size_t sub_solution_index =
                 antenna_and_solution_index % NSubSolutions();
             ThreadData& local_data = thread_data_[thread];

             // Flag channels where calibration yielded inf or nan
             double weight_sum = 0.0;
             for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
               const std::complex<double>& solution =
                   solutions(ch, antenna_index, sub_solution_index, 0);
               if (isfinite(solution)) {
                 local_data.data[ch] = std::arg(solution);
                 local_data.weights[ch] =
                     weights_[antenna_index * NChannelBlocks() + ch];
                 weight_sum += weights_[antenna_index * NChannelBlocks() + ch];
               } else {
                 local_data.data[ch] = 0.0;
                 local_data.weights[ch] = 0.0;
               }
             }

             TecOffsetDelayValues values = TecOffsetDelayGridSearch(
                 relative_frequencies_, local_data.data, local_data.weights,
                 include_offset_, max_wraps_);
             results_.back().weights[antenna_and_solution_index] = weight_sum;

             results_[0].vals[antenna_and_solution_index] =
                 values.a * reference_frequency_ / -8.44797245e9;
             results_[0].weights[antenna_and_solution_index] = weight_sum;
             if (include_offset_) {
               results_[1].vals[antenna_and_solution_index] = values.b;
               results_[1].weights[antenna_and_solution_index] = weight_sum;
             }
             results_.back().vals[antenna_and_solution_index] =
                 values.c / (reference_frequency_ * 2.0 * M_PI);
             results_.back().weights[antenna_and_solution_index] = weight_sum;

             EvaluateLinearTecOffsetValues(values, relative_frequencies_,
                                           local_data.data);

             for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
               solutions(ch, antenna_index, sub_solution_index, 0) =
                   std::polar<double>(1.0, local_data.data[ch]);
             }
           });
}

}  // namespace dp3::ddecal
