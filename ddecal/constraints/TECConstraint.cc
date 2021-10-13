// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "TECConstraint.h"

#include <aocommon/parallelfor.h>

namespace dp3 {
namespace ddecal {

TECConstraintBase::TECConstraintBase(Mode mode)
    : _mode(mode), _doPhaseReference(true), _phaseFitters() {}

void TECConstraintBase::Initialize(size_t nAntennas, size_t nDirections,
                                   const std::vector<double>& frequencies) {
  Constraint::Initialize(nAntennas, nDirections, frequencies);

  _phaseFitters.resize(NThreads());
  for (PhaseFitter& fitter : _phaseFitters) fitter.Initialize(frequencies);

  _weights.assign(NChannelBlocks() * NAntennas(), 1.0);
  initializeChild();
}

void TECConstraintBase::SetWeights(const std::vector<double>& weights) {
  _weights = weights;
}

void ApproximateTECConstraint::initializeChild() {
  _pwFitters.resize(NThreads());
  _threadData.resize(_pwFitters.size());
  _threadFittedData.resize(_pwFitters.size());
  _threadWeights.resize(_pwFitters.size());
  for (size_t threadId = 0; threadId != _pwFitters.size(); ++threadId) {
    _threadData[threadId].resize(NChannelBlocks());
    _threadFittedData[threadId].resize(NChannelBlocks());
    _threadWeights[threadId].resize(NChannelBlocks());
  }

  if (_fittingChunkSize == 0) {
    _fittingChunkSize = PieceWisePhaseFitter::CalculateChunkSize(
        _phaseFitters.front().GetFrequencies());
  }
  for (size_t i = 0; i != _pwFitters.size(); ++i)
    _pwFitters[i].SetChunkSize(_fittingChunkSize);
}

void TECConstraintBase::applyReferenceAntenna(
    std::vector<std::vector<dcomplex>>& solutions) const {
  // Choose reference antenna that has at least 20% channels unflagged
  size_t refAntenna = 0;
  for (; refAntenna != NAntennas(); ++refAntenna) {
    size_t nUnFlaggedChannels = 0;
    // Only check flagged state for first direction
    for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
      if (isfinite(solutions[ch][refAntenna * NDirections()]))
        nUnFlaggedChannels++;
    }
    if (nUnFlaggedChannels * 1.0 / NChannelBlocks() > 0.2)
      // Choose this refAntenna;
      break;
  }
  // All antennas are flagged, use first one (will lead to NaNs for this solint)
  if (refAntenna == NAntennas()) refAntenna = 0;

  for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
    for (size_t antennaIndex = 0; antennaIndex != NAntennas(); ++antennaIndex) {
      for (size_t d = 0; d != NDirections(); ++d) {
        size_t solutionIndex = antennaIndex * NDirections() + d;
        size_t refAntennaIndex = d + refAntenna * NDirections();
        if (antennaIndex != refAntenna) {
          solutions[ch][solutionIndex] =
              solutions[ch][solutionIndex] / solutions[ch][refAntennaIndex];
        }
      }
    }
    for (size_t d = 0; d != NDirections(); ++d)
      solutions[ch][refAntenna * NDirections() + d] = 1.0;
  }
}

std::vector<Constraint::Result> TECConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double,
    std::ostream* /*statStream*/) {
  size_t nRes = 3;
  if (_mode == TECOnlyMode) {
    nRes = 2;  // TEC and error
  } else {
    nRes = 3;  // TEC, phase and error
  }

  std::vector<Constraint::Result> res(nRes);
  res[0].vals.resize(NAntennas() * NDirections());
  res[0].weights.resize(NAntennas() * NDirections());
  res[0].axes = "ant,dir,freq";
  res[0].name = "tec";
  res[0].dims.resize(3);
  res[0].dims[0] = NAntennas();
  res[0].dims[1] = NDirections();
  res[0].dims[2] = 1;
  if (_mode == TECAndCommonScalarMode) {
    res[1] = res[0];
    res[1].name = "phase";
  }
  res.back() = res[0];
  res.back().name = "error";

  // Divide out the reference antenna
  if (_doPhaseReference) applyReferenceAntenna(solutions);

  aocommon::ParallelFor<size_t> loop(NThreads());
  loop.Run(
      0, NAntennas() * NDirections(), [&](size_t solutionIndex, size_t thread) {
        size_t antennaIndex = solutionIndex / NDirections();

        // Flag channels where calibration yielded inf or nan
        double weightSum = 0.0;
        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          if (isfinite(solutions[ch][solutionIndex])) {
            _phaseFitters[thread].PhaseData()[ch] =
                std::arg(solutions[ch][solutionIndex]);
            _phaseFitters[thread].WeightData()[ch] =
                _weights[antennaIndex * NChannelBlocks() + ch];
            weightSum += _weights[antennaIndex * NChannelBlocks() + ch];
          } else {
            _phaseFitters[thread].PhaseData()[ch] = 0.0;
            _phaseFitters[thread].WeightData()[ch] = 0.0;
          }
        }

        double alpha, beta = 0.0;
        if (_mode == TECOnlyMode) {
          res.back().vals[solutionIndex] =
              _phaseFitters[thread].FitDataToTEC1Model(alpha);
        } else {
          res.back().vals[solutionIndex] =
              _phaseFitters[thread].FitDataToTEC2Model(alpha, beta);
        }
        res.back().weights[solutionIndex] = weightSum;

        res[0].vals[solutionIndex] = alpha / -8.44797245e9;
        res[0].weights[solutionIndex] = weightSum;
        if (_mode == TECAndCommonScalarMode) {
          res[1].vals[solutionIndex] = beta;
          res[1].weights[solutionIndex] = weightSum;
        }

        for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
          solutions[ch][solutionIndex] =
              std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
        }
      });

  return res;
}

std::vector<Constraint::Result> ApproximateTECConstraint::Apply(
    std::vector<std::vector<dcomplex>>& solutions, double time,
    std::ostream* statStream) {
  if (_finishedApproximateStage)
    return TECConstraint::Apply(solutions, time, statStream);
  else {
    if (_doPhaseReference) applyReferenceAntenna(solutions);

    aocommon::ParallelFor<size_t> loop(NThreads());
    loop.Run(0, NAntennas() * NDirections(),
             [&](size_t solutionIndex, size_t thread) {
               size_t antennaIndex = solutionIndex / NDirections();
               std::vector<double>& data = _threadData[thread];
               std::vector<double>& fittedData = _threadFittedData[thread];
               std::vector<double>& weights = _threadWeights[thread];

               // Flag channels where calibration yielded inf or nan
               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 if (isfinite(solutions[ch][solutionIndex])) {
                   data[ch] = std::arg(solutions[ch][solutionIndex]);
                   weights[ch] = _weights[antennaIndex * NChannelBlocks() + ch];
                 } else {
                   data[ch] = 0.0;
                   weights[ch] = 0.0;
                 }
               }

               // TODO might be nice to make it a user option whether to break
               // or not
               _pwFitters[thread].SlidingFitWithBreak(
                   _phaseFitters[thread].GetFrequencies().data(), data.data(),
                   weights.data(), fittedData.data(), data.size());

               for (size_t ch = 0; ch != NChannelBlocks(); ++ch) {
                 solutions[ch][solutionIndex] =
                     std::polar<double>(1.0, fittedData[ch]);
               }
             });

    return std::vector<Constraint::Result>();
  }
}

}  // namespace ddecal
}  // namespace dp3
