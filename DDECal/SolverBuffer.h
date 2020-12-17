// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DDECAL_SOLVER_BUFFER_H
#define DDECAL_SOLVER_BUFFER_H

#include <complex>
#include <vector>

namespace DP3 {
namespace DPPP {

class SolverBuffer {
 public:
  typedef std::complex<float> Complex;

  SolverBuffer() : _nDirections(0), _nChannels(0), _nBaselines(0) {}

  SolverBuffer(size_t nDirections, size_t nChannels, size_t nBaselines)
      : _nDirections(nDirections),
        _nChannels(nChannels),
        _nBaselines(nBaselines) {}

  void SetDimensions(size_t nDirections, size_t nChannels, size_t nBaselines) {
    _nDirections = nDirections;
    _nChannels = nChannels;
    _nBaselines = nBaselines;
  }

  void CopyAndWeight(const std::vector<Complex*>& dataNoW,
                     const std::vector<float*>& weights,
                     const std::vector<std::vector<Complex*>>& modelDataNoW) {
    const size_t nTimes = dataNoW.size();
    _data.resize(nTimes);
    _modelData.resize(nTimes);

    for (size_t timestep = 0; timestep != nTimes; ++timestep) {
      _data[timestep].resize(_nBaselines * _nChannels * 4);
      _modelData[timestep].resize(_nDirections);
      for (size_t dir = 0; dir < modelDataNoW[0].size(); ++dir)
        _modelData[timestep][dir].resize(_nBaselines * _nChannels * 4);
      for (size_t bl = 0; bl < _nBaselines; ++bl) {
        for (size_t ch = 0; ch != _nChannels; ++ch) {
          bool isFlagged = false;
          for (size_t cr = 0; cr < 4; ++cr) {
            const size_t index = (bl * _nChannels + ch) * 4 + cr;

            if (!isfinite(dataNoW[timestep][index]))
              isFlagged = true;
            else {
              float wSqrt = sqrt(weights[timestep][index]);
              _data[timestep][index] = dataNoW[timestep][index] * wSqrt;
              for (size_t dir = 0; dir < _nDirections; ++dir) {
                if (!isfinite(modelDataNoW[timestep][dir][index]))
                  isFlagged = true;
                _modelData[timestep][dir][index] =
                    modelDataNoW[timestep][dir][index] * wSqrt;
              }
            }
          }

          if (isFlagged) {
            for (size_t cr = 0; cr < 4; ++cr) {
              const size_t index = (bl * _nChannels + ch) * 4 + cr;
              _data[timestep][index] = 0.0;
              for (size_t dir = 0; dir < _nDirections; ++dir)
                _modelData[timestep][dir][index] = 0.0;
            }
          }
        }
      }
    }
  }

  const std::vector<std::vector<Complex>>& Data() const { return _data; }

  const std::vector<std::vector<std::vector<Complex>>>& ModelData() const {
    return _modelData;
  }

 private:
  size_t _nDirections, _nChannels, _nBaselines;

  std::vector<std::vector<Complex>> _data;
  std::vector<std::vector<std::vector<Complex>>> _modelData;

  static bool isfinite(Complex c) {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DDECAL_SOLVER_BUFFER_H
