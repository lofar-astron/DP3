// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#ifndef MULTI_DIR_BUFFER_H
#define MULTI_DIR_BUFFER_H

#include <complex>
#include <vector>

class MultiDirBuffer
{
public:
  typedef std::complex<float> Complex;
  
  MultiDirBuffer() :
    _nDirections(0), _nChannels(0), _nBaselines(0)
  { }
  
  MultiDirBuffer(size_t nDirections, size_t nChannels, size_t nBaselines) :
    _nDirections(nDirections), _nChannels(nChannels), _nBaselines(nBaselines)
  { }
  
  void SetDimensions(size_t nDirections, size_t nChannels, size_t nBaselines)
  {
    _nDirections = nDirections;
    _nChannels = nChannels;
    _nBaselines = nBaselines;
  }
  
  void CopyAndWeight(const std::vector<Complex *>& dataNoW, const std::vector<float*>& weights, const std::vector<std::vector<Complex *> >& modelDataNoW)
  {
    const size_t nTimes = dataNoW.size();
    _data.resize(nTimes);
    _modelData.resize(nTimes);
    
    for(size_t timestep=0; timestep!=nTimes; ++timestep)
    {
      _data[timestep].resize(_nBaselines*_nChannels*4);
      _modelData[timestep].resize(_nDirections);
      for (size_t dir=0; dir<modelDataNoW[0].size(); ++dir)
        _modelData[timestep][dir].resize(_nBaselines*_nChannels*4);
      for (size_t bl=0; bl<_nBaselines; ++bl)
      {
        for(size_t ch=0; ch!=_nChannels; ++ch)
        {
          bool isFlagged = false;
          for (size_t cr=0; cr<4; ++cr)
          {
            const size_t index = (bl*_nChannels+ch)*4+cr;
            
            if(!isfinite(dataNoW[timestep][index]))
              isFlagged = true;
            else {
              float wSqrt = sqrt(weights[timestep][index]);
              _data[timestep][index] = dataNoW[timestep][index] * wSqrt;
              for (size_t dir=0; dir<_nDirections; ++dir)
              {
                if(!isfinite(modelDataNoW[timestep][dir][index]))
                  isFlagged = true;
                _modelData[timestep][dir][index] = modelDataNoW[timestep][dir][index] * wSqrt;
              }
            }
          }
          
          if(isFlagged)
          {
            for (size_t cr=0; cr<4; ++cr)
            {
              const size_t index = (bl*_nChannels+ch)*4+cr;
              _data[timestep][index] = 0.0;
              for (size_t dir=0; dir<_nDirections; ++dir)
                _modelData[timestep][dir][index] = 0.0;
            }
          }
        }
      }
    }
  }
  
  const std::vector<std::vector<Complex>>& Data() const
  {
    return _data;
  }

  const std::vector<std::vector<std::vector<Complex>>>& ModelData() const
  {
    return _modelData;
  }
private:
  size_t _nDirections, _nChannels, _nBaselines;
  
  std::vector<std::vector<Complex>> _data;
  std::vector<std::vector<std::vector<Complex>>> _modelData;
  
  static bool isfinite(Complex c)
  {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }  
};
  
#endif
