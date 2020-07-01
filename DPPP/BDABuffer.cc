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

#include "BDABuffer.h"

namespace DP3 {
  namespace DPPP {

    BDABuffer::Row::Row(const double time,
                        const double exposure,
                        const rownr_t rowNr,
                        const std::size_t baselineNr,
                        const std::size_t nChannels,
                        const std::size_t nCorrelations,
                        const std::vector<std::complex<float>>::iterator data,
                        const std::vector<bool>::iterator flags,
                        const std::vector<float>::iterator weights,
                        const std::vector<bool>::iterator fullResFlags,
                        const double *const UVW)
    : itsTime(time)
    , itsExposure(exposure)
    , itsRowNr(rowNr)
    , itsBaselineNr(baselineNr)
    , itsNChannels(nChannels)
    , itsNCorrelations(nCorrelations)
    , itsData(data)
    , itsFlags(flags)
    , itsWeights(weights)
    , itsFullResFlags(fullResFlags)
    , itsUVW{ UVW ? UVW[0] : std::nan(""),
              UVW ? UVW[1] : std::nan(""),
              UVW ? UVW[2] : std::nan("") }
    {}

    BDABuffer::BDABuffer(const std::size_t poolSize)
    : itsData(), itsFlags(), itsWeights(), itsFullResFlags(), itsRows()
    {
      itsData.reserve(poolSize);
      itsFlags.reserve(poolSize);
      itsWeights.reserve(poolSize);
      itsFullResFlags.reserve(poolSize);
    }

    BDABuffer::BDABuffer(const BDABuffer& other)
    : itsData(other.itsData)
    , itsFlags(other.itsFlags)
    , itsWeights(other.itsWeights)
    , itsFullResFlags(other.itsFullResFlags)
    , itsRows()
    {
      // Copy rows but set their iterators to the new memory pools.
      auto dataIt = itsData.begin();
      auto flagsIt = itsFlags.begin();
      auto weightsIt = itsWeights.begin();
      auto fullResFlagsIt = itsFullResFlags.begin();

      itsRows.reserve(other.itsRows.size());
      for (const auto& row : other.itsRows) {
        itsRows.emplace_back(row.itsTime,
                             row.itsExposure,
                             row.itsRowNr,
                             row.itsBaselineNr,
                             row.itsNChannels,
                             row.itsNCorrelations,
                             dataIt,
                             flagsIt,
                             weightsIt,
                             fullResFlagsIt,
                             row.itsUVW);

        const std::size_t dataSize = row.itsNChannels * row.itsNCorrelations;
        dataIt += dataSize;
        flagsIt += dataSize;
        weightsIt += dataSize;
        fullResFlagsIt += dataSize;
      }
    }

    bool BDABuffer::addRow(const double time,
                           const double exposure,
                           const rownr_t rowNr,
                           const std::size_t baselineNr,
                           const std::size_t nChannels,
                           const std::size_t nCorrelations,
                           const std::complex<float>* const data,
                           const bool* const flags,
                           const float* const weights,
                           const bool* const fullResFlags,
                           const double* const UVW)
    {
      const std::size_t dataSize = nChannels * nCorrelations;

      // Check if there is enough capacity left.
      if ((itsData.capacity() - itsData.size()) < dataSize) {
        return false;
      }

      itsRows.emplace_back(time,
                           exposure,
                           rowNr,
                           baselineNr,
                           nChannels,
                           nCorrelations,
                           itsData.end(),
                           itsFlags.end(),
                           itsWeights.end(),
                           itsFullResFlags.end(),
                           UVW);
      if (data) {
        itsData.insert(itsData.end(), data, data + dataSize);
      } else {
        const std::complex<float> nan(std::nanf(""), std::nanf(""));
        itsData.insert(itsData.end(), dataSize, nan);
      }

      if (flags) {
        itsFlags.insert(itsFlags.end(), flags, flags + dataSize);
      } else {
        itsFlags.insert(itsFlags.end(), dataSize, false);
      }

      if (weights) {
        itsWeights.insert(itsWeights.end(), weights, weights + dataSize);
      } else {
        const float nan = std::nanf("");
        itsWeights.insert(itsWeights.end(), dataSize, nan);
      }

      if (fullResFlags) {
        itsFullResFlags.insert(itsFullResFlags.end(),
                               fullResFlags,
                               fullResFlags + dataSize);
      } else {
        itsFullResFlags.insert(itsFullResFlags.end(), dataSize, false);
      }

      return true;
    }
  }
}
