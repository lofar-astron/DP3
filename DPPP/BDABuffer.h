// BDABuffer.h: Buffer holding BDA data
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

/// @file
/// @brief Buffer holding base dependent averaged (BDA) data.
/// @author Maik Nijhuis and Lars Krombeen

#ifndef DPPP_BDABUFFER_H
#define DPPP_BDABUFFER_H

#include "../Common/Types.h"

#include <complex>
#include <vector>

namespace DP3 {
  namespace DPPP {

  class BDABuffer
  {
  public:
    struct Row {
      Row(double time,
          double exposure,
          rownr_t rowNr,
          std::size_t nChannels,
          std::size_t nCorrelations,
          std::vector<std::complex<float>>::iterator data,
          std::vector<bool>::iterator flags,
          std::vector<float>::iterator weights,
          std::vector<bool>::iterator fullResFlags,
          const double *const UVW);

      double itsTime; ///< Start time for the measurements in MJD seconds.
      double itsExposure; ///< Exposure duration for the measurements in seconds.
      rownr_t itsRowNr;
      std::size_t itsNChannels;
      std::size_t itsNCorrelations;
      std::vector<std::complex<float>>::iterator itsData;
      std::vector<bool>::iterator itsFlags;
      std::vector<float>::iterator itsWeights;
      std::vector<bool>::iterator itsFullResFlags;
      double itsUVW[3];
    };
    
  public:
    /**
     * @param poolSize Size of the memory pool for this buffer
     *                 (number of complex values)
     */
    explicit BDABuffer(std::size_t poolSize);

    /**
     * Add a measurement line to the buffer.
     * @return True if the line is added.
     *         False if the buffer is full.
     */
    bool addRow(double time,
                double exposure,
                rownr_t rowNr,
                std::size_t nChannels,
                std::size_t nCorrelations,
                const std::complex<float>* data = nullptr,
                const bool* flags = nullptr,
                const float* weights = nullptr,
                const bool* fullResFlags = nullptr,
                const double* UVW = nullptr);

    const std::vector<std::complex<float>>& getData() const {
      return itsData;
    }
    std::vector<std::complex<float>>& getData() {
      return itsData;
    }

    const std::vector<bool>& getFlags() const {
      return itsFlags;
    }
    std::vector<bool>& getFlags() {
      return itsFlags;
    }

    const std::vector<float>& getWeights() const {
      return itsWeights;
    }
    std::vector<float>& getWeights() {
      return itsWeights;
    }

    const std::vector<bool>& getFullResFlags() const {
      return itsFlags;
    }
    std::vector<bool>& getFullResFlags() {
      return itsFlags;
    }

    std::vector<std::complex<float>>::const_iterator getData(const std::size_t row)  const {
      return itsRows[row].itsData;
    }
    std::vector<std::complex<float>>::iterator getData(const std::size_t row) {
      return itsRows[row].itsData;
    }

    std::vector<bool>::const_iterator getFlags(const std::size_t row) const {
      return itsRows[row].itsFlags;
    }
    std::vector<bool>::iterator getFlags(const std::size_t row) {
      return itsRows[row].itsFlags;
    }

    std::vector<float>::const_iterator getWeights(const std::size_t row) const {
      return itsRows[row].itsWeights;
    }
    std::vector<float>::iterator getWeights(const std::size_t row) {
      return itsRows[row].itsWeights;
    }

    std::vector<bool>::const_iterator getFullResFlags(const std::size_t row) const {
      return itsRows[row].itsFullResFlags;
    }
    std::vector<bool>::iterator getFullResFlags(const std::size_t row) {
      return itsRows[row].itsFullResFlags;
    }

    const std::vector<Row>& getRows() const {
      return itsRows;
    }

  private:
    /// Memory pools for the data in the rows.
    /// @{
    std::vector<std::complex<float>> itsData;
    std::vector<bool> itsFlags;
    std::vector<float> itsWeights;
    std::vector<bool> itsFullResFlags;
    /// @}

    /// The rows, which contain iterators to the memory pools above.
    std::vector<Row> itsRows;
  };

  }
}


#endif
