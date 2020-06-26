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
          std::complex<float>* data,
          std::vector<bool>::iterator flags,
          float* weights,
          std::vector<bool>::iterator fullResFlags);

      double itsTime; ///< Start time for the measurements in MJD seconds.
      double itsExposure; ///< Exposure duration for the measurements in seconds.
      rownr_t itsRowNr;
      std::size_t itsNChannels;
      std::size_t itsNCorrelations;
      std::complex<float>* itsData; ///< Visibility measurements.
      std::vector<bool>::iterator itsFlags;
      float* itsWeights;
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
                const std::complex<double> *const data = nullptr,
                const bool *const flags = nullptr,
                const float *const weights = nullptr,
                const bool *const fullResFlags = nullptr,
                const double *const UVW = nullptr);

    const std::vector<std::complex<float>>& getData() const {
      return itsData;
    }
    std::vector<std::complex<float>>& getData() {
      return itsData;
    }

    const std::vector<Row>& getRows() const {
      return itsRows;
    }

    std::complex<float>* getData(std::size_t row) {
      return itsRows[row].itsData;
    }

  private:
    /// Memory pools for the data in the rows.
    /// @{
    std::vector<std::complex<float>> itsData;
    std::vector<bool> itsFlags;
    std::vector<float> itsWeights;
    std::vector<bool> itsFullResFlags;
    /// @}

    /// The rows, which contain pointers to the memory pools above.
    std::vector<Row> itsRows;
  };

  }
}


#endif
