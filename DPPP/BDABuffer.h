// BDABuffer.h: Buffer holding BDA data
// Copyright (C) 2010
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

#include <complex>
#include <vector>

namespace DP3 {
  namespace DPPP {

  class BDABuffer
  {
  public:
    using Line = struct {
      double itsTime; ///< Start time for the measurements
      double itsExposure; ///< Exposure duration for the measurements
      std::complex<double>* itsMeasurements;
      bool* itsFlags;
      float* itsWeights;
      bool* itsFullResFlags;
      std::size_t itsNChannels;
      std::size_t itsNCorrelations;
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
    bool addLine(double time,
                 double exposure,
                 std::size_t nChannels,
                 std::size_t nCorrelations,
                 const std::complex<double>* measurements);

    const std::vector<std::complex<double>>& getMeasurements() const {
      return itsMeasurements;
    }
    std::vector<std::complex<double>>& getMeasurements() {
      return itsMeasurements;
    }

    const Line& getLine(std::size_t line) const {
      return itsLines[line];
    }

    std::complex<double>* getMeasurements( std::size_t line ) {
      return itsLines[line].itsMeasurements;
    }

  private:
    std::vector<std::complex<double>> itsMeasurements;
    std::vector<Line> itsLines;
  };

  }
}


#endif
