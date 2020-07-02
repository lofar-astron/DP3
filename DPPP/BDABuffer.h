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
            double interval,
            rownr_t row_nr,
            std::size_t baseline_nr,
            std::size_t n_channels,
            std::size_t n_correlations,
            std::vector<std::complex<float>>::iterator data,
            std::vector<bool>::iterator flags,
            std::vector<float>::iterator weights,
            std::vector<bool>::iterator fullResFlags,
            const double* uvw);

        const double time_; ///< Start time for the measurements in MJD seconds.
        const double interval_; ///< Duration time for the measurements in seconds.
        const rownr_t row_nr_;
        const std::size_t baseline_nr_;
        const std::size_t n_channels_;
        const std::size_t n_correlations_;
        const std::vector<std::complex<float>>::iterator data_;
        const std::vector<bool>::iterator flags_;
        const std::vector<float>::iterator weights_;
        const std::vector<bool>::iterator full_res_flags_;
        double uvw_[3];
      };
      
    public:
      /**
       * Create a new BDABuffer.
       * @param poolSize Size of the memory pool for this buffer
       *                 (number of complex values)
       */
      explicit BDABuffer(std::size_t pool_size);

      /**
       * Copy constructor.
       * This constructor sets the memory pool size to the
       * actual memory usage of the other buffer.
       * Adding new rows to the new buffer is not possible.
       * @param other An existing BDABuffer.
       */
      explicit BDABuffer(const BDABuffer& other);

      /**
       * Add a measurement line to the buffer.
       * @return True if the line is added.
       *         False if the buffer is full.
       */
      bool AddRow(double time,
                  double interval,
                  rownr_t row_nr,
                  std::size_t baseline_nr,
                  std::size_t n_channels,
                  std::size_t n_correlations,
                  const std::complex<float>* data = nullptr,
                  const bool* flags = nullptr,
                  const float* weights = nullptr,
                  const bool* full_res_flags = nullptr,
                  const double* uvw = nullptr);

      /**
       * Clears all data in the buffer.
       *
       * The memory pool capacity of the buffer remains unchanged, which allows
       * reusing the buffer.
       */
      void Clear();

      const std::vector<std::complex<float>>& GetData() const {
        return data_;
      }
      std::vector<std::complex<float>>& GetData() {
        return data_;
      }

      const std::vector<bool>& GetFlags() const {
        return flags_;
      }
      std::vector<bool>& GetFlags() {
        return flags_;
      }

      const std::vector<float>& GetWeights() const {
        return weights_;
      }
      std::vector<float>& GetWeights() {
        return weights_;
      }

      const std::vector<bool>& GetFullResFlags() const {
        return flags_;
      }
      std::vector<bool>& GetFullResFlags() {
        return flags_;
      }

      std::vector<std::complex<float>>::const_iterator GetData(std::size_t row)  const {
        return rows_[row].data_;
      }const
      std::vector<std::complex<float>>::iterator GetData(std::size_t row) {
        return rows_[row].data_;
      }

      std::vector<bool>::const_iterator GetFlags(std::size_t row) const {
        return rows_[row].flags_;
      }
      std::vector<bool>::iterator GetFlags(std::size_t row) {
        return rows_[row].flags_;
      }

      std::vector<float>::const_iterator GetWeights(std::size_t row) const {
        return rows_[row].weights_;
      }
      std::vector<float>::iterator GetWeights(std::size_t row) {
        return rows_[row].weights_;
      }

      std::vector<bool>::const_iterator GetFullResFlags(std::size_t row) const {
        return rows_[row].full_res_flags_;
      }
      std::vector<bool>::iterator GetFullResFlags(std::size_t row) {
        return rows_[row].full_res_flags_;
      }

      const std::vector<Row>& GetRows() const {
        return rows_;
      }

    private:
      /// Memory pools for the data in the rows.
      /// @{
      std::vector<std::complex<float>> data_;
      std::vector<bool> flags_;
      std::vector<float> weights_;
      std::vector<bool> full_res_flags_;
      /// @}

      /// The rows, which contain iterators to the memory pools above.
      std::vector<Row> rows_;
    };

  }
}


#endif
