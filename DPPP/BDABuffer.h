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

#include "../Common/EnumBitset.h"
#include "../Common/Types.h"
#include "../Common/UVector.h"

#include <complex>
#include <vector>

namespace DP3 {
  namespace DPPP {

    class BDABuffer
    {
    public:
      enum class Field { kData, kFlags, kWeights, kFullResFlags, kBitCount };
      using Fields = EnumBitset<Field>;

      struct Row {
        Row(double time,
            double interval,
            rownr_t row_nr,
            std::size_t baseline_nr,
            std::size_t n_channels,
            std::size_t n_correlations,
            std::complex<float>* data,
            bool* flags,
            float* weights,
            bool* fullResFlags,
            const double* uvw);

        std::size_t GetDataSize() const {
          return n_channels_ * n_correlations_;
        }

        const double time_; ///< Start time for the measurements in MJD seconds.
        const double interval_; ///< Duration time for the measurements in seconds.
        const rownr_t row_nr_;
        const std::size_t baseline_nr_;
        const std::size_t n_channels_;
        const std::size_t n_correlations_;
        std::complex<float>* const data_;
        bool* const flags_;
        float* const weights_;
        bool* const full_res_flags_;
        double uvw_[3];
      };

    public:
      /**
       * Create a new BDABuffer.
       * @param pool_size Size of the memory pool for this buffer.
       *                  (number of items)
       * @param fields The fields that should be enabled in this buffer.
       */
      explicit BDABuffer(std::size_t pool_size, Fields fields = Fields().Set());

      /**
       * Copy constructor.
       * This constructor sets the memory pool size of the new buffer to the
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

      /**
       * Determine the number of stored elements in all rows.
       * @return The total number of elements in this buffer.
       */
      std::size_t GetNumberOfElements() const;

      const std::complex<float>* GetData() const {
        return data_.empty() ? nullptr : data_.data();
      }
      std::complex<float>* GetData() {
        return data_.empty() ? nullptr : data_.data();
      }

      const bool* GetFlags() const {
        return flags_.empty() ? nullptr : flags_.data();
      }
      bool* GetFlags() {
        return flags_.empty() ? nullptr : flags_.data();
      }

      const float* GetWeights() const {
        return weights_.empty() ? nullptr : weights_.data();
      }
      float* GetWeights() {
        return weights_.empty() ? nullptr : weights_.data();
      }

      const bool* GetFullResFlags() const {
        return full_res_flags_.empty() ? nullptr : full_res_flags_.data();
      }
      bool* GetFullResFlags() {
        return full_res_flags_.empty() ? nullptr : full_res_flags_.data();
      }

      const std::complex<float>* GetData(std::size_t row)  const {
        return rows_[row].data_;
      }
      std::complex<float>* GetData(std::size_t row) {
        return rows_[row].data_;
      }

      const bool* GetFlags(std::size_t row) const {
        return rows_[row].flags_;
      }
      bool* GetFlags(std::size_t row) {
        return rows_[row].flags_;
      }

      const float* GetWeights(std::size_t row) const {
        return rows_[row].weights_;
      }
      float* GetWeights(std::size_t row) {
        return rows_[row].weights_;
      }

      const bool* GetFullResFlags(std::size_t row) const {
        return rows_[row].full_res_flags_;
      }
      bool* GetFullResFlags(std::size_t row) {
        return rows_[row].full_res_flags_;
      }

      const std::vector<Row>& GetRows() const {
        return rows_;
      }

    private:
      static constexpr double kTimeEpsilon = 1.0e-8; // For comparing measurement timestamps.

    public:
      static constexpr bool TimeIsLess(double x, double y) {
        return x < (y - kTimeEpsilon);
      }

      static constexpr bool TimeIsGreaterEqual(double x, double y) {
        return x > (y - kTimeEpsilon);
      }

      static constexpr bool TimeIsEqual(double x, double y) {
        return abs(x - y) < kTimeEpsilon;
      }

    private:
      /// Memory pools for the data in the rows. Since std::vector<bool>
      /// does not support pointers to its elements, use ao::uvector instead.
      /// @{
      ao::uvector<std::complex<float>> data_;
      ao::uvector<bool> flags_;
      ao::uvector<float> weights_;
      ao::uvector<bool> full_res_flags_;
      /// @}

      /// The rows, which contain iterators to the memory pools above.
      std::vector<Row> rows_;

      std::size_t original_capacity_; ///< Original capacity (number of items)
      std::size_t remaining_capacity_; ///< Remaining capacity (number of items)
    };

  }
}


#endif
