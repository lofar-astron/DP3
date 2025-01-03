// BdaBuffer.h: Buffer holding BDA data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Buffer holding baseline-dependently averaged (BDA) data.
/// @author Maik Nijhuis and Lars Krombeen

#ifndef DP3_BASE_BDABUFFER_H_
#define DP3_BASE_BDABUFFER_H_

#include <dp3/common/Types.h>

#include <aocommon/uvector.h>
#include <complex>
#include <vector>

namespace dp3 {
namespace base {

class BdaBuffer {
 public:
  /**
   * Parameter structure for indicating which buffer elements are enabled.
   */
  struct Fields {
    /**
     * This constructor is necessary because of bugs in gcc and clang.
     * See
     * https://stackoverflow.com/questions/53408962/try-to-understand-compiler-error-message-default-member-initializer-required-be
     */
    explicit Fields(bool default_value = true)
        : data(default_value), flags(default_value), weights(default_value) {}
    bool data;     ///< Enable/Disable visibilities.
    bool flags;    ///< Enable/Disable flags.
    bool weights;  ///< Enable/Disable weights.
  };

  struct Row {
    Row(double time, double interval, double exposure, common::rownr_t row_nr,
        std::size_t baseline_nr, std::size_t n_channels,
        std::size_t n_correlations, std::size_t offset, const double* uvw);
    std::size_t GetDataSize() const { return n_channels * n_correlations; }
    bool IsMetadataEqual(const BdaBuffer::Row& other) const;
    const double time;  ///< Centroid time for the measurements in MJD seconds.
    const double interval;  ///< Duration time for the measurements in seconds.
    const double exposure;  ///< Exposure time for the measurements in seconds.
    common::rownr_t row_nr;
    const std::size_t baseline_nr;
    const std::size_t n_channels;
    const std::size_t n_correlations;
    const std::size_t offset;  ///< Relative position in BdaBuffer vectors.
    double uvw[3];
  };

  /**
   * Create a new BdaBuffer.
   * @param pool_size Size of the memory pool for this buffer.
   *                  (number of items)
   * @param fields The fields that should be enabled in this buffer.
   */
  explicit BdaBuffer(std::size_t pool_size, const Fields& fields = Fields());

  /**
   * Disabling the default copy constructor and copy assignment operator avoids
   * expensive implicit copies.
   * @{
   */
  BdaBuffer(const BdaBuffer& other) = delete;
  BdaBuffer& operator=(const BdaBuffer& other) = delete;
  /** @} */

  /**
   * Custom copy constructor.
   * This constructor sets the memory pool size of the new buffer to the
   * actual memory usage of the other buffer.
   * Adding new rows to the new buffer is therefore not possible.
   * @param other An existing BdaBuffer.
   * @param fields The fields that are enabled in the new buffer.
   * If 'other' does not have the field, memory is allocated for the field but
   * the content are not initialized.
   * @param copy_fields The fields that are copied to the new buffer.
   * Fields that are not enabled or are never copied.
   * For enabled, but not copied fields, the memory is allocated but not
   * initialized.
   */
  BdaBuffer(const BdaBuffer& other, const Fields& fields,
            const Fields& copy_fields = Fields());

  /**
   * Add a measurement line to the buffer.
   *
   * Measurement lines have to obey the following ordering constraint:
   * If a row starts at time T, all rows that end before or at T must be
   * added before this row. A new row thus may not have an end time before
   * or equal to the start time of the last row.
   *
   * Use GetRemainingCapacity() for checking if the buffer has enough space.
   *
   * @return True if the line is added.
   *         False if the buffer is full.
   * @throw std::invalid_argument If the row ordering is incorrect.
   */
  bool AddRow(double time, double interval, double exposure,
              std::size_t baseline_nr, std::size_t n_channels,
              std::size_t n_correlations,
              const std::complex<float>* data = nullptr,
              const bool* flags = nullptr, const float* weights = nullptr,
              const double* uvw = nullptr);

  /**
   * Update the row numbers of the rows in this buffer.
   * Does nothing if the buffer is empty.
   * @param base_rownr The row number for the first row in this buffer.
   *        The following rows get base_rownr + 1, base_rownr + 2, etc. as
   *        their row number.
   */
  void SetBaseRowNr(common::rownr_t base_rownr);

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

  /**
   * Determine the remaining capacity.
   * @return The remaining capacity (in number of elements) for this buffer.
   */
  std::size_t GetRemainingCapacity() const { return remaining_capacity_; }

  const std::complex<float>* GetData() const {
    return data_.empty() ? nullptr : data_.data();
  }
  std::complex<float>* GetData() {
    return data_.empty() ? nullptr : data_.data();
  }
  const bool* GetFlags() const {
    return flags_.empty() ? nullptr : flags_.data();
  }
  bool* GetFlags() { return flags_.empty() ? nullptr : flags_.data(); }
  const float* GetWeights() const {
    return weights_.empty() ? nullptr : weights_.data();
  }
  float* GetWeights() { return weights_.empty() ? nullptr : weights_.data(); }

  const std::complex<float>* GetData(std::size_t row) const {
    return data_.empty() ? nullptr : data_.data() + rows_[row].offset;
  }
  std::complex<float>* GetData(std::size_t row) {
    return data_.empty() ? nullptr : data_.data() + rows_[row].offset;
  }
  const bool* GetFlags(std::size_t row) const {
    return flags_.empty() ? nullptr : flags_.data() + rows_[row].offset;
  }
  bool* GetFlags(std::size_t row) {
    return flags_.empty() ? nullptr : flags_.data() + rows_[row].offset;
  }
  const float* GetWeights(std::size_t row) const {
    return weights_.empty() ? nullptr : weights_.data() + rows_[row].offset;
  }
  float* GetWeights(std::size_t row) {
    return weights_.empty() ? nullptr : weights_.data() + rows_[row].offset;
  }
  std::vector<Row>& GetRows() { return rows_; }
  const std::vector<Row>& GetRows() const { return rows_; }

  static constexpr bool TimeIsLess(double x, double y) {
    return x < (y - kTimeEpsilon);
  }
  static constexpr bool TimeIsLessEqual(double x, double y) {
    return x < (y + kTimeEpsilon);
  }
  static constexpr bool TimeIsGreaterEqual(double x, double y) {
    return x > (y - kTimeEpsilon);
  }
  static constexpr bool TimeIsEqual(double x, double y) {
    // Don't use std::fabs, since it is not a constexpr.
    return ((x > y) ? (x - y) : (y - x)) < kTimeEpsilon;
  }

  bool IsMetadataEqual(const BdaBuffer& other) const;

 private:
  static constexpr double kTimeEpsilon =
      1.0e-8;  // For comparing measurement timestamps.

  /// Memory pools for the data in the rows. Since std::vector<bool>
  /// does not support pointers to its elements, use aocommon::UVector instead.
  /// @{
  aocommon::UVector<std::complex<float>> data_;
  aocommon::UVector<bool> flags_;
  aocommon::UVector<float> weights_;
  /// @}
  /// The rows, which contain pointers to the memory pools above.
  std::vector<Row> rows_;
  std::size_t original_capacity_;   ///< Original capacity (number of items)
  std::size_t remaining_capacity_;  ///< Remaining capacity (number of items)
};

}  // namespace base
}  // namespace dp3

#endif  // DP3_BASE_BDABUFFER_H_
