// BdaBuffer.h: Buffer holding BDA data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Buffer holding baseline-dependently averaged (BDA) data.
/// @author Maik Nijhuis and Lars Krombeen

#ifndef DP3_BASE_BDABUFFER_H_
#define DP3_BASE_BDABUFFER_H_

#include <complex>
#include <map>
#include <string>
#include <vector>

#include <aocommon/uvector.h>

#include <dp3/common/Fields.h>
#include <dp3/common/Types.h>

namespace dp3 {
namespace base {

class BdaBuffer {
 public:
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
   *               The 'data' field enables the main visibility buffer.
   *               The 'uvw' field is always enabled in a BdaBuffer.
   */
  explicit BdaBuffer(std::size_t pool_size, const common::Fields& fields);

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
   * If 'other' does not have the field, memory is allocated for the field and
   * the content is not initialized.
   */
  BdaBuffer(const BdaBuffer& other, const common::Fields& fields);

  /**
   * Adds a visibility buffer.
   * If the buffer already exists, nothing happens.
   * If the buffer is created, its visibility values are not initialized.
   * @param name Name for the new buffer.
   *        If empty, adds the main data buffer. The effect is then equal to
   *        enabling the 'data' field.
   */
  void AddData(const std::string& name = "") {
    aocommon::UVector<std::complex<float>>& new_data = data_[name];
    new_data.reserve(original_capacity_);
    new_data.resize(GetNumberOfElements());
  }

  /**
   * Moves a visibility buffer into the current BdaBuffer.
   * If the target visibility buffer already exists, it is overwritten.
   * @param source The source BdaBuffer. It may equal the current BdaBuffer.
   *               In that case, rename a visibility buffer.
   * @param source_name Name of a visibility buffer in 'source'.
   * @param target_name Name of the target visibility buffer.
   */
  void MoveData(BdaBuffer& source, const std::string& source_name = "",
                const std::string& target_name = "");

  /**
   * Removes a visibility buffer.
   * If the buffer does not exist, nothing happens.
   * @param name Name of the buffer to remove.
   *        If empty, removes the main data buffer. The effect is then equal to
   *        disabling the 'data' field.
   */
  void RemoveData(const std::string& name = "") { data_.erase(name); }

  /**
   * @return If the BdaBuffer has a visibility buffer for the given name.
   */
  bool HasData(const std::string& name = "") const {
    return data_.find(name) != data_.end();
  }

  /** @return The sorted names of all visibility buffers. */
  std::vector<std::string> GetDataNames() const {
    std::vector<std::string> names;
    names.reserve(data_.size());
    for (const auto& name_vector : data_) {
      names.push_back(name_vector.first);
    }
    return names;
  }

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
   * @param data Pointer to visibilities for the row.
   *             If null, initializes the visibilities to zero.
   * @param flags Pointer to flags for the row.
   *              If null, initializes the flags to false.
   * @param weights Pointer to weights for the row.
   *                If null, initializes the weights to zero.
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
   * Removes all rows from the buffer.
   *
   * The memory pool capacity of the buffer remains unchanged. All data buffers
   * also remain. Keeping all internal buffers allows reusing the BdaBuffer.
   */
  void Clear();

  /**
   * Determine the number of stored elements in all rows.
   * @return The total number of elements in this buffer.
   */
  std::size_t GetNumberOfElements() const {
    return original_capacity_ - remaining_capacity_;
  }

  /**
   * Determine the remaining capacity.
   * @return The remaining capacity (in number of elements) for this buffer.
   */
  std::size_t GetRemainingCapacity() const { return remaining_capacity_; }

  const std::complex<float>* GetData(const std::string& name = "") const;
  std::complex<float>* GetData(const std::string& name = "");

  const std::complex<float>* GetData(std::size_t row,
                                     const std::string& name = "") const;
  std::complex<float>* GetData(std::size_t row, const std::string& name = "");

  const bool* GetFlags() const {
    return flags_.empty() ? nullptr : flags_.data();
  }
  bool* GetFlags() { return flags_.empty() ? nullptr : flags_.data(); }

  const bool* GetFlags(std::size_t row) const {
    return flags_.empty() ? nullptr : flags_.data() + rows_[row].offset;
  }
  bool* GetFlags(std::size_t row) {
    return flags_.empty() ? nullptr : flags_.data() + rows_[row].offset;
  }

  const float* GetWeights() const {
    return weights_.empty() ? nullptr : weights_.data();
  }
  float* GetWeights() { return weights_.empty() ? nullptr : weights_.data(); }

  const float* GetWeights(std::size_t row) const {
    return weights_.empty() ? nullptr : weights_.data() + rows_[row].offset;
  }
  float* GetWeights(std::size_t row) {
    return weights_.empty() ? nullptr : weights_.data() + rows_[row].offset;
  }

  std::vector<Row>& GetRows() { return rows_; }
  const std::vector<Row>& GetRows() const { return rows_; }

  static constexpr bool TimeIsLess(double x, double y) {
    return x < (y - kTimeTolerance);
  }
  static constexpr bool TimeIsLessEqual(double x, double y) {
    return x < (y + kTimeTolerance);
  }
  static constexpr bool TimeIsGreaterEqual(double x, double y) {
    return x > (y - kTimeTolerance);
  }
  static constexpr bool TimeIsEqual(double x, double y) {
    // Don't use std::fabs, since it is not a constexpr.
    return ((x > y) ? (x - y) : (y - x)) < kTimeTolerance;
  }

  bool IsMetadataEqual(const BdaBuffer& other) const;

  // Tolerance value for comparing measurement timestamps. It should be larger
  // than the rounding errors that occur due to BDA averaging and expanding.
  static constexpr double kTimeTolerance = 2.0e-6;

 private:
  /// Memory pools for the data in the rows. Since std::vector<bool>
  /// does not support pointers to its elements, use aocommon::UVector instead.
  /// @{
  std::map<std::string, aocommon::UVector<std::complex<float>>> data_;
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
