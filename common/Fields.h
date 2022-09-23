// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_FIELDS_H_
#define DP3_COMMON_FIELDS_H_

#include <bitset>

namespace dp3 {
namespace common {

/**
 * Specifies which data types a Step uses.
 */
class Fields {
 public:
  /**
   * Values for specifying a single need.
   */
  enum class Single {
    kData,          ///< Is the visibility data needed?
    kFlags,         ///< Are the flags needed?
    kWeights,       ///< Are the weights needed?
    kFullResFlags,  ///< Are the full res flags needed?
    kUvw,           ///< Are the uvw needed?
    kCount          ///< Number of need values. Must be last.
  };

  /**
   * Creates a Fields object with no set fields.
   */
  constexpr Fields() : value_() {}

  /**
   * Creates a Fields object with a single field.
   * @param field The single field that must be set. Cannot be Single::kCount.
   */
  constexpr explicit Fields(Single field)
      : value_(1 << static_cast<int>(field)) {}

  /**
   * @return True if visibility data is needed, false if not.
   */
  constexpr bool Data() const {
    return value_[static_cast<int>(Single::kData)];
  }

  /**
   * @return True if flags are needed, false if not.
   */
  constexpr bool Flags() const {
    return value_[static_cast<int>(Single::kFlags)];
  }

  /**
   * @return True if weights are needed, false if not.
   */
  constexpr bool Weights() const {
    return value_[static_cast<int>(Single::kWeights)];
  }

  /**
   * @return True if full res flags are needed, false if not.
   */
  constexpr bool FullResFlags() const {
    return value_[static_cast<int>(Single::kFullResFlags)];
  }

  /**
   * @return True if uvw values are needed, false if not.
   */
  constexpr bool Uvw() const { return value_[static_cast<int>(Single::kUvw)]; }

  /**
   * Adds the fields of another Fields object to the current object.
   * @param other Other fields value.
   * @return A reference to the current Fields object.
   */
  Fields& operator|=(const Fields& other) {
    value_ |= other.value_;
    return *this;
  }

  /**
   * Combines two Fields objects.
   * @return A Fields object with the fields of both arguments.
   */
  friend Fields operator|(const Fields& left, const Fields& right) {
    Fields fields(left);
    fields |= right;
    return fields;
  }

 private:
  // Using a bitset instead of individual booleans simplifies adding more
  // fields in the future.
  std::bitset<static_cast<int>(Single::kCount)> value_;
};

}  // namespace common
}  // namespace dp3

#endif
