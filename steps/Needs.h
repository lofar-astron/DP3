// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_NEEDS_H_
#define DP3_STEPS_NEEDS_H_

#include <bitset>

namespace dp3 {
namespace steps {

/**
 * Specifies which data types a Step uses.
 */
class Needs {
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
   * Creates a Needs object with no set needs.
   */
  constexpr Needs() : value_() {}

  /**
   * Creates a Needs object with a single need.
   * @param need The single need that must be set. Cannot be Single::kCount.
   */
  constexpr explicit Needs(Single need) : value_(1 << static_cast<int>(need)) {}

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
   * Adds the needs of another Needs object to the current object.
   * @param other Other needs value.
   * @return A reference to the current Needs object.
   */
  Needs& operator|=(const Needs& other) {
    value_ |= other.value_;
    return *this;
  }

  /**
   * Combines two Needs objects.
   * @return A Needs object with the needs of both arguments.
   */
  friend Needs operator|(const Needs& left, const Needs& right) {
    Needs needs(left);
    needs |= right;
    return needs;
  }

 private:
  // Using a bitset instead of individual booleans simplifies adding more
  // needs in the future.
  std::bitset<static_cast<int>(Single::kCount)> value_;
};

}  // namespace steps
}  // namespace dp3

#endif
