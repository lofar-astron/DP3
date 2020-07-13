// EnumBitset.h: Bitset which uses an enumeration type as indices.
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
/// @brief Bitset which uses an enumeration type as indices.
/// @author Maik Nijhuis

#include <bitset>
#include <initializer_list>

namespace DP3
{
  /**
   * Helper class which simplifies creating a bitset from an enumeration type,
   * by using enumeration type values as bit positions.
   * Methods that do not use bit positions have the same interface as std::bitset.
   * @tparam E An enumeration type. It should have a kBitCount member, which
   *         holds the number of representable values in the bitset.
   */
  template<class E>
  class EnumBitset
  {
  private:
    using BitSetType = std::bitset<static_cast<std::size_t>(E::kBitCount)>;

  public:
    /**
     * Constructor. Initializes all bits to false.
     */
    constexpr EnumBitset()
    : bitset_()
    {}

    /**
     * Constructor using an initializer list.
     * @param list A list with values for the bitset. If the list size is
     *        smaller than the bitset size, remaining bits become false.
     */
    explicit EnumBitset(std::initializer_list<bool> list)
    : bitset_()
    {
      std::size_t i = 0;
      auto list_it = list.begin();
      while(i < static_cast<std::size_t>(E::kBitCount) &&
            list_it != list.end()) {
        bitset_.set(i, *list_it);;
        ++i;
        ++list_it;
      }
    }

    constexpr bool operator[](E bit) const
    {
      return bitset_[static_cast<std::size_t>(bit)];
    }

    typename BitSetType::reference operator[](E bit)
    {
      return bitset_[static_cast<std::size_t>(bit)];
    }

    /**
     * Set all bits to true.
     * @return A reference to the bitset.
     */
    EnumBitset& Set() { bitset_.set(); return *this; }

  private:
    BitSetType bitset_;
  };

} // end namespace DP3
