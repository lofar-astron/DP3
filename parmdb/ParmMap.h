// ParmMap.h: A map of parameter name to value set
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief A map of parameter name to value set.
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMMAP_H
#define LOFAR_PARMDB_PARMMAP_H

#include "ParmValue.h"

#include <map>

namespace dp3 {
namespace parmdb {

class ParmDB;

/// @ingroup ParmDB
/// @{

/// @brief A map of parameter name to value set.

/// ParmMap holds a map of name to ParmValueSet.
/// It is meant to hold the default values, but could be used for
/// other purposes as well.
class ParmMap {
 public:
  /// Set up a map for the given domain in the ParmDB.
  ParmMap() {}

  /// Add or replace a ParmValueSet.
  void define(const std::string& name, const ParmValueSet& pset) {
    itsValueSets[name] = pset;
  }

  /// Is the map empty?
  bool empty() const { return itsValueSets.empty(); }

  /// Return the size of the map.
  unsigned int size() const { return itsValueSets.size(); }

  /// Get the value belonging to the name.
  /// An exception is thrown if the value does not exist.
  const ParmValueSet& operator[](const std::string& name) const;

  /// Iterator functionality.
  ///@{
  typedef std::map<std::string, ParmValueSet>::iterator iterator;
  typedef std::map<std::string, ParmValueSet>::const_iterator const_iterator;
  typedef std::map<std::string, ParmValueSet>::value_type valueType;
  iterator begin() { return itsValueSets.begin(); }
  const_iterator begin() const { return itsValueSets.begin(); }
  iterator end() { return itsValueSets.end(); }
  const_iterator end() const { return itsValueSets.end(); }
  iterator find(const std::string& name) { return itsValueSets.find(name); }
  const_iterator find(const std::string& name) const {
    return itsValueSets.find(name);
  }
  ///@}

  /// Clear the map.
  void clear() { itsValueSets.clear(); }

 private:
  /// Cannot copy.
  ///@{
  ParmMap(const ParmMap&);
  ParmMap& operator=(const ParmMap&);
  ///@}

  std::map<std::string, ParmValueSet> itsValueSets;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
