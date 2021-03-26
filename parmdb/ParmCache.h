// ParmCache.h: A class dealing with caching and handling ParmDB entries
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief A class dealing with caching and handling ParmDB entries
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMCACHE_H
#define LOFAR_PARMDB_PARMCACHE_H

#include "ParmSet.h"
#include "ParmValue.h"
#include "AxisMapping.h"

namespace dp3 {
namespace parmdb {

class ParmDB;

/// @ingroup ParmDB
/// @{

/// @brief A class dealing with caching and handling ParmDB entries

/// ParmCache caches the parm records for a given work domain.
class ParmCache {
 public:
  /// Set up a cache. Attach it to he given ParmSet.
  /// No work domain is set yet; this has to be done later using the
  /// reset function.
  ParmCache(ParmSet&);

  /// Set up a cache for the given work domain.
  /// Attach it to the given ParmSet.
  /// Note that nameIds in the ParmSet might get changed when new values are
  /// written.
  ParmCache(ParmSet&, const Box& workDomain);

  /// Get access to the underlying ParmSet.
  ///@{
  ParmSet& getParmSet() { return *itsParmSet; }
  const ParmSet& getParmSet() const { return *itsParmSet; }
  ///@}

  /// Clear the cache.
  void clear();

  /// Reset the work domain which will clear the cache.
  /// A new prefetch should be done the get the values for the new work domain.
  void reset(const Box& workDomain);

  /// Cache the values of the parameters in the attached ParmSet for the
  /// current work domain.
  /// It will only do it for the parameters not prefetched yet.
  /// It makes it possible to prefetch the values of parameters added to
  /// the ParmSet since the last cacheValues.
  void cacheValues();

  /// Get the value set for the given parm.
  ///@{
  ParmValueSet& getValueSet(ParmId parmid) { return itsValueSets[parmid]; }
  const ParmValueSet& getValueSet(ParmId parmid) const {
    return itsValueSets[parmid];
  }
  ///@}

  /// Get the AxisMappingCache object.
  AxisMappingCache& getAxisMappingCache() { return itsAxisCache; }

  /// Check for a solvable parm if the domains in the value set match the
  /// given solve domains. The solve domains can exceed the work domain.
  /// If they exceed, they are limited to the work domain.
  /// If the value set has only the default value, it is replaced by a
  /// value set with values for each solve domain.
  void setSolveGrid(ParmId parmid, const Grid& solveGrid);

  /// Flush the cache which means that all updates are written into the
  /// appropriate ParmDB tables.
  void flush();

 private:
  /// Cannot copy.
  ///@{
  ParmCache(const ParmCache&);
  ParmCache& operator=(const ParmCache&);
  ///@}

  ParmSet* itsParmSet;
  Box itsWorkDomain;
  std::vector<ParmValueSet> itsValueSets;
  AxisMappingCache itsAxisCache;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
