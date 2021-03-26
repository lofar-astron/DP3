// ParmDBLocker.h: Class to hold a read or write lock on ParmDBs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to hold a read or write lock on ParmDBs
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBLOCKER_H
#define LOFAR_PARMDB_PARMDBLOCKER_H

#include <memory>
#include <vector>

namespace dp3 {
namespace parmdb {

class ParmSet;
class ParmDB;

/// @ingroup ParmDB
/// @{

/// @brief Class to hold a read or write lock on ParmDBs

/// This class locks a single ParmDB or all ParmDBs used by a ParmSet.
/// Because the destructor does the unlocking, this class is very well
/// suited for automatically managing the locks. Even in case of an
/// exception, the locks are automatically released.
class ParmDBLocker {
 public:
  /// Define a shared pointer for this type.
  typedef std::shared_ptr<ParmDBLocker> ShPtr;

  /// Create a read or write lock on all ParmDBs in the ParmSet.
  explicit ParmDBLocker(const ParmSet& parmSet, bool write = false);

  /// Create a lock on a specific ParmDB.
  explicit ParmDBLocker(ParmDB& parmdb, bool write = false);

  /// The destructor unlocks the ParmDBs locked by the constructor.
  ~ParmDBLocker();

 private:
  /// Cannot copy.
  ///@{
  ParmDBLocker(const ParmDBLocker&);
  ParmDBLocker& operator=(const ParmDBLocker&);
  ///@}

  /// The locked DBs.
  std::vector<ParmDB*> itsDBs;
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
