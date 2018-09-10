//# ParmDBLocker.h: Class to hold a read or write lock on ParmDBs
//#
//# Copyright (C) 2008
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: ParmDBLocker.h 14038 2009-09-17 13:59:12Z diepen $

// @file
// @brief Class to hold a read or write lock on ParmDBs
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBLOCKER_H
#define LOFAR_PARMDB_PARMDBLOCKER_H

#include <Common/lofar_smartptr.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
namespace BBS {

  //# Forward Declarations
  class ParmSet;
  class ParmDB;

  // @ingroup ParmDB
  // @{

  // @brief Class to hold a read or write lock on ParmDBs
  // This class locks a single ParmDB or all ParmDBs used by a ParmSet.
  // Because the destructor does the unlocking, this class is very well
  // suited for automatically managing the locks. Even in case of an
  // exception, the locks are automatically released.
  class ParmDBLocker
  {
  public:
    // Define a shared pointer for this type.
    typedef shared_ptr<ParmDBLocker> ShPtr;

    // Create a read or write lock on all ParmDBs in the ParmSet.
    explicit ParmDBLocker (const ParmSet& parmSet, bool write=false);

    // Create a lock on a specific ParmDB.
    explicit ParmDBLocker (ParmDB& parmdb, bool write=false);

    // The destructor unlocks the ParmDBs locked by the constructor.
    ~ParmDBLocker();

  private:
    // Cannot copy.
    // <group>
    ParmDBLocker (const ParmDBLocker&);
    ParmDBLocker& operator= (const ParmDBLocker&);
    // </group>

    //# The locked DBs.
    vector<ParmDB*> itsDBs;
  };

  // @}

} //# end namespace BBS
} //# end namspace LOFAR

#endif
