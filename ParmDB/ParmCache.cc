//# ParmCache.cc: A class to cache ParmDB entries for a given work domain
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
//# $Id: ParmCache.cc 14038 2009-09-17 13:59:12Z diepen $

#include "ParmCache.h"
#include "ParmValue.h"
#include "ParmDBLocker.h"
#include "ParmDB.h"

namespace DP3 {
namespace BBS {

  ParmCache::ParmCache (ParmSet& parmSet)
    : itsParmSet (&parmSet)
  {}

  ParmCache::ParmCache (ParmSet& parmSet, const Box& workDomain)
    : itsParmSet    (&parmSet),
      itsWorkDomain (workDomain)
  {
    cacheValues();
  }

  void ParmCache::clear()
  {
    itsValueSets.clear();
    itsAxisCache.clear();
  }

  void ParmCache::reset (const Box& workDomain)
  {
    clear();
    itsWorkDomain = workDomain;
    cacheValues();
  }

  void ParmCache::cacheValues()
  {
    if (itsParmSet->size() > itsValueSets.size()) {
      itsParmSet->getValues (itsValueSets, itsWorkDomain);
    }
  }

  void ParmCache::setSolveGrid (ParmId parmId, const Grid& solveGrid)
  {
    assert (parmId < itsValueSets.size());
    itsValueSets[parmId].setSolveGrid (solveGrid);
  }

  void ParmCache::flush()
  {
    ParmDBLocker (*itsParmSet, true);
    for (unsigned int i=0; i<itsValueSets.size(); ++i) {
      ParmValueSet& pvset = itsValueSets[i];
      if (pvset.isDirty()) {
        itsParmSet->write (i, pvset);
        pvset.setDirty (false);
      }
    }
  }

} //# end namespace BBS
} //# end namspace LOFAR
