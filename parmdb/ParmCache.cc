// ParmCache.cc: A class to cache ParmDB entries for a given work domain
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmCache.h"
#include "ParmValue.h"
#include "ParmDBLocker.h"
#include "ParmDB.h"

namespace dp3 {
namespace parmdb {

ParmCache::ParmCache(ParmSet& parmSet) : itsParmSet(&parmSet) {}

ParmCache::ParmCache(ParmSet& parmSet, const Box& workDomain)
    : itsParmSet(&parmSet), itsWorkDomain(workDomain) {
  cacheValues();
}

void ParmCache::clear() {
  itsValueSets.clear();
  itsAxisCache.clear();
}

void ParmCache::reset(const Box& workDomain) {
  clear();
  itsWorkDomain = workDomain;
  cacheValues();
}

void ParmCache::cacheValues() {
  if (itsParmSet->size() > itsValueSets.size()) {
    itsParmSet->getValues(itsValueSets, itsWorkDomain);
  }
}

void ParmCache::setSolveGrid(ParmId parmId, const Grid& solveGrid) {
  assert(parmId < itsValueSets.size());
  itsValueSets[parmId].setSolveGrid(solveGrid);
}

void ParmCache::flush() {
  ParmDBLocker(*itsParmSet, true);
  for (unsigned int i = 0; i < itsValueSets.size(); ++i) {
    ParmValueSet& pvset = itsValueSets[i];
    if (pvset.isDirty()) {
      itsParmSet->write(i, pvset);
      pvset.setDirty(false);
    }
  }
}

}  // namespace parmdb
}  // namespace dp3
