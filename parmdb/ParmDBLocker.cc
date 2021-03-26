// ParmDBLocker.cc: Class to hold a read or write lock on ParmDBs
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmDBLocker.h"
#include "ParmDB.h"
#include "ParmSet.h"

namespace dp3 {
namespace parmdb {

ParmDBLocker::ParmDBLocker(const ParmSet& parmSet, bool lockForWrite)
    : itsDBs(parmSet.getDBs()) {
  for (unsigned int i = 0; i < itsDBs.size(); ++i) {
    itsDBs[i]->lock(lockForWrite);
  }
}

ParmDBLocker::ParmDBLocker(ParmDB& parmdb, bool lockForWrite)
    : itsDBs(1, &parmdb) {
  parmdb.lock(lockForWrite);
}

ParmDBLocker::~ParmDBLocker() {
  for (unsigned int i = 0; i < itsDBs.size(); ++i) {
    itsDBs[i]->unlock();
  }
}

}  // namespace parmdb
}  // namespace dp3
