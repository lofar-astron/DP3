// SourceDB.cc: Object to hold parameters in a table.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceDB.h"
#include "SourceDBCasa.h"
#include "SourceDBBlob.h"
#include "ParmDB.h"

#include <casacore/casa/OS/File.h>

namespace dp3 {
namespace parmdb {

SourceDBRep::SourceDBRep(const ParmDBMeta& ptm, bool forceNew)
    : itsCount(0), itsParmDB(ptm, forceNew) {}

SourceDBRep::~SourceDBRep() {}

void SourceDBRep::lock(bool) {}

void SourceDBRep::unlock() {}

SourceDB::SourceDB(const ParmDBMeta& ptm, bool mustExist, bool forceNew)
    : itsSourceDBPath(boost::filesystem::path(ptm.getTableName())) {
  if (mustExist && !casacore::File(ptm.getTableName()).exists())
    throw std::runtime_error("The sourcedb '" + ptm.getTableName() +
                             "' does not exist");
  ParmDBMeta pm(ptm);
  // Determine type if not given.
  // Default is casa, but an existing regular file is blob.
  if (pm.getType().empty()) {
    pm = ParmDBMeta("casa", pm.getTableName());
    if (!forceNew) {
      // Check if an existing DB is stored as a file (thus as SourceDBBlob).
      // This is for compatibility reasons.
      casacore::File file(ptm.getTableName());
      if (file.exists() && file.isRegular()) {
        pm = ParmDBMeta("blob", pm.getTableName());
      }
    }
  }
  if (pm.getType() == "casa") {
    itsRep = new SourceDBCasa(pm, forceNew);
  } else if (pm.getType() == "blob") {
    itsRep = new SourceDBBlob(pm, forceNew);
  } else {
    throw std::runtime_error("unknown sourceTableType: " + pm.getType());
  }
  itsRep->link();
}

SourceDB::~SourceDB() { decrCount(); }

SourceDB::SourceDB(SourceDBRep* rep) : itsRep(rep) { itsRep->link(); }

SourceDB::SourceDB(const SourceDB& that) : itsRep(that.itsRep) {
  itsRep->link();
}

SourceDB& SourceDB::operator=(const SourceDB& that) {
  if (this != &that) {
    decrCount();
    itsRep = that.itsRep;
    itsRep->link();
  }
  return *this;
}

void SourceDB::decrCount() {
  if (itsRep->unlink() == 0) {
    delete itsRep;
    itsRep = 0;
  }
}

}  // namespace parmdb
}  // namespace dp3
