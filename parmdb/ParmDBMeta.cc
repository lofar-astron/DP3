//  ParmDBMeta.cc: one line description
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmDBMeta.h"

#include "../blob/BlobOStream.h"
#include "../blob/BlobIStream.h"

namespace dp3 {
namespace parmdb {

using dp3::blob::BlobIStream;

ParmDBMeta::ParmDBMeta() {}

ParmDBMeta::ParmDBMeta(const std::string& type, const std::string& tableName)
    : itsType(type), itsTableName(tableName) {}

void ParmDBMeta::setSQLMeta(const std::string& dbName,
                            const std::string& userName,
                            const std::string& dbPwd,
                            const std::string& hostName) {
  itsDBName = dbName;
  itsUserName = userName;
  itsDBPwd = dbPwd;
  itsHostName = hostName;
}

blob::BlobOStream& operator<<(blob::BlobOStream& bos, const ParmDBMeta& pdm) {
  bos << pdm.itsType << pdm.itsTableName << pdm.itsDBName << pdm.itsUserName
      << pdm.itsDBPwd << pdm.itsHostName;
  return bos;
}

blob::BlobIStream& operator>>(blob::BlobIStream& bis, ParmDBMeta& pdm) {
  bis >> pdm.itsType >> pdm.itsTableName >> pdm.itsDBName >> pdm.itsUserName >>
      pdm.itsDBPwd >> pdm.itsHostName;
  return bis;
}

}  // namespace parmdb
}  // namespace dp3
