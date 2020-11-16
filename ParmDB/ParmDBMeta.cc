//  ParmDBMeta.cc: one line description
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ParmDBMeta.h"

#include "../Blob/BlobOStream.h"
#include "../Blob/BlobIStream.h"

namespace DP3 {
namespace BBS {

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

BlobOStream& operator<<(BlobOStream& bos, const ParmDBMeta& pdm) {
  bos << pdm.itsType << pdm.itsTableName << pdm.itsDBName << pdm.itsUserName
      << pdm.itsDBPwd << pdm.itsHostName;
  return bos;
}

BlobIStream& operator>>(BlobIStream& bis, ParmDBMeta& pdm) {
  bis >> pdm.itsType >> pdm.itsTableName >> pdm.itsDBName >> pdm.itsUserName >>
      pdm.itsDBPwd >> pdm.itsHostName;
  return bis;
}

}  // namespace BBS
}  // namespace DP3
