//  ParmDBMeta.h: Meta information for the name and type of a ParmDB
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Meta information for the name and type of a ParmDB
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBMETA_H
#define LOFAR_PARMDB_PARMDBMETA_H

#include <string>

namespace dp3 {
namespace blob {
class BlobOStream;
class BlobIStream;
}  // namespace blob

namespace parmdb {

/// @ingroup ParmDB
/// @{

/// @brief Meta information for the name and type of a ParmDB
class ParmDBMeta {
 public:
  ParmDBMeta();

  /// Construct from a given type and file/table name.
  /// The type can be empty, casa or blob.
  /// If empty, the code doing the open will detect the exact type. At the
  /// moment that can only be used for SourceDB.
  ParmDBMeta(const std::string& type, const std::string& tableName);

  void setSQLMeta(const std::string& dbName, const std::string& userName,
                  const std::string& dbPwd, const std::string& hostName);

  const std::string& getType() const { return itsType; }

  const std::string& getTableName() const { return itsTableName; }

  const std::string& getDBName() const { return itsDBName; }

  const std::string& getUserName() const { return itsUserName; }

  const std::string& getDBPwd() const { return itsDBPwd; }

  const std::string& getHostName() const { return itsHostName; }

  /// Write the object into a blob.
  friend blob::BlobOStream& operator<<(blob::BlobOStream&, const ParmDBMeta&);

  /// Read the object from a blob.
  friend blob::BlobIStream& operator>>(blob::BlobIStream&, ParmDBMeta&);

 private:
  std::string itsType;
  std::string itsTableName;
  /// these options are used for sql databases
  ///@{
  std::string itsDBName;
  std::string itsUserName;
  std::string itsDBPwd;
  std::string itsHostName;
  ///@}
};

/// @}

}  // namespace parmdb
}  // namespace dp3

#endif
