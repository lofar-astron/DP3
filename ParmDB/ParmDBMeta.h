//#  ParmDBMeta.h: Meta information for the name and type of a ParmDB
//#
//#  Copyright (C) 2005
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
//#  $Id: ParmDBMeta.h 21605 2012-07-16 14:40:36Z diepen $

// @file
// @brief Meta information for the name and type of a ParmDB
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PARMDBMETA_H
#define LOFAR_PARMDB_PARMDBMETA_H

#include <string>

namespace DP3 {
//# Forward Declarations.
class BlobOStream;
class BlobIStream;

namespace BBS {

  // @ingroup ParmDB
  // @{

  // @brief Meta information for the name and type of a ParmDB
  class ParmDBMeta
  {
  public:
    ParmDBMeta();

    // Construct from a given type and file/table name.
    // The type can be empty, casa or blob.
    // If empty, the code doing the open will detect the exact type. At the
    // moment that can only be used for SourceDB.
    ParmDBMeta (const std::string& type, const std::string& tableName);

    void setSQLMeta (const std::string& dbName, const std::string& userName,
        const std::string& dbPwd, const std::string& hostName);

    const std::string& getType() const
      { return itsType; }

    const std::string& getTableName() const
      { return itsTableName; }

    const std::string& getDBName() const
      { return itsDBName; }
      
    const std::string& getUserName() const
      { return itsUserName; }
      
    const std::string& getDBPwd() const
      { return itsDBPwd; }
    
    const std::string& getHostName() const
      { return itsHostName; }
    
    // Write the object into a blob.
    friend BlobOStream& operator<< (BlobOStream&, const ParmDBMeta&);

    // Read the object from a blob.
    friend BlobIStream& operator>> (BlobIStream&, ParmDBMeta&);

  private:
    //# Datamembers
    std::string itsType;
    std::string itsTableName;
    // these options are used for sql databases
    std::string itsDBName;
    std::string itsUserName;
    std::string itsDBPwd;
    std::string itsHostName;
  };

  // @}

} // namespace BBS
} // namespace LOFAR

#endif
