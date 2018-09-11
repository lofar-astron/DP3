//#  ParmDBMeta.cc: one line description
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
//#  $Id: ParmDBMeta.cc 14038 2009-09-17 13:59:12Z diepen $

#include "ParmDBMeta.h"

#include "../Blob/BlobOStream.h"
#include "../Blob/BlobIStream.h"

namespace DP3 {
namespace BBS {

  ParmDBMeta::ParmDBMeta()
  {}

  ParmDBMeta::ParmDBMeta (const std::string& type,
                          const std::string& tableName)
    : itsType      (type),
      itsTableName (tableName)
  {}

  void ParmDBMeta::setSQLMeta (const std::string& dbName,
                               const std::string& userName,
                               const std::string& dbPwd,
                               const std::string& hostName)
  {
    itsDBName   = dbName;
    itsUserName = userName;
    itsDBPwd    = dbPwd;
    itsHostName = hostName;
  }

  BlobOStream& operator<< (BlobOStream& bos, const ParmDBMeta& pdm)
  {
    bos << pdm.itsType   << pdm.itsTableName
        << pdm.itsDBName << pdm.itsUserName
        << pdm.itsDBPwd  << pdm.itsHostName;
    return bos;
  }
    
  BlobIStream& operator>> (BlobIStream& bis, ParmDBMeta& pdm)
  {
    bis >> pdm.itsType >> pdm.itsTableName
        >> pdm.itsDBName >> pdm.itsUserName
        >> pdm.itsDBPwd >> pdm.itsHostName;
    return bis;
  }   

} // namespace BBS
} // namespace LOFAR
