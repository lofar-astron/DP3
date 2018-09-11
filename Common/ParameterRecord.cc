//# ParameterRecord.cc: A record of parameter values
//#
//# Copyright (C) 2012
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
//# $Id: ParameterRecord.cc 20264 2012-02-28 07:22:45Z diepen $

#include "ParameterRecord.h"

#include <cstdio>
#include <ostream>
#include <string>

namespace DP3 { 

  std::ostream& operator<< (std::ostream& os, const ParameterRecord& prec)
  {
    bool first = true;
    os << '{';
    for (ParameterRecord::const_iterator iter=prec.begin();
         iter!=prec.end(); ++iter) {
      if (first) {
        first = false;
      } else {
        os << ", ";
      }
      os << '\'' << iter->first << "': " << iter->second;
    }
    os << '}';
    return os;
  }

  bool ParameterRecord::getRecursive (const std::string& key,
                                      ParameterValue& value) const
  {
    const_iterator iter = find(key);
    if (iter != end()) {
      value = iter->second;
      return true;
    }
    // Try to find the key in possible ParmRecords.
    // Strip the last part from the key.
    std::string::size_type pos = key.rfind ('.');
    while (pos != std::string::npos) {
      const_iterator iter = find(key.substr(0, pos));
      if (iter != end()) {
        ParameterValue pv (iter->second);
        if (pv.isRecord()  &&
            pv.getRecord().getRecursive(key.substr(pos+1), value)) {
          return true;
        }
      }
      if (pos == 0) return false;
      pos = key.rfind ('.', pos-1);
    }
    return false;
  }

}
