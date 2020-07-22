// DPLogger.h: Log on cout/cerr or through the logging system
// Copyright (C) 2010
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Log on cout/cerr or through the logging system
/// @author Ger van Diepen

#ifndef DPPP_DPLOGGER_H
#define DPPP_DPLOGGER_H

#include <iostream>

namespace DP3 {
namespace DPPP {

/// @brief This class contains the flag to choose between cout/cerr and logging
/// system.
class DPLogger {
 public:
  static bool useLogger;
};
}  // namespace DPPP
}  // namespace DP3

/// Log an informational message.
#define DPLOG_INFO_STR(stream) std::cout << stream << std::endl;

/// Log a fatal message.
#define DPLOG_WARN_STR(stream) std::cerr << stream << std::endl;

/// Log an informational message.
#define DPLOG_INFO(msg, removeEndl)                                  \
  std::string str(msg);                                              \
  if (removeEndl && str.size() > 0 && str[str.size() - 1] == '\n') { \
    str = str.substr(0, str.size() - 1);                             \
  }                                                                  \
  std::cout << str << std::endl;

#define LOGCOUT(msg)                    \
  {                                     \
    std::ostringstream ostr;            \
    ostr << msg;                        \
    printf("%s\n", ostr.str().c_str()); \
  }

#endif
