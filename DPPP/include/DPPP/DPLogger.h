//# DPLogger.h: Log on cout/cerr or through the logging system
//# Copyright (C) 2010
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
//# $Id$
//#
//# @author Ger van Diepen

#ifndef DPPP_DPLOGGER_H
#define DPPP_DPLOGGER_H

// @file
// @brief Log on cout/cerr or through the logging system

#include <Common/LofarLogger.h>
#include <iostream>

namespace LOFAR {
  namespace DPPP {

    // This class contains the flag to choose between cout/cerr and logging
    // system.
    class DPLogger
    {
    public:
      static bool useLogger;
    };
  }
}

// Log an informational message.
#define DPLOG_INFO_STR(stream)       \
  if (DPLogger::useLogger) {     \
    LOG_INFO_STR (stream);       \
  } else {                       \
    std::cout << stream << endl; \
  }

// Log a fatal message.
#define DPLOG_WARN_STR(stream)       \
  if (DPLogger::useLogger) {     \
    LOG_WARN_STR (stream);       \
  } else {                       \
    std::cerr << stream << endl; \
  }

// Log an informational message.
#define DPLOG_INFO(msg, removeEndl)    \
  std::string str(msg);                \
  if (removeEndl  &&  str.size() > 0  &&  str[str.size()-1] == '\n') { \
    str = str.substr(0, str.size()-1); \
  }                                    \
  if (DPLogger::useLogger) {           \
    LOG_INFO (str);                    \
  } else {                             \
    std::cout << str << endl;          \
  }

#endif
