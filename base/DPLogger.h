// DPLogger.h: Log on cout/cerr or through the logging system
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Log on cout/cerr or through the logging system
/// @author Ger van Diepen

#ifndef DPPP_DPLOGGER_H
#define DPPP_DPLOGGER_H

#include <iostream>

namespace dp3 {
namespace base {

/// @brief This class contains the flag to choose between cout/cerr and logging
/// system.
class DPLogger {
 public:
  static bool useLogger;
};
}  // namespace base
}  // namespace dp3

/// Log an informational message.
#define DPLOG_INFO_STR(stream) std::cout << stream << std::endl;

/// Log a warning message.
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
