// DP3.cc: Program to execute steps like averaging and flagging on an MS
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include <aocommon/logger.h>
#include <H5Cpp.h>

#include <dp3/base/DP3.h>

int main(int argc, char* argv[]) {
  try {
    dp3::base::ExecuteFromCommandLine(argc, argv);
    return 0;
  } catch (const std::exception& err) {
    aocommon::Logger::Error << "\nstd exception detected: " << err.what()
                            << '\n';
    return 1;
  } catch (const H5::Exception& err) {
    // Since H5::Exception is not derived from std::exception, DP3 crashes
    // if it does not catch them -> Print an error message instead.
    aocommon::Logger::Error << "\nH5 exception detected: " << err.getDetailMsg()
                            << "\nStack trace:\n";
    err.printErrorStack();
    aocommon::Logger::Error << '\n';
    return 1;
  }
}
