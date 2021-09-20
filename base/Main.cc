// DPPP.cc: Program to execute steps like averaging and flagging on an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DP3.h"

#include <Version.h>

#include <aocommon/checkblas.h>

#include <boost/filesystem/operations.hpp>

#include <iostream>
#include <stdexcept>

// Define handler that tries to print a backtrace.
// Exception::TerminateHandler t(Exception::terminate);

void showUsage() {
  std::cout
      << "Usage: DP3 [-v] [parsetfile] [parsetkeys...]\n"
         "  parsetfile: a file containing one parset key=value pair per line\n"
         "  parsetkeys: any number of parset key=value pairs, e.g. "
         "msin=my.MS\n\n"
         "If both a file and command-line keys are specified, the keys on the "
         "command\n"
         "line override those in the file.\n"
         "If no arguments are specified, the program tries to read "
         "\"DP3.parset\",\n"
         "\"NDPPP.parset\" or \"DPPP.parset\" as a default.\n"
         "-v will show version info and exit.\n"
         "Documentation is at:\n"
         "https://www.astron.nl/citt/DP3\n";
}

int main(int argc, char* argv[]) {
  try {
    check_openblas_multithreading();

    // Get the name of the parset file.
    if (argc > 1) {
      string param = argv[1];
      if (param == "--help" || param == "-help" || param == "-h" ||
          param == "--usage" || param == "-usage") {
        showUsage();
        return 0;
      } else if (param == "-v" || param == "--version") {
        std::cout << DP3Version::AsString() << '\n';
        return 0;
      }
    }

    string parsetName;
    if (argc > 1 && string(argv[1]).find('=') == string::npos) {
      // First argument is parset name (except if it's a key-value pair)
      parsetName = argv[1];
    } else if (argc == 1) {
      // No arguments given: try to load [N]DPPP.parset
      if (boost::filesystem::exists("DP3.parset")) {
        parsetName = "DP3.parset";
      } else if (boost::filesystem::exists("DPPP.parset")) {
        parsetName = "DPPP.parset";
      } else if (boost::filesystem::exists("NDPPP.parset")) {
        parsetName = "NDPPP.parset";
      } else {  // No default file, show usage and exit
        showUsage();
        return 0;
      }
    }

    // Execute the parset file.
    dp3::base::DP3::execute(parsetName, argc, argv);
  } catch (std::exception& err) {
    std::cerr << "\nstd exception detected: " << err.what() << '\n';
    return 1;
  }
  return 0;
}
