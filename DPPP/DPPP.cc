//# DPPP.cc: Program to execute steps like averaging and flagging on an MS
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

#include "DPRun.h"
#include "Version.h"

#include <boost/filesystem/operations.hpp>

#include <iostream>
#include <stdexcept>

// Define handler that tries to print a backtrace.
//Exception::TerminateHandler t(Exception::terminate);

void showUsage() {
  std::cout <<
    "Usage: DPPP [-v] [parsetfile] [parsetkeys...]\n"
    "  parsetfile: a file containing one parset key=value pair per line\n"
    "  parsetkeys: any number of parset key=value pairs, e.g. msin=my.MS\n\n"
    "If both a file and command-line keys are specified, the keys on the command\n"
    "line override those in the file.\n"
    "If no arguments are specified, the program tries to read \"NDPPP.parset\"\n"
    "or \"DPPP.parset\" as a default.\n"
    "-v will show version info and exit.\n"
    "Documentation is at:\n"
    "http://www.lofar.org/wiki/doku.php?id=public:user_software:documentation:ndppp\n";
}

int main(int argc, char *argv[])
{
  try
  {
    // Get the name of the parset file.
    if (argc>1) {
      string param = argv[1];
      if(param=="--help" ||
        param=="-help" || param=="-h" || 
        param=="--usage" || param=="-usage") {
        showUsage();
        return 0;
      }
      else if(param == "-v" || param == "--version")
      {
        std::cout << DPPPVersion::AsString() << '\n';
        return 0;
      }
    }

    string parsetName;
    if (argc > 1  &&  string(argv[1]).find('=') == string::npos) {
      // First argument is parset name (except if it's a key-value pair)
      parsetName = argv[1];
    } else if (argc==1) {
      // No arguments given: try to load [N]DPPP.parset
      if (boost::filesystem::exists("DPPP.parset")) {
        parsetName="DPPP.parset";
      } else if (boost::filesystem::exists("NDPPP.parset")) {
        parsetName="NDPPP.parset";
      } else { // No default file, show usage and exit
        showUsage();
        return 0;
      }
    }

    // Execute the parset file.
    DP3::DPPP::DPRun::execute (parsetName, argc, argv);
  } catch (std::exception& err) {
    std::cerr << "\nstd exception detected: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
