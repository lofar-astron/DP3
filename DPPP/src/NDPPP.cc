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

#include <lofar_config.h>
#include <DPPP/DPRun.h>
#include <DPPP/Package__Version.h>
#include <Common/LofarLogger.h>
#include <Common/SystemUtil.h>
#include <Common/Exception.h>
#include <iostream>
#include <stdexcept>

#include <casacore/casa/OS/File.h>

using namespace LOFAR::DPPP;
using namespace LOFAR;

// Define handler that tries to print a backtrace.
Exception::TerminateHandler t(Exception::terminate);

void showUsage() {
  std::cout<<"Usage: DPPP [-v] [parsetfile] [parsetkeys...]"<<std::endl;
  std::cout<<"  parsetfile: a file containing one parset key=value pair "<<
    "per line"<<std::endl;
  std::cout<<"  parsetkeys: any number of parset key=value pairs, e.g. "<<
    "msin=my.MS"<<std::endl<<std::endl;
  std::cout<<"If both a file and command-line keys are specified, the "<<
    "keys on the command line override those in the file."<<std::endl;
  std::cout<<"If no arguments are specified, the program tries to read "<<
    "\"NDPPP.parset\" or \"DPPP.parset\" as a default."<<std::endl;
  std::cout<<"-v will show version info and exit."<<std::endl;
  std::cout<<"Documentation is at http://www.lofar.org/wiki/doku.php?id="<<
    "public:user_software:documentation:ndppp"<<std::endl;
}

int main(int argc, char *argv[])
{
  try
  {
    TEST_SHOW_VERSION (argc, argv, DPPP);
    INIT_LOGGER("DPPP");
    // Get the name of the parset file.
    if (argc>1 && (
          string(argv[1])=="--help" ||
          string(argv[1])=="-help"   || string(argv[1])=="-h" || 
          string(argv[1])=="--usage" || string(argv[1])=="-usage")) {
      showUsage();
      return 0;
    }

    string parsetName;
    if (argc > 1  &&  string(argv[1]).find('=') == string::npos) {
      // First argument is parset name (except if it's a key-value pair)
      parsetName = argv[1];
    } else if (argc==1) {
      // No arguments given: try to load [N]DPPP.parset
      if (casacore::File("NDPPP.parset").exists()) {
        parsetName="NDPPP.parset";
      } else if (casacore::File("DPPP.parset").exists()) {
        parsetName="DPPP.parset";
      } else { // No default file, show usage and exit
        showUsage();
        return 0;
      }
    }

    // Execute the parset file.
    DPRun::execute (parsetName, argc, argv);
  } catch (LOFAR::APSException& err) {
    // just send err.what() to the error stream
    // this is just the error message, not a full backtrace
    std::cerr << std::endl;
    std::cerr << "ParameterSet Exception detected: "<< err.what() << std::endl;
    return 1;
  } catch (LOFAR::Exception& err) {
    std::cerr << "LOFAR Exception detected: " << err << std::endl;
    return 1;
#ifdef __clang__
  } catch (std::exception& err) {
    std::cerr << std::endl;
    std::cerr << "std exception detected: " << err.what() << std::endl;
    return 1;
#endif
  }
  return 0;
}
