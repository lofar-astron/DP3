//# showsourcedb.cc: Show contents of a SourceDB catalog
//#
//# Copyright (C) 2013
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
//# $Id: mergesourcedb.cc 24953 2013-05-17 11:34:36Z diepen $

// This program shows the contents of a SourceDB catalogs.
// It can show only patches or sources or both.
//
// The program can be run as:
//    showsourcedb  in=inname mode=all|patch|source
// in            name of the SourceDB.
// mode          all    = show patches and sources
//               patch  = only show patches
//               source = only show sources

#include "SourceDB.h"

#include "../Common/StringUtil.h"
#include "../Common/StreamUtil.h"

#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Inputs/Input.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <vector>

using namespace std;
using namespace casacore;

void show (const string& name, const string& mode, const string& patt)
{
  // Open the input SourceDB.
  DP3::BBS::SourceDB in ((DP3::BBS::ParmDBMeta(string(), name)));
  // Read all patches from the SourceDB and write them.
  vector<DP3::BBS::PatchInfo> patch (in.getPatchInfo(-1, patt));
  for (size_t i=0; i<patch.size(); ++i)
  {
    if (mode != "source")
    {
      cout << patch[i] << '\n';
    }
    if (mode != "patch")
    {
      vector<DP3::BBS::SourceData> sources(in.getPatchSourceData (patch[i].getName()));
      for (DP3::BBS::SourceData& s : sources)
        s.print (cout);
    }
  }
}


int main (int argc, char* argv[])
{
  try {
    // Define the input parameters.
    Input inputs(1);
    inputs.version ("GvD 2013-Jun-12");
    inputs.create("in", "",
                  "Input SourceDB", "string");
    inputs.create ("mode", "all",
                   "patch=show all patches, "
                   "source=show all sources, "
                   "all=show patches and sources",
                   "string");
    inputs.create ("patches", "*",
                   "Pattern for names of patches to show", "string");
    // Read and check the input parameters.
    inputs.readArguments(argc, argv);
    string in = inputs.getString("in");
    if(in.empty())
      throw std::runtime_error("no input sourcedb name given");
    string mode = boost::to_lower_copy(inputs.getString("mode"));
    if(mode!="patch" && mode!="source" && mode!="all")
      throw std::runtime_error("incorrect mode given");
    string patt = inputs.getString("patches");
    show (in, mode, patt);
  } catch (AipsError& x) {
    cerr << "Caught AIPS error: " << x.what() << endl;
    return 1;
  }
  catch (std::exception& x) {
    cerr << "Caught LOFAR exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
