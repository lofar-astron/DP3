//# tparse.cc: Test program for function PreFlagger::PSet::exprToRpn
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
#include <DPPP/PreFlagger.h>
#include <Common/StringUtil.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <iostream>

using namespace LOFAR;
using namespace LOFAR::DPPP;
using namespace std;

namespace LOFAR {
  namespace DPPP {
    // This class name should match the friend in PreFlagger.
    class TestPSet
    {
    public:
      static void testParse (const string&);
    };
  }
}

void TestPSet::testParse (const string& expr)
{
  PreFlagger::PSet pset;
  vector<string> names = pset.exprToRpn (expr);
  cout << pset.itsRpn << endl;
  cout << names << endl;
}


int main(int argc, char* argv[])
{
  INIT_LOGGER ("tPSet");
  try {
    if (argc > 1) {
      TestPSet::testParse (argv[1]);
    } else {
      TestPSet::testParse ("(s1&s_1)|!(!!s2&s2)");
    }
  } catch (std::exception& x) {
    cout << x.what() << endl;
    return 1;
  }
  return 0;
}
