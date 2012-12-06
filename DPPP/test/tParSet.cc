//# tParSet.cc: Test for class ParSet
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
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>

using namespace LOFAR;
using namespace LOFAR::DPPP;

void doTest()
{
  ParameterSet parset;
  parset.add ("key1", "abc");
  parset.add ("key2", "def");
  parset.add ("key3", "g");
  ParSet pset(parset);
  ASSERT (pset.unusedKeys().size() == 3);
  ASSERT (pset.getString("key1") == "abc");
  ASSERT (pset.getString("key1", "") == "abc");
  ASSERT (pset.getString("key1a", "12") == "12");
  ASSERT (pset.getString("key3") == "g");
  vector<string> unused = pset.unusedKeys();
  ASSERT (unused.size() == 1);
  ASSERT (unused[0] == "key2");
}

int main()
{
  try {
    doTest();
  } catch (std::exception& x) {
    cout << "Unexpected exception: " << x.what() << endl;
    return 1;
  }
  return 0;
}
