//# tMedian.cc: Program to test performance of kthLargest and nth_element
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
#include <casa/OS/Timer.h>
#include <casa/Utilities/GenSort.h>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace casa;
using namespace std;

void testCasa (size_t sz, size_t n)
{
  Timer timer;
  for (uint i=0; i<n; ++i) {
    vector<float> vec(sz);
    uint j=0;
    for (vector<float>::iterator iter=vec.begin(); iter!=vec.end(); ++iter) {
      *iter = j++;
    }
    GenSort<float>::kthLargest (&(vec[0]), sz, sz/2);
  }
  timer.show ("casa");
}

void testStl (size_t sz, size_t n)
{
  Timer timer;
  for (uint i=0; i<n; ++i) {
    vector<float> vec(sz);
    uint j=0;
    for (vector<float>::iterator iter=vec.begin(); iter!=vec.end(); ++iter) {
      *iter = j++;
    }
    nth_element (vec.begin(), vec.end(), vec.begin()+sz/2);
  }
  timer.show ("stl ");
}

int main (int argc, char* argv[])
{
  size_t window = 155;
  size_t ntimes = 1000;
  if (argc > 1) window = atoi(argv[1]);
  if (argc > 2) ntimes = atoi(argv[2]);
  testCasa (window, ntimes);
  testStl  (window, ntimes);
  return 0;
}
