//# tMirror.cc: Test if the way of mirroring done in MedFlagger is fine
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
#include <Common/LofarTypes.h>
#include <Common/LofarLogger.h>
#include <iostream>

void doChan (int windowSize, int nchan, int chan)
{
  // At the beginning or end of the window the values are wrapped.
  // So we might need to move in two parts.
  int hw = windowSize/2;
  int s1 = chan - hw;
  int e1 = chan + hw + 1;
  int s2 = 1;
  int e2 = 1;
  if (s1 < 0) {
    e2 = -s1 + 1;
    s1 = 0;
  } else if (e1 > nchan) {
    s2 = nchan + nchan - e1 - 1; // e1-nchan+1 too far, so go back that amount
    e2 = nchan-1;
    e1 = nchan;
  }
  std::cout <<"wdw,nch=" << windowSize << ',' << nchan << " chan=" << chan
            << ' ' << s1 << '-' << e1 << ' ' << s2 << '-' << e2 << std::endl;
  ASSERT (e1-s1 + e2-s2 == windowSize);
}


int main (int argc, char* argv[])
{
  uint windowSize = 5;
  uint nchan = 8;
  if (argc > 1) {
    windowSize = atoi(argv[1]);
  }
  if (argc > 2) {
    nchan = atoi(argv[2]);
  }
  if (windowSize == 0) windowSize = 1;
  if (windowSize > nchan) windowSize = nchan;
  if (windowSize%2 == 0) windowSize--;

  for (uint i=0; i<nchan; ++i) {
    doChan (windowSize, nchan, i);
  }
}
