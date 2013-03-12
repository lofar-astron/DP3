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

#include <casa/Arrays/Matrix.h>
#include <casa/Utilities/LinearSearch.h>
#include <Common/StreamUtil.h>
#include <Common/lofar_vector.h>

using namespace casa;
using namespace std;
using namespace LOFAR;

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

void testAdd()
{
  uint nrold = 5;
  uint nrnew = 7;
  vector<int> itsAnt1(15);
  vector<int> itsAnt2(15);
  int inx=0;
  for (int a1=0; a1<5; ++a1) {
    for (int a2=a1; a2<5; ++a2) {
      itsAnt1[inx] = a1;
      itsAnt2[inx] = a2;
      ++inx;
    }
  }
  vector<Vector<int> > itsParts(2);
  itsParts[0].resize (2);
  itsParts[0][0] = 0;
  itsParts[0][1] = 1;
  itsParts[1].resize (2);
  itsParts[1][0] = 3;
  itsParts[1][1] = 4;
  
  
      vector<int> newbl(nrnew);
      vector<vector<int> > itsBufRows;
      bool itsMakeAutoCorr = true;
      // Loop over the superstations.
      // Note that by making this the outer loop, the baselines between
      // superstations are also formed.
      // At the end the new baselines are added to itsAnt1 and itsAnt2.
      // itsBufRows contains for each new baseline the rownrs in the DPBuffer
      // to be added for the new baseline. If rownr<0, the conjugate has to be
      // added (1 is added to rownr, otherwise 0 is ambiguous).
      // Note that a rownr can be the rownr of a new baseline.
      for (uint j=0; j<itsParts.size(); ++j) {
        std::fill (newbl.begin(), newbl.end(), -1);
        vector<int> newAnt1;
        vector<int> newAnt2;
        // Loop through all baselines and find out if a baseline should
        // be used for a superstation.
        for (uint i=0; i<itsAnt1.size(); ++i) {
          bool havea1 = linearSearch1 (itsParts[j], itsAnt1[i]) >= 0;
          bool havea2 = linearSearch1 (itsParts[j], itsAnt2[i]) >= 0;
          int  ant    = nrold+j;
          int  take   = 0;
          if (havea1) {
            // If both stations are in same superstation, only use them
            // if it is an autocorrelation.
            if (havea2) {
              if (itsMakeAutoCorr  &&  itsAnt1[i] == itsAnt2[i]) {
                take = 1;
              }
            } else {
              ant  = itsAnt2[i];
              take = -1;            // conjugate has to be added
            }
          } else if (havea2) {
            ant  = itsAnt1[i];
            take = 1;
          }
          if (take != 0) {
            // We have a baseline for the superstation.
            // Get its index; create it if not used before.
            int blinx = newbl[ant];
            if (blinx < 0) {
              blinx = newbl[ant] = itsBufRows.size();
              itsBufRows.push_back (vector<int>());
              newAnt1.push_back (ant);
              newAnt2.push_back (nrold+j);
            }
            itsBufRows[blinx].push_back (take*(i+1));
          }
        }
        // Copy the new baselines for this superstation to the baseline list.
        // Give a warning if nothing found.
        if (newAnt1.empty()) {
          //          DPLOG_WARN_STR ("StationAdder: no baseline found for superstation");
          cout << "StationAdder: no baseline found for superstation" << endl;
        } else {
          uint oldsz = itsAnt1.size();
          itsAnt1.resize (oldsz + newAnt1.size());
          itsAnt2.resize (oldsz + newAnt1.size());
          for (uint i=0; i<newAnt1.size(); ++i) {
            itsAnt1[oldsz+i] = newAnt1[i];
            itsAnt2[oldsz+i] = newAnt2[i];
          }
        }
      }
      cout << itsAnt1<<endl<<itsAnt2<<endl;
      writeVector (cout, itsBufRows);
      cout<<endl;
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

  testAdd();
}
