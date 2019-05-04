//# Simulate.cc: Simulate visibilities for a patch of sources.
//#
//# Copyright (C) 2012
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

#include "Simulate.h"
#include "PhaseShift.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/MatrixMath.h>

using namespace casacore;

namespace DP3
{
namespace DPPP
{

vector<int> nsetupSplitUVW (unsigned int nant, const Vector<int>& ant1,
                           const Vector<int>& ant2)
{
  // Get the indices of the baselines needed to split the baseline UVWs into
  // station UVWs. They are in such an order that the UVW of a station is known
  // before used in another baseline to derive the UVW of the other station.
  // It can handle cases where baselines occur in disjoint station groups
  // like 0-1, 0-2, 1-2 and 3-4, 4-5, 5-6.
  // Note that the first station of a group gets UVW=0. All other station UVWs
  // are relative to it using the baseline UVWs.
  // Also note that nr of groups can be derived from the size of the returned
  // vector (because it contains no entry for the first antenna in a group).
  vector<int> uvwbl;
  uvwbl.reserve (nant);
  Block<bool> known(nant, false);
  unsigned int nset = 0;
  // Loop until all stations are set.
  while (nset < nant) {
    // Disjoint groups might exist, so keep a vector containing related antennae
    // which are members of the same group.
    vector<unsigned int> members(1);
    // Set first unset station as the reference station (which gets UVW=0).
    for (unsigned int i=0; i<nant; ++i) {
      if (!known[i]) {
        members[0] = i;
        known[i] = true;
        ++nset;
        break;
      }
    }
    // Loop through all members in the group.
    // Note that new members can be appended in this loop.
    for (unsigned int j=0; j<members.size(); ++j) {
      int refst = members[j];
      // Find all stations having a baseline with the reference station.
      for (unsigned int i=0; i<ant1.size(); ++i) {
        int a1 = ant1[i];
        int a2 = ant2[i];
        // Only take baselines into account for which one station is known,
        // so the other can be derived from it.
        // The unknown station becomes member of the group.
        if (known[a1] != known[a2]) {
          if (a1 == refst) {
            uvwbl.push_back (i);
            members.push_back (a2);
            known[a2] = true;
            ++nset;
          } else if (a2 == refst) {
            uvwbl.push_back (-(i+1));
            members.push_back (a1);
            known[a1] = true;
            ++nset;
          }
        }
      }
    }
  }
  return uvwbl;
}

void nsplitUVW (const vector<int>& blindex,
               const vector<Baseline>& baselines,
               const Matrix<double>& uvwbl,
               Matrix<double>& uvwant)
{
  uvwant = 0.;
  double* uvwl;
  double* uvwr;
  const double* uvwb;
  for (unsigned int i=0; i<blindex.size(); ++i) {
    int inx = blindex[i];
    if (inx < 0) {
      // Ant2 is known.
      inx = -inx - 1;
      //ASSERT((unsigned int)(uvwant.shape()[1])>baselines[inx].first);
      //ASSERT((unsigned int)(uvwant.shape()[1])>baselines[inx].second);
      uvwl = uvwant.data() + 3 * baselines[inx].first;
      uvwr = uvwant.data() + 3 * baselines[inx].second;
      uvwb = uvwbl.data() + 3*inx;
      for (int j=0; j<3; ++j) {
        uvwl[j] = uvwr[j] - uvwb[j];
      }
    } else {
      //ASSERT((unsigned int)(uvwant.shape()[1])>baselines[inx].first);
      //ASSERT((unsigned int)(uvwant.shape()[1])>baselines[inx].second);
      // Ant1 is known.
      uvwl = uvwant.data() + 3 * baselines[inx].first;
      uvwr = uvwant.data() + 3 * baselines[inx].second;
      uvwb = uvwbl.data() + 3*inx;
      for (int j=0; j<3; ++j) {
        uvwr[j] = uvwl[j] + uvwb[j];
      }
    }
  }
}


void rotateUVW(const Position &from, const Position &to, size_t nUVW,
    double *uvw)
{
    casacore::Matrix<double> oldUVW(3,3);
    casacore::Matrix<double> newUVW(3,3);
    PhaseShift::fillTransMatrix(oldUVW, from[0], from[1]);
    PhaseShift::fillTransMatrix(newUVW, to[0], to[1]);

    casacore::Matrix<double> tmp(casacore::product(casacore::transpose(newUVW), oldUVW));
    const double *R = tmp.data();

    for(size_t i = 0; i < 3*nUVW; i+=3)
    {
        // Compute rotated UVW.
        double u = uvw[i+0] * R[0] + uvw[i+1] * R[3] + uvw[i+2] * R[6];
        double v = uvw[i+0] * R[1] + uvw[i+1] * R[4] + uvw[i+2] * R[7];
        double w = uvw[i+0] * R[2] + uvw[i+1] * R[5] + uvw[i+2] * R[8];

        uvw[i+0] = u;
        uvw[i+1] = v;
        uvw[i+2] = w;

        // Move to next station.
    } // Stations.
}

} //# namespace DPPP
} //# namespace LOFAR
