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

#include <lofar_config.h>
#include <DPPP/Simulate.h>
#include <DPPP/Simulator.h>
#include <Common/LofarLogger.h>

// Only required for rotateUVW().
#include <DPPP/PhaseShift.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixMath.h>

using namespace casa;

namespace LOFAR
{
namespace DPPP
{

vector<int> nsetupSplitUVW (uint nant, const Vector<int>& ant1,
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
  uint nset = 0;
  // Loop until all stations are set.
  while (nset < nant) {
    // Disjoint groups might exist, so keep a vector containing related antennae
    // which are members of the same group.
    vector<uint> members(1);
    // Set first unset station as the reference station (which gets UVW=0).
    for (uint i=0; i<nant; ++i) {
      if (!known[i]) {
        members[0] = i;
        known[i] = true;
        ++nset;
        break;
      }
    }
    // Loop through all members in the group.
    // Note that new members can be appended in this loop.
    for (uint j=0; j<members.size(); ++j) {
      int refst = members[j];
      // Find all stations having a baseline with the reference station.
      for (uint i=0; i<ant1.size(); ++i) {
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
  for (uint i=0; i<blindex.size(); ++i) {
    int inx = blindex[i];
    if (inx < 0) {
      // Ant2 is known.
      inx = -inx - 1;
      uvwl = uvwant.data() + 3 * baselines[inx].first;
      uvwr = uvwant.data() + 3 * baselines[inx].second;
      uvwb = uvwbl.data() + 3*inx;
      for (int j=0; j<3; ++j) {
        uvwl[j] = uvwr[j] - uvwb[j];
      }
    } else {
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

void splitUVW(size_t nStation, size_t nBaseline,
    const_cursor<Baseline> baselines, const_cursor<double> uvw,
    cursor<double> split)
{
    // If the number of stations is zero, then the number of baselines should be
    // zero as well because no valid baselines exist in this case.
    ASSERT(nStation > 0 || nBaseline == 0);

    // Flags to keep track of the stations for which the UVW coordinates are
    // known.
    vector<bool> flag(nStation, true);

    // Find reachable stations, i.e. stations that participate in at least one
    // (cross) baseline.
    size_t nReachable = 0;
    const_cursor<Baseline> tmp(baselines);
    for(size_t i = 0; i < nBaseline; ++i)
    {
        const size_t p = tmp->first;
        const size_t q = tmp->second;

        if(p != q)
        {
            if(flag[p])
            {
                flag[p] = false;
                ++nReachable;
            }

            if(flag[q])
            {
                flag[q] = false;
                ++nReachable;
            }
        }

        if(nReachable == nStation)
        {
            break;
        }

        // Move to next baseline.
        ++tmp;
    }

    // Zero the UVW coordinates of all unreachable stations and flag them as
    // known.
    if(nReachable < nStation)
    {
        for(size_t i = 0; i < nStation; ++i)
        {
            if(flag[i])
            {
                split.forward(1, i);
                split[0] = 0.0;
                split[1] = 0.0;
                split[2] = 0.0;
                split.backward(1, i);
            }
        }
    }

    // If no stations are reachable, nothing more can be done.
    if(nReachable == 0)
    {
        return;
    }

    size_t ref = 0;
    cursor<double> known(split);
    while(true)
    {
        // Find first station for which UVW coordinates are unknown.
        while(ref < nStation && flag[ref])
        {
            ++ref;
        }

        // UVW coordinates known for all stations, done.
        if(ref == nStation)
        {
            break;
        }

        // Set UVW coordinates of the reference station to zero. Note that each
        // isolated group of stations will have such a reference station. This
        // is OK because the groups are isolated.
        split.forward(1, ref);
        split[0] = 0.0;
        split[1] = 0.0;
        split[2] = 0.0;
        split.backward(1, ref);
        flag[ref] = true;

        // Single isolated station left, no use iterating over baselines.
        if(ref == nStation - 1)
        {
            break;
        }

        bool done = false;
        while(!done)
        {
            // Find UVW coordinates for other stations linked to the reference
            // station via baselines.
            size_t nFound = 0;
            for(size_t i = 0; i < nBaseline; ++i)
            {
                const size_t p = baselines->first;
                const size_t q = baselines->second;
                if(p != q && flag[p] != flag[q])
                {
                    if(flag[p])
                    {
                        known.forward(1, p);
                        split.forward(1, q);
                        split[0] = uvw[0] + known[0];
                        split[1] = uvw[1] + known[1];
                        split[2] = uvw[2] + known[2];
                        split.backward(1, q);
                        known.backward(1, p);
                        flag[q] = true;
                    }
                    else
                    {
                        known.forward(1, q);
                        split.forward(1, p);
                        split[0] = -uvw[0] + known[0];
                        split[1] = -uvw[1] + known[1];
                        split[2] = -uvw[2] + known[2];
                        split.backward(1, p);
                        known.backward(1, q);
                        flag[p] = true;
                    }

                    ++nFound;
                }

                // Move to next baseline.
                uvw.forward(1);
                ++baselines;
            } // Baselines.

            // Reset cursors.
            uvw.backward(1, nBaseline);
            baselines -= nBaseline;

            // Depending on the baseline order it may be possible to derive the
            // UVW coordinates of additional stations when UVW coordinates were
            // found for at leat one station (i.e. nFound larger than 0).
            // Another iteration is required in this case, unless UVW
            // coordinates were found for all (remaining) stations (i.e. nFound
            // equals nStation - ref - 1).
            done = (nFound == 0 || nFound == nStation - ref - 1);
        }
    }

    ASSERT(static_cast<size_t>(std::count(flag.begin(), flag.end(), true))
        == flag.size());
}

void rotateUVW(const Position &from, const Position &to, size_t nUVW,
    cursor<double> uvw)
{
    casa::Matrix<double> oldUVW(3,3);
    casa::Matrix<double> newUVW(3,3);
    PhaseShift::fillTransMatrix(oldUVW, from[0], from[1]);
    PhaseShift::fillTransMatrix(newUVW, to[0], to[1]);

    casa::Matrix<double> tmp(casa::product(casa::transpose(newUVW), oldUVW));
    const double *R = tmp.data();

    for(size_t i = 0; i < nUVW; ++i)
    {
        // Compute rotated UVW.
        double u = uvw[0] * R[0] + uvw[1] * R[3] + uvw[2] * R[6];
        double v = uvw[0] * R[1] + uvw[1] * R[4] + uvw[2] * R[7];
        double w = uvw[0] * R[2] + uvw[1] * R[5] + uvw[2] * R[8];

        uvw[0] = u;
        uvw[1] = v;
        uvw[2] = w;

        // Move to next station.
        uvw.forward(1);
    } // Stations.
}

void simulate(const Position &reference, const Patch::ConstPtr &patch,
    size_t nStation, size_t nBaseline, size_t nChannel,
    const_cursor<Baseline> baselines, const_cursor<double> freq,
    const_cursor<double> uvw, cursor<dcomplex> vis)
{
    Simulator simulator(reference, nStation, nBaseline, nChannel, baselines,
        freq, uvw, vis);
    for(size_t i = 0; i < patch->nComponents(); ++i)
    {
      simulator.simulate(patch->component(i));
    }
}

} //# namespace DPPP
} //# namespace LOFAR
