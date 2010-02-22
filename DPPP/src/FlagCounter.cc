//# FlagCounter.cc: Class to keep counts of nr of flagged points
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
#include <DPPP/FlagCounter.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    void FlagCounter::init (uint nbaselines, uint nchan, uint ncorr)
    {
      itsBLCounts.resize (nbaselines);
      itsChanCounts.resize (nchan);
      itsCorrCounts.resize (ncorr);
      std::fill (itsBLCounts.begin(), itsBLCounts.end(), 0);
      std::fill (itsChanCounts.begin(),itsChanCounts.end(), 0);
      std::fill (itsCorrCounts.begin(),itsCorrCounts.end(), 0);
    }

    void FlagCounter::showBaseline (ostream& os, const casa::Vector<int>& ant1,
                                    const casa::Vector<int>& ant2,
                                    int64 npoints) const
    {
      os << endl << "Percentage of points flagged per baseline:" << endl;
      uint nrant = 1 + std::max(max(ant1), max(ant2));
      // Collect counts per baseline and antenna.
      Vector<int64> nusedAnt(nrant, 0);
      Vector<int64> countAnt(nrant, 0);
      Matrix<int64> nusedBL (nrant, nrant, 0);
      Matrix<int64> countBL (nrant, nrant, 0);
      for (uint i=0; i<itsBLCounts.size(); ++i) {
        countBL(ant1[i], ant2[i]) += itsBLCounts[i];
        countBL(ant2[i], ant1[i]) += itsBLCounts[i];
        nusedBL(ant1[i], ant2[i])++;
        nusedBL(ant2[i], ant1[i])++;
        countAnt[ant1[i]] += itsBLCounts[i];
        countAnt[ant2[i]] += itsBLCounts[i];
        nusedAnt[ant1[i]]++;
        nusedAnt[ant2[i]]++;
      }
      // Print the header for the antennae being used.
      for (uint i=0; i<nrant; ++i) {
        if (nusedAnt[i] > 0) {
          os << std::setw(5) << i;
        }
      }
      // Print the percentages per antenna pair.
      for (uint i=0; i<nrant; ++i) {
        if (nusedAnt[i] > 0) {
          os << std::setw(4) << i << "  ";
          for (uint j=0; j<nrant; ++j) {
            if (nusedBL(j,i) > 0) {
              os << std::setw(4)
                 << int(100. * countBL(j,i) / (nusedBL(j,i) * npoints)) << '%';
            } else {
              os << "     ";
            }
          }
        }
      }
      // Print the percentages per antenna.
      os << "TOTAL";
      for (uint i=0; i<nrant; ++i) {
        if (nusedAnt[i] > 0) {
          os << std::setw(4)
             << int(100. * countAnt[i] / (nusedAnt[i] * npoints)) << '%';
        } else {
          os << "     ";
        }
      }
    }

    void FlagCounter::showChannel (ostream& os, int64 npoints) const
    {
      os << endl << "Percentage of points flagged per channel:" << endl;
      for (uint i=0; i<itsChanCounts.size(); ++i) {
        os << "Channel " << std::setw(4) << i << ":   "
           << std::setw(4) << int(100. * itsChanCounts[i] / npoints)
           << '%' << endl;
      }
    }

    void FlagCounter::showCorrelation (ostream& os) const
    {
      os << endl << "Points flagged by correlation:" << endl;
      for (uint i=0; i<itsCorrCounts.size(); ++i) {
        os << "Correlation " << i << ": " << itsCorrCounts[i] << endl;
      }
    }

  } //# end namespace
}
