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
#include <Common/StreamUtil.h>
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
                                    int64 ntimes) const
    {
      int64 npoints = ntimes * itsChanCounts.size();
      os << endl << "Percentage of points flagged per baseline"
         " (antenna pair):" << endl;
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
      // Determine nr of antennae used.
      int nrused = 0;
      for (uint i=0; i<nrant; ++i) {
        if (nusedAnt[i] > 0) {
          nrused++;
        }
      }
      // Print 15 antennae per line.
      const int nantpl = 15;
      int nrl = (nrused + nantpl - 1) / nantpl;
      int ant = 0;
      int ia = ant;
      // Loop over nr of lines needed for the antennae.
      for (int i=0; i<nrl; ++i) {
        // Determine nrant per line
        int nra = std::min(nantpl, nrused - i*nantpl);
        // Print the header for the antennae being used.
        os << endl << " ant";
        for (int j=0; j<nra;) {
          if (nusedAnt[ia] > 0) {
            os << std::setw(5) << ia;
            j++;
          }
          ia++;
        }
        os << endl;
        // Print the percentages per antenna pair.
        for (uint k=0; k<nrant; ++k) {
          if (nusedAnt[k] > 0) {
            os << std::setw(4) << k << " ";
            ia = ant;
            for (int j=0; j<nra;) {
              if (nusedAnt[ia] > 0) {
                if (nusedBL(k,ia) > 0) {
                  os << std::setw(4)
                     << int((100. * countBL(k,ia)) /
                            (nusedBL(k,ia) * npoints) + 0.5)
                     << '%';
                } else {
                  os << "     ";
                }
                j++;
              }
              ia++;
            }
            os << endl;
          }
        }
        // Print the percentages per antenna.
        os << "TOTAL";
        ia = ant;
        for (int j=0; j<nra;) {
          if (nusedAnt[ia] > 0) {
            os << std::setw(4)
               << int((100. * countAnt[ia]) /
                      (nusedAnt[ia] * npoints) + 0.5) << '%';
            j++;
          }
          ia++;
        }
        os << endl;
      }
      ant = ia;
    }

    void FlagCounter::showChannel (ostream& os, int64 ntimes) const
    {
      int64 npoints = ntimes * itsBLCounts.size();
      os << endl << "Percentage of points flagged per channel:" << endl;
      if (npoints == 0) {
        return;
      }
      // Print 10 channels per line.
      const int nchpl = 10;
      os << " channels    ";
      for (int i=0; i<std::min(nchpl, int(itsChanCounts.size())); ++i) {
        os << std::setw(5) << i;
      }
      os << endl;
      int nrl = (itsChanCounts.size() + nchpl - 1) / nchpl;
      int ch = 0;
      for (int i=0; i<nrl; ++i) {
        int nrc = std::min(nchpl, int(itsChanCounts.size() - i*nchpl));
        os << std::setw(4) << ch << '-' << std::setw(4) << ch+nrc-1 << ":    ";
        for (int j=0; j<nrc; ++j) {
          os << std::setw(4) << int((100. * itsChanCounts[ch]) / npoints + 0.5)
             << '%';
          ch++;
        }
        os << endl;
      }
    }

    void FlagCounter::showCorrelation (ostream& os, int64 ntimes) const
    {
      int64 ntotal = ntimes * itsBLCounts.size() * itsChanCounts.size();
      os << endl << "Percentage of points flagged per correlation:" << endl;
      os << "  " << itsCorrCounts << " out of " << ntotal << " points   [";
      for (uint i=0; i<itsCorrCounts.size(); ++i) {
        if (i > 0) {
          os << ", ";
        }
        os << int(100. * itsCorrCounts[i] / ntotal + 0.5) << '%';
      }
      os << ']' << endl;
    }

    void FlagCounter::showPerc1 (ostream& os, double value, double total)
    {
      int perc = int(1000. * value / total + 0.5);
      os << std::setw(3) << perc/10 << '.' << perc%10 << '%';
    }


  } //# end namespace
}
