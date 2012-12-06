//# FlagCounter.cc: Class to keep counts of nr of flagged visibilities
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
#include <DPPP/DPInput.h>
#include <Common/ParameterSet.h>
#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ScalarColumn.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayIO.h>
#include <vector>
#include <map>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    FlagCounter::FlagCounter()
      : itsInfo   (0),
        itsShowFF (false)
    {}

    FlagCounter::FlagCounter (const string& msName,
                              const ParameterSet& parset,
                              const string& prefix)
    {
      itsWarnPerc = parset.getDouble (prefix+"warnperc", 0);
      itsShowFF   = parset.getBool   (prefix+"showfullyflagged", false);
      bool save   = parset.getBool   (prefix+"save", false);
      if (save) {
        // Percentages have to be saved, so form the table name to use.
        string path = parset.getString (prefix+"path", "");
        // Use the step name (without dot) as a name suffix.
        string suffix = prefix;
        string::size_type pos = suffix.find ('.');
        if (pos != string::npos) {
          suffix = suffix.substr(0, pos);
        }
        // Use the MS name as the name.
        // If no path is given, use the path of the MS (use . if no path).
        string name = msName;
        pos = name.rfind ('/');
        if (path.empty()) {
          if (pos == string::npos) {
            path = '.';
          } else {
            path = name.substr(0, pos);
          }
        }
        name = name.substr(pos+1);
        pos = name.find ('.');
        if (pos != string::npos) {
          name = name.substr(0, pos);
        }
        itsSaveName = path + '/' + name + '_' + suffix + ".flag";
      }
    }

    void FlagCounter::init (const DPInfo& info)
    {
      itsInfo = &info;
      itsBLCounts.resize (info.nbaselines());
      itsChanCounts.resize (info.nchan());
      itsCorrCounts.resize (info.ncorr());
      std::fill (itsBLCounts.begin(), itsBLCounts.end(), 0);
      std::fill (itsChanCounts.begin(),itsChanCounts.end(), 0);
      std::fill (itsCorrCounts.begin(),itsCorrCounts.end(), 0);
    }

    void FlagCounter::add (const FlagCounter& that)
    {
      // Add that to this after checking for equal sizes.
      ASSERT (itsBLCounts.size()   == that.itsBLCounts.size());
      ASSERT (itsChanCounts.size() == that.itsChanCounts.size());
      ASSERT (itsCorrCounts.size() == that.itsCorrCounts.size());
      std::transform (itsBLCounts.begin(), itsBLCounts.end(),
                      that.itsBLCounts.begin(), itsBLCounts.begin(),
                      std::plus<int64>());
      std::transform (itsChanCounts.begin(), itsChanCounts.end(),
                      that.itsChanCounts.begin(), itsChanCounts.begin(),
                      std::plus<int64>());
      std::transform (itsCorrCounts.begin(), itsCorrCounts.end(),
                      that.itsCorrCounts.begin(), itsCorrCounts.begin(),
                      std::plus<int64>());
    }

    void FlagCounter::showBaseline (ostream& os, int64 ntimes) const
    {
      const Vector<Int>& ant1 = itsInfo->getAnt1();
      const Vector<Int>& ant2 = itsInfo->getAnt2();
      const Vector<String>& antNames = itsInfo->antennaNames();
      // Keep track of fully flagged baselines.
      std::vector<std::pair<int,int> > fullyFlagged;
      int64 npoints = ntimes * itsChanCounts.size();
      os << endl << "Percentage of visibilities flagged per baseline"
         " (antenna pair):";
      uint nrant = 1 + std::max(max(ant1), max(ant2));
      // Collect counts per baseline and antenna.
      Vector<int64> nusedAnt(nrant, 0);
      Vector<int64> countAnt(nrant, 0);
      Matrix<int64> nusedBL (nrant, nrant, 0);
      Matrix<int64> countBL (nrant, nrant, 0);
      for (uint i=0; i<itsBLCounts.size(); ++i) {
        countBL(ant1[i], ant2[i]) += itsBLCounts[i];
        nusedBL(ant1[i], ant2[i])++;
        countAnt[ant1[i]] += itsBLCounts[i];
        nusedAnt[ant1[i]]++;
        if (ant1[i] != ant2[i]) {
          countBL(ant2[i], ant1[i]) += itsBLCounts[i];
          nusedBL(ant2[i], ant1[i])++;
          countAnt[ant2[i]] += itsBLCounts[i];
          nusedAnt[ant2[i]]++;
        }
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
      // Loop over nr of lines needed for the antennae.
      for (int i=0; i<nrl; ++i) {
        int oldant = ant;
        // Determine nrant per line
        int nra = std::min(nantpl, nrused - i*nantpl);
        // Print the header for the antennae being used.
        // It also updates ant for the next iteration.
        os << endl << " ant";
        for (int j=0; j<nra;) {
          if (nusedAnt[ant] > 0) {
            os << std::setw(5) << ant;
            j++;
          }
          ant++;
        }
        os << endl;
        // Print the percentages per antenna pair.
        for (uint k=0; k<nrant; ++k) {
          if (nusedAnt[k] > 0) {
            os << std::setw(4) << k << " ";
            int ia = oldant;
            for (int j=0; j<nra;) {
              if (nusedAnt[ia] > 0) {
                if (nusedBL(k,ia) > 0) {
                  os << std::setw(4)
                     << int((100. * countBL(k,ia)) /
                            (nusedBL(k,ia) * npoints) + 0.5)
                     << '%';
                  // Determine if baseline is fully flagged.
                  // Do it only for ANT1<=ANT2
                  if (int(k) <= ia) {
                    if (countBL(k,ia) == nusedBL(k,ia) * npoints) {
                      fullyFlagged.push_back (std::pair<int,int>(k,ia));
                    }
                  }
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
        int ia = oldant;
        for (int j=0; j<nra;) {
          if (nusedAnt[ia] > 0) {
            double perc = 100. * countAnt[ia] / (nusedAnt[ia] * npoints);
            os << std::setw(4) << int(perc + 0.5) << '%';
            j++;
          }
          ia++;
        }
        os << endl;
      }
      if (itsWarnPerc > 0) {
        for (uint i=0; i<nrant; ++i) {
          if (nusedAnt[i] > 0) {
            double perc = (100. * countAnt[i]) / (nusedAnt[i] * npoints);
            if (perc >= itsWarnPerc) {
              os << "** NOTE: ";
              showPerc1 (os, perc, 100);
              os << " of data are flagged for station " << i
                 << " (" << antNames[i] << ')' << endl;
            }
          }
        }
      }
      if (itsShowFF) {
        os << "Fully flagged baselines: ";
        for (uint i=0; i<fullyFlagged.size(); ++i) {
          if (i>0) os << "; ";
          os << fullyFlagged[i].first << '&' << fullyFlagged[i].second;
        }
        os << endl;
      }
      if (! itsSaveName.empty()) {
        saveStation (npoints, nusedAnt, countAnt);
      }
    }

    void FlagCounter::showChannel (ostream& os, int64 ntimes) const
    {
      int64 npoints = ntimes * itsBLCounts.size();
      int64 nflagged = 0;
      os << endl << "Percentage of visibilities flagged per channel:" << endl;
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
          nflagged += itsChanCounts[ch];
          os << std::setw(4) << int((100. * itsChanCounts[ch]) / npoints + 0.5)
             << '%';
          ch++;
        }
        os << endl;
      }
      int64 totalnpoints = npoints * itsChanCounts.size();
      os << "Total flagged: ";
      showPerc3 (os, nflagged, totalnpoints);
      os << "   (" << nflagged << " out of " << totalnpoints
         << " visibilities)" << endl;
      if (itsWarnPerc > 0) {
        for (uint i=0; i<itsChanCounts.size(); ++i) {
          double perc = (100. * itsChanCounts[i]) / npoints;
          if (perc >= itsWarnPerc) {
            os << "** NOTE: ";
            showPerc1 (os, perc, 100);
            os << " of data are flagged for channel " << i << endl;
          }
        }
      }
      if (! itsSaveName.empty()) {
        saveChannel (npoints, itsChanCounts);
      }
    }

    void FlagCounter::showCorrelation (ostream& os, int64 ntimes) const
    {
      int64 ntotal = ntimes * itsBLCounts.size() * itsChanCounts.size();
      os << endl
         << "Percentage of flagged visibilities detected per correlation:"
         << endl;
      os << "  " << itsCorrCounts << " out of " << ntotal
         << " visibilities   [";
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
      int perc = (total==0  ?  0 : int(1000. * value / total + 0.5));
      os << std::setw(3) << perc/10 << '.' << perc%10 << '%';
    }

    void FlagCounter::showPerc3 (ostream& os, double value, double total)
    {
      int perc = (total==0  ?  0 : int(100000. * value / total + 0.5));
      os << std::setw(5) << perc/1000 << '.';
      // It looks as if std::setfill keeps the fill character, so use
      // ios.fill to be able to reset it.
      char prev = os.fill ('0');
      os << std::setw(3) << perc%1000 << '%';
      os.fill (prev);
    }

    void FlagCounter::saveStation (int64 npoints, const Vector<int64>& nused,
                                   const Vector<int64>& count) const
    {
      // Create the table.
      TableDesc td;
      td.addColumn (ScalarColumnDesc<Int>   ("Station"));
      td.addColumn (ScalarColumnDesc<String>("Name"));
      td.addColumn (ScalarColumnDesc<float> ("Percentage"));
      SetupNewTable newtab(itsSaveName+"stat", td, Table::New);
      Table tab(newtab);
      ScalarColumn<Int>    statCol(tab, "Station");
      ScalarColumn<String> nameCol(tab, "Name");
      ScalarColumn<float>  percCol(tab, "Percentage");
      const Vector<String>& antNames = itsInfo->antennaNames();
      // Write if an antenna is used.
      for (uint i=0; i<nused.size(); ++i) {
        if (nused[i] > 0) {
          int rownr = tab.nrow();
          tab.addRow();
          statCol.put (rownr, i);
          nameCol.put (rownr, antNames[i]);
          percCol.put (rownr, (100. * count[i]) / (nused[i] * npoints));
        }
      }
    }

    void FlagCounter::saveChannel (int64 npoints,
                                   const Vector<int64>& count) const
    {
      // Create the table.
      TableDesc td;
      td.addColumn (ScalarColumnDesc<double>("Frequency"));
      td.addColumn (ScalarColumnDesc<float> ("Percentage"));
      SetupNewTable newtab(itsSaveName+"freq", td, Table::New);
      Table tab(newtab);
      ScalarColumn<double> freqCol(tab, "Frequency");
      ScalarColumn<float>  percCol(tab, "Percentage");
      // Get the channel frequencies.
      const Vector<double>& chanFreqs = itsInfo->chanFreqs();
      for (uint i=0; i<count.size(); ++i) {
        int rownr = tab.nrow();
        tab.addRow();
        freqCol.put (rownr, chanFreqs[i]);
        percCol.put (rownr, (100. * count[i]) / npoints);
      }
    }

  } //# end namespace
}
