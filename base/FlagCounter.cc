// FlagCounter.cc: Class to keep counts of nr of flagged visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "FlagCounter.h"

#include "../steps/InputStep.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <cassert>
#include <vector>
#include <map>
#include <iomanip>
#include <ostream>

using namespace casacore;

using dp3::common::operator<<;

namespace dp3 {
namespace base {

FlagCounter::FlagCounter() : itsInfo(0), itsShowFF(false) {}

FlagCounter::FlagCounter(const string& msName,
                         const common::ParameterSet& parset,
                         const string& prefix) {
  itsWarnPerc = parset.getDouble(prefix + "warnperc", 0);
  itsShowFF = parset.getBool(prefix + "showfullyflagged", false);
  bool save = parset.getBool(prefix + "save", false);
  if (save) {
    // Percentages have to be saved, so form the table name to use.
    string path = parset.getString(prefix + "path", "");
    // Use the step name (without dot) as a name suffix.
    string suffix = prefix;
    string::size_type pos = suffix.find('.');
    if (pos != string::npos) {
      suffix = suffix.substr(0, pos);
    }
    // Use the MS name as the name.
    // If no path is given, use the path of the MS (use . if no path).
    string name = msName;
    pos = name.rfind('/');
    if (path.empty()) {
      if (pos == string::npos) {
        path = '.';
      } else {
        path = name.substr(0, pos);
      }
    }
    name = name.substr(pos + 1);
    pos = name.find('.');
    if (pos != string::npos) {
      name = name.substr(0, pos);
    }
    itsSaveName = path + '/' + name + '_' + suffix + ".flag";
  }
}

void FlagCounter::init(const DPInfo& info) {
  itsInfo = &info;
  itsBLCounts.resize(info.nbaselines());
  itsChanCounts.resize(info.nchan());
  itsCorrCounts.resize(info.ncorr());
  std::fill(itsBLCounts.begin(), itsBLCounts.end(), 0);
  std::fill(itsChanCounts.begin(), itsChanCounts.end(), 0);
  std::fill(itsCorrCounts.begin(), itsCorrCounts.end(), 0);
}

void FlagCounter::add(const FlagCounter& that) {
  // Add that to this after checking for equal sizes.
  assert(itsBLCounts.size() == that.itsBLCounts.size());
  assert(itsChanCounts.size() == that.itsChanCounts.size());
  assert(itsCorrCounts.size() == that.itsCorrCounts.size());
  std::transform(itsBLCounts.begin(), itsBLCounts.end(),
                 that.itsBLCounts.begin(), itsBLCounts.begin(),
                 std::plus<int64_t>());
  std::transform(itsChanCounts.begin(), itsChanCounts.end(),
                 that.itsChanCounts.begin(), itsChanCounts.begin(),
                 std::plus<int64_t>());
  std::transform(itsCorrCounts.begin(), itsCorrCounts.end(),
                 that.itsCorrCounts.begin(), itsCorrCounts.begin(),
                 std::plus<int64_t>());
}

void FlagCounter::showStation(std::ostream& os, int64_t ntimes) const {
  const Vector<Int>& ant1 = itsInfo->getAnt1();
  const Vector<Int>& ant2 = itsInfo->getAnt2();
  const Vector<String>& antNames = itsInfo->antennaNames();

  const int64_t nPoints = ntimes * itsChanCounts.size();
  const size_t nrAnt = antNames.size();

  // Collect counts per antenna.
  std::vector<size_t> nUsedAnt(nrAnt, 0);
  std::vector<size_t> countAnt(nrAnt, 0);

  // Collect ratio of flagged visibilities per antenna
  std::vector<double> flagRatiosPerAntenna(nrAnt, 0);

  for (size_t i = 0; i < itsBLCounts.size(); ++i) {
    countAnt[ant1[i]] += itsBLCounts[i];
    nUsedAnt[ant1[i]]++;
    if (ant1[i] != ant2[i]) {
      countAnt[ant2[i]] += itsBLCounts[i];
      nUsedAnt[ant2[i]]++;
    }
  }

  // Determine nr of antennae used.
  int nrused = 0;
  for (size_t i = 0; i < nrAnt; ++i) {
    if (nUsedAnt[i] > 0) {
      nrused++;
    }
  }

  // Calculate the flag ratios per antenna.
  for (size_t ant = 0; ant != nrAnt; ++ant) {
    if (nUsedAnt[ant] > 0) {
      flagRatiosPerAntenna[ant] =
          double(countAnt[ant]) / double(nUsedAnt[ant] * nPoints);
    }
  }

  // Print antenna names and ratios
  os << "{\"flagged_fraction_dict\": \"{";
  for (size_t i = 0; i < antNames.size(); i++) {
    if (i > 0) {
      os << ", ";
    }
    os << "'" << antNames[i] << "': ";
    os << flagRatiosPerAntenna[i];
  }
  os << "}\"}";
}

void FlagCounter::showBaseline(std::ostream& os, int64_t ntimes) const {
  const Vector<Int>& ant1 = itsInfo->getAnt1();
  const Vector<Int>& ant2 = itsInfo->getAnt2();
  const Vector<String>& antNames = itsInfo->antennaNames();
  // Keep track of fully flagged baselines.
  std::vector<std::pair<int, int>> fullyFlagged;
  int64_t nPoints = ntimes * itsChanCounts.size();
  os << std::endl
     << "Percentage of visibilities flagged per baseline"
        " (antenna pair):";
  unsigned int nrAnt = 1 + std::max(max(ant1), max(ant2));
  // Collect counts per baseline and antenna.
  Vector<int64_t> nUsedAnt(nrAnt, 0);
  Vector<int64_t> countAnt(nrAnt, 0);
  Matrix<int64_t> nusedBL(nrAnt, nrAnt, 0);
  Matrix<int64_t> countBL(nrAnt, nrAnt, 0);
  for (unsigned int i = 0; i < itsBLCounts.size(); ++i) {
    countBL(ant1[i], ant2[i]) += itsBLCounts[i];
    nusedBL(ant1[i], ant2[i])++;
    countAnt[ant1[i]] += itsBLCounts[i];
    nUsedAnt[ant1[i]]++;
    if (ant1[i] != ant2[i]) {
      countBL(ant2[i], ant1[i]) += itsBLCounts[i];
      nusedBL(ant2[i], ant1[i])++;
      countAnt[ant2[i]] += itsBLCounts[i];
      nUsedAnt[ant2[i]]++;
    }
  }
  // Determine nr of antennae used.
  int nrused = 0;
  for (unsigned int i = 0; i < nrAnt; ++i) {
    if (nUsedAnt[i] > 0) {
      nrused++;
    }
  }
  // Print 15 antennae per line.
  const int nantpl = 15;
  int nrl = (nrused + nantpl - 1) / nantpl;
  int ant = 0;
  // Loop over nr of lines needed for the antennae.
  for (int i = 0; i < nrl; ++i) {
    int oldant = ant;
    // Determine nrAnt per line
    int nra = std::min(nantpl, nrused - i * nantpl);
    // Print the header for the antennae being used.
    // It also updates ant for the next iteration.
    os << std::endl << " ant";
    for (int j = 0; j < nra;) {
      if (nUsedAnt[ant] > 0) {
        os << std::setw(5) << ant;
        j++;
      }
      ant++;
    }
    os << std::endl;
    // Print the percentages per antenna pair.
    for (unsigned int k = 0; k < nrAnt; ++k) {
      if (nUsedAnt[k] > 0) {
        os << std::setw(4) << k << " ";
        int ia = oldant;
        for (int j = 0; j < nra;) {
          if (nUsedAnt[ia] > 0) {
            if (nusedBL(k, ia) > 0) {
              os << std::setw(4)
                 << int((100. * countBL(k, ia)) / (nusedBL(k, ia) * nPoints) +
                        0.5)
                 << '%';
              // Determine if baseline is fully flagged.
              // Do it only for ANT1<=ANT2
              if (int(k) <= ia) {
                if (countBL(k, ia) == nusedBL(k, ia) * nPoints) {
                  fullyFlagged.push_back(std::pair<int, int>(k, ia));
                }
              }
            } else {
              os << "     ";
            }
            j++;
          }
          ia++;
        }
        os << std::endl;
      }
    }
    // Print the percentages per antenna.
    os << "TOTAL";
    int ia = oldant;
    for (int j = 0; j < nra;) {
      if (nUsedAnt[ia] > 0) {
        double perc = 100. * countAnt[ia] / (nUsedAnt[ia] * nPoints);
        os << std::setw(4) << int(perc + 0.5) << '%';
        j++;
      }
      ia++;
    }
    os << std::endl;
  }
  if (itsWarnPerc > 0) {
    for (unsigned int i = 0; i < nrAnt; ++i) {
      if (nUsedAnt[i] > 0) {
        double perc = (100. * countAnt[i]) / (nUsedAnt[i] * nPoints);
        if (perc >= itsWarnPerc) {
          os << "** NOTE: ";
          showPerc1(os, perc, 100);
          os << " of data are flagged for station " << i << " (" << antNames[i]
             << ')' << std::endl;
        }
      }
    }
  }
  if (itsShowFF) {
    os << "Fully flagged baselines: ";
    for (unsigned int i = 0; i < fullyFlagged.size(); ++i) {
      if (i > 0) os << "; ";
      os << fullyFlagged[i].first << '&' << fullyFlagged[i].second;
    }
    os << std::endl;
  }
  if (!itsSaveName.empty()) {
    saveStation(nPoints, nUsedAnt, countAnt);
  }
}

void FlagCounter::showChannel(std::ostream& os, int64_t ntimes) const {
  int64_t nPoints = ntimes * itsBLCounts.size();
  int64_t nflagged = 0;
  os << std::endl
     << "Percentage of visibilities flagged per channel:" << std::endl;
  if (nPoints == 0) {
    return;
  }
  // Print 10 channels per line.
  const int nchpl = 10;
  os << " channels    ";
  for (int i = 0; i < std::min(nchpl, int(itsChanCounts.size())); ++i) {
    os << std::setw(5) << i;
  }
  os << std::endl;
  int nrl = (itsChanCounts.size() + nchpl - 1) / nchpl;
  int ch = 0;
  for (int i = 0; i < nrl; ++i) {
    int nrc = std::min(nchpl, int(itsChanCounts.size() - i * nchpl));
    os << std::setw(4) << ch << '-' << std::setw(4) << ch + nrc - 1 << ":    ";
    for (int j = 0; j < nrc; ++j) {
      nflagged += itsChanCounts[ch];
      os << std::setw(4) << int((100. * itsChanCounts[ch]) / nPoints + 0.5)
         << '%';
      ch++;
    }
    os << std::endl;
  }
  int64_t totalnPoints = nPoints * itsChanCounts.size();
  // Prevent division by zero
  if (totalnPoints == 0) {
    totalnPoints = 1;
  }
  os << "Total flagged: ";
  showPerc3(os, nflagged, totalnPoints);
  os << "   (" << nflagged << " out of " << totalnPoints << " visibilities)"
     << std::endl;
  if (itsWarnPerc > 0) {
    for (unsigned int i = 0; i < itsChanCounts.size(); ++i) {
      double perc = (100. * itsChanCounts[i]) / nPoints;
      if (perc >= itsWarnPerc) {
        os << "** NOTE: ";
        showPerc1(os, perc, 100);
        os << " of data are flagged for channel " << i << std::endl;
      }
    }
  }
  if (!itsSaveName.empty()) {
    saveChannel(nPoints, itsChanCounts);
  }
}

void FlagCounter::showCorrelation(std::ostream& os, int64_t ntimes) const {
  int64_t ntotal = ntimes * itsBLCounts.size() * itsChanCounts.size();
  // Prevent division by zero
  if (ntotal == 0) {
    ntotal = 1;
  }
  os << '\n'
     << "Percentage of flagged visibilities detected per correlation:" << '\n'
     << "  " << itsCorrCounts << " out of " << ntotal << " visibilities   [";
  for (unsigned int i = 0; i < itsCorrCounts.size(); ++i) {
    if (i > 0) {
      os << ", ";
    }
    os << int(100. * itsCorrCounts[i] / ntotal + 0.5) << '%';
  }
  os << ']' << std::endl;
}

void FlagCounter::showPerc1(ostream& os, double value, double total) {
  int perc = (total == 0 ? 0 : int(1000. * value / total + 0.5));
  os << std::setw(3) << perc / 10 << '.' << perc % 10 << '%';
}

void FlagCounter::showPerc3(ostream& os, double value, double total) {
  int perc = (total == 0 ? 0 : int(100000. * value / total + 0.5));
  os << std::setw(5) << perc / 1000 << '.';
  // It looks as if std::setfill keeps the fill character, so use
  // ios.fill to be able to reset it.
  char prev = os.fill('0');
  os << std::setw(3) << perc % 1000 << '%';
  os.fill(prev);
}

void FlagCounter::saveStation(int64_t nPoints, const Vector<int64_t>& nused,
                              const Vector<int64_t>& count) const {
  // Create the table.
  TableDesc td;
  td.addColumn(ScalarColumnDesc<Int>("Station"));
  td.addColumn(ScalarColumnDesc<String>("Name"));
  td.addColumn(ScalarColumnDesc<float>("Percentage"));
  SetupNewTable newtab(itsSaveName + "stat", td, Table::New);
  Table tab(newtab);
  ScalarColumn<Int> statCol(tab, "Station");
  ScalarColumn<String> nameCol(tab, "Name");
  ScalarColumn<float> percCol(tab, "Percentage");
  const Vector<String>& antNames = itsInfo->antennaNames();
  // Write if an antenna is used.
  for (unsigned int i = 0; i < nused.size(); ++i) {
    if (nused[i] > 0) {
      int rownr = tab.nrow();
      tab.addRow();
      statCol.put(rownr, i);
      nameCol.put(rownr, antNames[i]);
      percCol.put(rownr, (100. * count[i]) / (nused[i] * nPoints));
    }
  }
}

void FlagCounter::saveChannel(int64_t nPoints,
                              const Vector<int64_t>& count) const {
  // Create the table.
  TableDesc td;
  td.addColumn(ScalarColumnDesc<double>("Frequency"));
  td.addColumn(ScalarColumnDesc<float>("Percentage"));
  SetupNewTable newtab(itsSaveName + "freq", td, Table::New);
  Table tab(newtab);
  ScalarColumn<double> freqCol(tab, "Frequency");
  ScalarColumn<float> percCol(tab, "Percentage");
  // Get the channel frequencies.
  const Vector<double>& chanFreqs = itsInfo->chanFreqs();
  for (unsigned int i = 0; i < count.size(); ++i) {
    int rownr = tab.nrow();
    tab.addRow();
    freqCol.put(rownr, chanFreqs[i]);
    percCol.put(rownr, (100. * count[i]) / nPoints);
  }
}

}  // namespace base
}  // namespace dp3
