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

static std::string MakeSaveFilename(std::string path,
                                    const std::string& ms_name,
                                    std::string suffix) {
  // Use the step name (without dot) as a name suffix.
  std::string::size_type pos = suffix.find('.');
  if (pos != std::string::npos) {
    suffix.resize(pos);
  }
  // If no path is given, use the path of the name (use . if no path).
  pos = ms_name.rfind('/');
  if (path.empty()) {
    if (pos == std::string::npos) {
      path = '.';
    } else {
      path = ms_name.substr(0, pos);
    }
  }
  std::string name = ms_name.substr(pos + 1);
  pos = name.find('.');
  if (pos != std::string::npos) {
    name = name.substr(0, pos);
  }
  return path + '/' + name + '_' + suffix + ".flag";
}

FlagCounter::FlagCounter(const common::ParameterSet& parset,
                         const std::string& prefix)
    : warning_percentage_(parset.getDouble(prefix + "warnperc", 0)),
      show_fully_flagged_(parset.getBool(prefix + "showfullyflagged", false)),
      save_(parset.getBool(prefix + "save", false)),
      path_(parset.getString(prefix + "path", "")),
      name_(prefix) {}

void FlagCounter::init(const DPInfo& info) {
  info_ = &info;

  if (save_)
    // Percentages have to be saved, so form the table name to use.
    save_filename_ = MakeSaveFilename(path_, info.msName(), name_);
  base_line_counts_.resize(info.nbaselines());
  channel_counts_.resize(info.nchan());
  correlation_counts_.resize(info.ncorr());
  std::fill(base_line_counts_.begin(), base_line_counts_.end(), 0);
  std::fill(channel_counts_.begin(), channel_counts_.end(), 0);
  std::fill(correlation_counts_.begin(), correlation_counts_.end(), 0);
}

void FlagCounter::add(const FlagCounter& that) {
  // Add that to this after checking for equal sizes.
  assert(base_line_counts_.size() == that.base_line_counts_.size());
  assert(channel_counts_.size() == that.channel_counts_.size());
  assert(correlation_counts_.size() == that.correlation_counts_.size());
  std::transform(base_line_counts_.begin(), base_line_counts_.end(),
                 that.base_line_counts_.begin(), base_line_counts_.begin(),
                 std::plus<int64_t>());
  std::transform(channel_counts_.begin(), channel_counts_.end(),
                 that.channel_counts_.begin(), channel_counts_.begin(),
                 std::plus<int64_t>());
  std::transform(correlation_counts_.begin(), correlation_counts_.end(),
                 that.correlation_counts_.begin(), correlation_counts_.begin(),
                 std::plus<int64_t>());
}

void FlagCounter::showStation(std::ostream& os, int64_t ntimes) const {
  const std::vector<int>& ant1 = info_->getAnt1();
  const std::vector<int>& ant2 = info_->getAnt2();
  const std::vector<std::string>& antNames = info_->antennaNames();

  const int64_t nPoints = ntimes * channel_counts_.size();
  const size_t nrAnt = antNames.size();

  // Collect counts per antenna.
  std::vector<size_t> nUsedAnt(nrAnt, 0);
  std::vector<size_t> countAnt(nrAnt, 0);

  // Collect ratio of flagged visibilities per antenna
  std::vector<double> flagRatiosPerAntenna(nrAnt, 0);

  for (size_t i = 0; i < base_line_counts_.size(); ++i) {
    countAnt[ant1[i]] += base_line_counts_[i];
    nUsedAnt[ant1[i]]++;
    if (ant1[i] != ant2[i]) {
      countAnt[ant2[i]] += base_line_counts_[i];
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
  const std::vector<int>& ant1 = info_->getAnt1();
  const std::vector<int>& ant2 = info_->getAnt2();
  const std::vector<std::string>& antNames = info_->antennaNames();
  // Keep track of fully flagged baselines.
  std::vector<std::pair<int, int>> fullyFlagged;
  int64_t nPoints = ntimes * channel_counts_.size();
  os << std::endl
     << "Percentage of visibilities flagged per baseline"
        " (antenna pair):";
  std::vector<int>::const_iterator ant1_max =
      std::max_element(ant1.begin(), ant1.end());
  std::vector<int>::const_iterator ant2_max =
      std::max_element(ant2.begin(), ant2.end());
  const unsigned int nrAnt =
      (ant1_max == ant1.end()) ? 0 : 1 + std::max(*ant1_max, *ant2_max);
  // Collect counts per baseline and antenna.
  Vector<int64_t> nUsedAnt(nrAnt, 0);
  Vector<int64_t> countAnt(nrAnt, 0);
  Matrix<int64_t> nusedBL(nrAnt, nrAnt, 0);
  Matrix<int64_t> countBL(nrAnt, nrAnt, 0);
  for (unsigned int i = 0; i < base_line_counts_.size(); ++i) {
    countBL(ant1[i], ant2[i]) += base_line_counts_[i];
    nusedBL(ant1[i], ant2[i])++;
    countAnt[ant1[i]] += base_line_counts_[i];
    nUsedAnt[ant1[i]]++;
    if (ant1[i] != ant2[i]) {
      countBL(ant2[i], ant1[i]) += base_line_counts_[i];
      nusedBL(ant2[i], ant1[i])++;
      countAnt[ant2[i]] += base_line_counts_[i];
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
                 << int((100.0 * countBL(k, ia)) / (nusedBL(k, ia) * nPoints) +
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
        double perc = 100.0 * countAnt[ia] / (nUsedAnt[ia] * nPoints);
        os << std::setw(4) << int(perc + 0.5) << '%';
        j++;
      }
      ia++;
    }
    os << std::endl;
  }
  if (warning_percentage_ > 0) {
    for (unsigned int i = 0; i < nrAnt; ++i) {
      if (nUsedAnt[i] > 0) {
        double perc = (100.0 * countAnt[i]) / (nUsedAnt[i] * nPoints);
        if (perc >= warning_percentage_) {
          os << "** NOTE: ";
          showPerc1(os, perc, 100);
          os << " of data are flagged for station " << i << " (" << antNames[i]
             << ')' << std::endl;
        }
      }
    }
  }
  if (show_fully_flagged_) {
    os << "Fully flagged baselines: ";
    for (unsigned int i = 0; i < fullyFlagged.size(); ++i) {
      if (i > 0) os << "; ";
      os << fullyFlagged[i].first << '&' << fullyFlagged[i].second;
    }
    os << std::endl;
  }
  if (!save_filename_.empty()) {
    saveStation(nPoints, nUsedAnt, countAnt);
  }
}

void FlagCounter::showChannel(std::ostream& os, int64_t ntimes) const {
  int64_t nPoints = ntimes * base_line_counts_.size();
  int64_t nflagged = 0;
  os << std::endl
     << "Percentage of visibilities flagged per channel:" << std::endl;
  if (nPoints == 0) {
    return;
  }
  // Print 10 channels per line.
  const int nchpl = 10;
  os << " channels    ";
  for (int i = 0; i < std::min(nchpl, int(channel_counts_.size())); ++i) {
    os << std::setw(5) << i;
  }
  os << std::endl;
  int nrl = (channel_counts_.size() + nchpl - 1) / nchpl;
  int ch = 0;
  for (int i = 0; i < nrl; ++i) {
    int nrc = std::min(nchpl, int(channel_counts_.size() - i * nchpl));
    os << std::setw(4) << ch << '-' << std::setw(4) << ch + nrc - 1 << ":    ";
    for (int j = 0; j < nrc; ++j) {
      nflagged += channel_counts_[ch];
      os << std::setw(4) << int((100.0 * channel_counts_[ch]) / nPoints + 0.5)
         << '%';
      ch++;
    }
    os << std::endl;
  }
  int64_t totalnPoints = nPoints * channel_counts_.size();
  // Prevent division by zero
  if (totalnPoints == 0) {
    totalnPoints = 1;
  }
  os << "Total flagged: ";
  showPerc3(os, nflagged, totalnPoints);
  os << "   (" << nflagged << " out of " << totalnPoints << " visibilities)"
     << std::endl;
  if (warning_percentage_ > 0) {
    for (unsigned int i = 0; i < channel_counts_.size(); ++i) {
      double perc = (100.0 * channel_counts_[i]) / nPoints;
      if (perc >= warning_percentage_) {
        os << "** NOTE: ";
        showPerc1(os, perc, 100);
        os << " of data are flagged for channel " << i << std::endl;
      }
    }
  }
  if (!save_filename_.empty()) {
    saveChannel(nPoints, channel_counts_);
  }
}

void FlagCounter::showCorrelation(std::ostream& os, int64_t ntimes) const {
  int64_t ntotal = ntimes * base_line_counts_.size() * channel_counts_.size();
  // Prevent division by zero
  if (ntotal == 0) {
    ntotal = 1;
  }
  os << '\n'
     << "Percentage of flagged visibilities detected per correlation:" << '\n'
     << "  " << correlation_counts_ << " out of " << ntotal
     << " visibilities   [";
  for (unsigned int i = 0; i < correlation_counts_.size(); ++i) {
    if (i > 0) {
      os << ", ";
    }
    os << int(100.0 * correlation_counts_[i] / ntotal + 0.5) << '%';
  }
  os << ']' << std::endl;
}

static void WriteDuration(ostream& os, double duration) {
  if (duration < 10) {
    duration *= 1000;
    os << std::setw(5) << int(duration) << " ms";
  } else {
    os << std::setw(5) << int(duration) << "  s";
  }
}

void FlagCounter::showPerc1(ostream& os, double value, double total) {
  int perc = (total == 0 ? 0 : int(1000.0 * value / total + 0.5));
  os << std::setw(3) << perc / 10 << '.' << perc % 10 << "% (";
  WriteDuration(os, value);
  os << ')';
}

void FlagCounter::showPerc3(ostream& os, double value, double total) {
  int perc = (total == 0 ? 0 : int(100000.0 * value / total + 0.5));
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
  SetupNewTable newtab(save_filename_ + "stat", td, Table::New);
  Table tab(newtab);
  ScalarColumn<Int> statCol(tab, "Station");
  ScalarColumn<String> nameCol(tab, "Name");
  ScalarColumn<float> percCol(tab, "Percentage");
  const std::vector<std::string>& antNames = info_->antennaNames();
  // Write if an antenna is used.
  for (unsigned int i = 0; i < nused.size(); ++i) {
    if (nused[i] > 0) {
      int rownr = tab.nrow();
      tab.addRow();
      statCol.put(rownr, i);
      nameCol.put(rownr, antNames[i]);
      percCol.put(rownr, (100.0 * count[i]) / (nused[i] * nPoints));
    }
  }
}

void FlagCounter::saveChannel(int64_t nPoints,
                              const std::vector<int64_t>& count) const {
  // Create the table.
  TableDesc td;
  td.addColumn(ScalarColumnDesc<double>("Frequency"));
  td.addColumn(ScalarColumnDesc<float>("Percentage"));
  SetupNewTable newtab(save_filename_ + "freq", td, Table::New);
  Table tab(newtab);
  ScalarColumn<double> freqCol(tab, "Frequency");
  ScalarColumn<float> percCol(tab, "Percentage");
  // Get the channel frequencies.
  const std::vector<double>& chanFreqs = info_->chanFreqs();
  for (unsigned int i = 0; i < count.size(); ++i) {
    int rownr = tab.nrow();
    tab.addRow();
    freqCol.put(rownr, chanFreqs[i]);
    percCol.put(rownr, (100.0 * count[i]) / nPoints);
  }
}

}  // namespace base
}  // namespace dp3
