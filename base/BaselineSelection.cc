// BaselineSelection.cc: Class to handle the baseline selection
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "BaselineSelection.h"

#include <boost/algorithm/string.hpp>

#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/Arrays/Matrix.h>

#include "../common/BaselineSelect.h"
#include "../common/ParameterSet.h"
#include "../common/ParameterValue.h"
#include "../common/StreamUtil.h"

#include <vector>

#include <aocommon/logger.h>

using casacore::IPosition;
using casacore::Matrix;
using dp3::common::operator<<;

using aocommon::Logger;

namespace dp3 {
namespace base {

BaselineSelection::BaselineSelection() {}

BaselineSelection::BaselineSelection(const common::ParameterSet& parset,
                                     const string& prefix, bool minmax,
                                     const string& defaultCorrType,
                                     const string& defaultBaseline)
    : itsStrBL(parset.getString(prefix + "baseline", defaultBaseline)),
      itsCorrType(parset.getString(prefix + "corrtype", defaultCorrType)),
      itsRangeBL(
          parset.getDoubleVector(prefix + "blrange", std::vector<double>())) {
  if (minmax) {
    double minbl = parset.getDouble(prefix + "blmin", -1);
    double maxbl = parset.getDouble(prefix + "blmax", -1);
    if (minbl > 0) {
      itsRangeBL.push_back(0.);
      itsRangeBL.push_back(minbl);
    }
    if (maxbl > 0) {
      itsRangeBL.push_back(maxbl);
      itsRangeBL.push_back(1e30);
    }
  }
  if (itsRangeBL.size() % 2 != 0)
    throw std::runtime_error(
        "DP3 error: uneven number of lengths in baseline range");
}

bool BaselineSelection::hasSelection() const {
  return !((itsStrBL.empty() || itsStrBL == "[]") && itsCorrType.empty() &&
           itsRangeBL.empty());
}

void BaselineSelection::show(std::ostream& os,
                             const std::string& blanks) const {
  os << "  Baseline selection:" << '\n';
  os << "    baseline:     " << blanks << itsStrBL << '\n';
  os << "    corrtype:     " << blanks << itsCorrType << '\n';
  os << "    blrange:      " << blanks << itsRangeBL << '\n';
}

Matrix<bool> BaselineSelection::apply(const DPInfo& info) const {
  // Size and initialize the selection matrix.
  int nant = info.antennaNames().size();
  Matrix<bool> selectBL(nant, nant, true);
  // Apply the various parts if given.
  if (!itsStrBL.empty() && itsStrBL != "[]") {
    handleBL(selectBL, info);
  }
  if (!itsCorrType.empty()) {
    handleCorrType(selectBL);
  }
  if (!itsRangeBL.empty()) {
    handleLength(selectBL, info);
  }
  return selectBL;
}

casacore::Vector<bool> BaselineSelection::applyVec(const DPInfo& info) const {
  Matrix<bool> sel = apply(info);
  casacore::Vector<bool> vec;
  vec.resize(info.nbaselines());
  for (unsigned int i = 0; i < info.nbaselines(); ++i) {
    vec[i] = sel(info.getAnt1()[i], info.getAnt2()[i]);
  }
  return vec;
}

void BaselineSelection::handleBL(Matrix<bool>& selectBL,
                                 const DPInfo& info) const {
  // Handle the value(s) in the baseline selection string.
  common::ParameterValue pvBL(itsStrBL);
  // The value can be a vector or an MSSelection string.
  // Alas the ParameterValue vector test cannot be used, because
  // the first character of a MSSelection string can also be [.
  // So if the first is [ and a ] is found before the end and before
  // another [, it must be a MSSelection string.
  bool mssel = true;
  if (itsStrBL[0] == '[') {
    std::string::size_type rb = itsStrBL.find(']');
    if (rb == string::npos)
      throw std::runtime_error("Baseline selection " + itsStrBL +
                               " has no ending ]");
    if (rb == itsStrBL.size() - 1) {
      mssel = false;
    } else {
      std::string::size_type lb = itsStrBL.find('[', 1);
      mssel = (lb == string::npos || lb > rb);
    }
  }
  if (!mssel) {
    // Specified as a vector of antenna name patterns.
    selectBL =
        selectBL && handleBLVector(pvBL, casacore::Vector<casacore::String>(
                                             info.antennaNames()));
  } else {
    // Specified in casacore's MSSelection format.
    string msName = info.msName();
    if (msName.empty()) throw std::runtime_error("Empty measurement set name");
    std::ostringstream os;
    Matrix<bool> sel(common::BaselineSelect::convert(msName, itsStrBL, os));
    // Show possible messages about unknown stations.
    if (!os.str().empty()) {
      Logger::Warn << os.str();
    }
    // The resulting matrix can be smaller because new stations might have
    // been added that are not present in the MS's ANTENNA table.
    if (sel.nrow() == selectBL.nrow()) {
      selectBL = selectBL && sel;
    } else {
      // Only and the subset.
      Matrix<bool> selBL =
          selectBL(IPosition(2, 0), IPosition(2, sel.nrow() - 1));
      selBL = selBL && sel;
    }
  }
}

Matrix<bool> BaselineSelection::handleBLVector(
    const common::ParameterValue& pvBL,
    const casacore::Vector<casacore::String>& antNames) const {
  Matrix<bool> sel(antNames.size(), antNames.size());
  sel = false;
  std::vector<common::ParameterValue> pairs = pvBL.getVector();
  // Each ParameterValue can be a single value (antenna) or a pair of
  // values (a baseline).
  // Note that [ant1,ant2] is somewhat ambiguous; it means two antennae,
  // but one might think it means a baseline [[ant1,ant2]].
  if (pairs.size() == 2 && !(pairs[0].isVector() || pairs[1].isVector())) {
    Logger::Warn << "PreFlagger baseline " << pvBL.get()
                 << " means two antennae, but is somewhat ambigious; "
                    "it's more clear to use [[ant1],[ant2]]\n";
  }
  for (unsigned int i = 0; i < pairs.size(); ++i) {
    std::vector<string> bl = pairs[i].getStringVector();
    if (bl.size() == 1) {
      // Turn the given antenna name pattern into a regex.
      casacore::Regex regex(casacore::Regex::fromPattern(bl[0]));
      int nmatch = 0;
      // Loop through all antenna names and set matrix for matching ones.
      for (unsigned int i2 = 0; i2 < antNames.size(); ++i2) {
        if (casacore::String(antNames[i2]).matches(regex)) {
          nmatch++;
          // Antenna matches, so set all corresponding flags.
          for (unsigned int j = 0; j < antNames.size(); ++j) {
            sel(i2, j) = true;
            sel(j, i2) = true;
          }
        }
      }
      if (nmatch == 0) {
        Logger::Warn << "PreFlagger: no matches for antenna name pattern ["
                     << bl[0] << "]\n";
      }
    } else {
      if (bl.size() != 2)
        throw std::runtime_error(
            "PreFlagger baseline "
            " should contain 1 or 2 antenna name patterns");
      // Turn the given antenna name pattern into a regex.
      casacore::Regex regex1(casacore::Regex::fromPattern(bl[0]));
      casacore::Regex regex2(casacore::Regex::fromPattern(bl[1]));
      int nmatch = 0;
      // Loop through all antenna names and set matrix for matching ones.
      for (unsigned int i2 = 0; i2 < antNames.size(); ++i2) {
        if (casacore::String(antNames[i2]).matches(regex2)) {
          // Antenna2 matches, now try Antenna1.
          for (unsigned int i1 = 0; i1 < antNames.size(); ++i1) {
            if (casacore::String(antNames[i1]).matches(regex1)) {
              nmatch++;
              sel(i1, i2) = true;
              sel(i2, i1) = true;
            }
          }
        }
      }
      if (nmatch == 0) {
        Logger::Warn << "PreFlagger: no matches for baseline name pattern ["
                     << bl[0] << ',' << bl[1] << "]\n";
      }
    }
  }
  return sel;
}

void BaselineSelection::handleCorrType(Matrix<bool>& selectBL) const {
  // Process corrtype if given.
  string corrType = boost::to_lower_copy(itsCorrType);
  if (corrType != "auto" && corrType != "cross")
    throw std::runtime_error(
        "DP3 corrType " + corrType +
        " is invalid; must be auto, cross, or empty string");
  if (corrType == "auto") {
    casacore::Vector<bool> diag = selectBL.diagonal().copy();
    selectBL = false;
    selectBL.diagonal() = diag;
  } else {
    selectBL.diagonal() = false;
  }
}

void BaselineSelection::handleLength(Matrix<bool>& selectBL,
                                     const DPInfo& info) const {
  // Get baseline lengths.
  const std::vector<double>& blength = info.getBaselineLengths();
  const std::vector<int>& ant1 = info.getAnt1();
  const std::vector<int>& ant2 = info.getAnt2();
  for (unsigned int i = 0; i < ant1.size(); ++i) {
    // Clear selection if no range matches.
    bool match = false;
    for (unsigned int j = 0; j < itsRangeBL.size(); j += 2) {
      if (blength[i] >= itsRangeBL[j] && blength[i] <= itsRangeBL[j + 1]) {
        match = true;
        break;
      }
    }
    if (!match) {
      int a1 = ant1[i];
      int a2 = ant2[i];
      selectBL(a1, a2) = false;
      selectBL(a2, a1) = false;
    }
  }
}

}  // namespace base
}  // namespace dp3
