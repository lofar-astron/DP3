// BaselineSelection.h: Class to handle the baseline selection
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to handle the baseline selection
/// @author Ger van Diepen

#ifndef DPPP_BASELINESELECTION_H
#define DPPP_BASELINESELECTION_H

#include "DPInfo.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>

namespace dp3 {
namespace common {
class ParameterSet;
class ParameterValue;
}  // namespace common

namespace base {

/// \brief Class containing a few static functions to parse a baseline selection
/// string.
class BaselineSelection {
 public:
  /// Default constructor has no selection.
  BaselineSelection();

  /// Construct from the parset using the given prefix.
  /// The keys used are:
  /// <ul>
  ///  <li> baseline: for a baseline selection (e.g. CS*&)
  ///  <li> corrtype: for correlation selection (auto, cross, or empty)
  ///  <li> blrange:  ranges of baseline lengths (in m)
  ///  <li> minbl:    minimum baseline length (in m); only if minmax=true
  ///  <li> maxbl:    maximum baseline length (in m); only if minmax=true
  /// </ul>
  BaselineSelection(const common::ParameterSet&, const string& prefix,
                    bool minmax = false,
                    const string& defaultCorrType = string(),
                    const string& defaultBaseline = string());

  /// Is there any selection?
  bool hasSelection() const;

  /// Show the parameters.
  /// Optional extra blanks can be put before the value.
  void show(std::ostream& os, const std::string& blanks = string()) const;

  /// Form the selection matrix telling for each baseline if it is selected.
  /// If no selection is made, all values in the matrix are true.
  casacore::Matrix<bool> apply(const DPInfo& info) const;

  /// Form the selection vector telling if a baseline in the DPInfo object
  /// is selected.
  casacore::Vector<bool> applyVec(const DPInfo& info) const;

 private:
  /// Convert the baseline selection string.
  void handleBL(casacore::Matrix<bool>& selectBL, const DPInfo& info) const;

  /// Handle a vector of baseline specifications.
  casacore::Matrix<bool> handleBLVector(
      const common::ParameterValue& pvBL,
      const casacore::Vector<casacore::String>&) const;

  /// Handle the correlation type selection.
  void handleCorrType(casacore::Matrix<bool>& selectBL) const;

  /// Handle the baseline length selection.
  void handleLength(casacore::Matrix<bool>& selectBL, const DPInfo& info) const;

  string itsStrBL;
  string itsCorrType;
  std::vector<double> itsRangeBL;
};

}  // namespace base
}  // namespace dp3

#endif
