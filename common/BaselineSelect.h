// BaselineSelect.h: Convert MSSelection baseline string to a Matrix
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// Convert MSSelection baseline string to a Matrix
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef MS_BASELINESELECT_H
#define MS_BASELINESELECT_H

#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>
#include <casacore/casa/Arrays/Array.h>

#include <ostream>

namespace casacore {
class Table;
class TableExprNode;
class MPosition;
}  // namespace casacore

namespace dp3 {
namespace common {

/// @ingroup MS
/// @brief Convert MSSelection baseline string to a Matrix
/// @{

/// Class with a static function to convert a casacore MSSelection baseline
/// string to a Matrix<Bool> telling which baselines are selected.

class BaselineSelect {
 public:
  /// Parse the MSSelection baseline string and create a Matrix telling
  /// which baselines are selected.
  /// Possible messages from the parser are written to the ostream.
  static casacore::Matrix<bool> convert(const string& msName,
                                        const string& baselineSelection,
                                        std::ostream&);

  /// Parse the MSSelection baseline string and create a Matrix telling
  /// which baselines are selected.
  /// The input is a vector of station names and positions.
  /// Possible messages from the parser are written to the ostream.
  static casacore::Matrix<bool> convert(
      const casacore::Vector<casacore::String>& names,
      const std::vector<casacore::MPosition>& pos,
      const casacore::Vector<casacore::Int>& ant1,
      const casacore::Vector<casacore::Int>& ant2,
      const string& baselineSelection, std::ostream&);

 private:
  static casacore::Matrix<bool> convert(casacore::Table& anttab,
                                        casacore::TableExprNode& a1,
                                        casacore::TableExprNode& a2,
                                        const string& baselineSelection,
                                        std::ostream& os);
};

/// @brief This class handles an error from the Casacore's MSAntennaParse.
/// It adds the message to the message list of the parent BaselineSelect.
class BaselineSelectErrorHandler : public casacore::MSSelectionErrorHandler {
 public:
  BaselineSelectErrorHandler(std::ostream& os) : itsStream(os) {}
  virtual ~BaselineSelectErrorHandler();
  virtual void reportError(const char* token, const casacore::String message);

 private:
  std::ostream& itsStream;
};

/// @}

}  // namespace common
}  // namespace dp3

#endif
