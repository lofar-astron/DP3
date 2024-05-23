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
  /// Parse an MSSelection baseline string and create a Matrix telling
  /// which baselines are selected.
  /// @param names Name for each station/antenna.
  /// @param positions Position for each station/antenna.
  /// @param antenna1 For each baseline, the index of the first antenna.
  /// @param antenna2 For each baseline, the index of the second antenna.
  /// @param baseline_selection Selection string in MSSelection format.
  /// @param os Possible messages from the parser are written to this stream.
  /// @return An n_stations x n_stations Matrix, which holds true for selected
  //          baselines and false for the other baselines.
  static casacore::Matrix<bool> convert(
      const std::vector<std::string>& names,
      const std::vector<casacore::MPosition>& positions,
      const std::vector<int>& antenna1, const std::vector<int>& antenna2,
      const std::string& baseline_selection, std::ostream&);

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
  ~BaselineSelectErrorHandler() override;
  void reportError(const char* token, const casacore::String message) override;

 private:
  std::ostream& itsStream;
};

/// @}

}  // namespace common
}  // namespace dp3

#endif
