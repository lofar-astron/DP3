//# BaselineSelect.h: Convert MSSelection baseline string to a Matrix
//#
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
//#  $Id: BaselineSelect.h 34753 2016-06-20 10:43:42Z schaap $

#ifndef MS_BASELINESELECT_H
#define MS_BASELINESELECT_H

// @file
// Convert MSSelection baseline string to a Matrix
// @author Ger van Diepen (diepen AT astron nl)

//# Includes
#include <casacore/ms/MSSel/MSSelectionErrorHandler.h>

#include <ostream>

//# Forward Declarations
namespace casacore
{
  template<class T> class Matrix;
  template<class T> class Vector;
  class Table;
  class TableExprNode;
  class MPosition;
}

namespace DP3
{

// @ingroup MS
// @brief Convert MSSelection baseline string to a Matrix
// @{

// Class with a static function to convert a casacore MSSelection baseline
// string to a Matrix<Bool> telling which baselines are selected.

class BaselineSelect
{
public:
  // Parse the MSSelection baseline string and create a Matrix telling
  // which baselines are selected.
  // Possible messages from the parser are written to the ostream.
  static casacore::Matrix<bool> convert (const string& msName,
                                     const string& baselineSelection,
                                     std::ostream&);

  // Parse the MSSelection baseline string and create a Matrix telling
  // which baselines are selected.
  // The input is a vector of station names and positions.
  // Possible messages from the parser are written to the ostream.
  static casacore::Matrix<bool> convert (const casacore::Vector<casacore::String>& names,
                                     const std::vector<casacore::MPosition>& pos,
                                     const casacore::Vector<casacore::Int>& ant1,
                                     const casacore::Vector<casacore::Int>& ant2,
                                     const string& baselineSelection,
                                     std::ostream&);

private:
  static casacore::Matrix<bool> convert (casacore::Table& anttab,
                                     casacore::TableExprNode& a1,
                                     casacore::TableExprNode& a2,
                                     const string& baselineSelection,
                                     std::ostream& os);

};



// This class handles an error from the Casacore's MSAntennaParse.
// It adds the message to the message list of the parent BaselineSelect.
class BaselineSelectErrorHandler : public casacore::MSSelectionErrorHandler
{
public:
  BaselineSelectErrorHandler (std::ostream& os)
    : itsStream (os)
  {}
  virtual ~BaselineSelectErrorHandler();
  virtual void reportError (const char *token, const casacore::String message);
private:
  std::ostream& itsStream;
};


// @}

} // end namespace

#endif
