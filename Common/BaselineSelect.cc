//# BaselineSelect.cc: Convert MSSelection baseline string to a Matrix
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
//#  $Id: BaselineSelect.cc 34753 2016-06-20 10:43:42Z schaap $

#include "BaselineSelect.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MSSel/MSAntennaGram.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>

using namespace casacore;

namespace DP3 {

  Matrix<bool> BaselineSelect::convert (const string& msName,
                                        const string& baselineSelection,
                                        std::ostream& os)
  {
    // Find the unique baselines in the MS.
    // Do not use unique sort, because that is slow for a large MS.
    // Simply go through all baselines.
    Table bltab;
    {
      Table tab(msName);
      Vector<Int> a1 = ROScalarColumn<Int>(tab, "ANTENNA1").getColumn();
      Vector<Int> a2 = ROScalarColumn<Int>(tab, "ANTENNA2").getColumn();
      int nant = 1 + std::max(max(a1), max(a2));
      Matrix<bool> bl(nant, nant, false);
      vector<uInt> rows;
      rows.reserve (nant*nant);
      for (unsigned int i=0; i<a1.size(); ++i) {
        if (! bl(a1[i], a2[i])) {
          rows.push_back (i);
          bl(a1[i], a2[i]) = true;
        }
      }
      bltab = tab(Vector<uInt>(rows));
    }
    TableExprNode a1 (bltab.col("ANTENNA1"));
    TableExprNode a2 (bltab.col("ANTENNA2"));
    Table anttab (bltab.keywordSet().asTable("ANTENNA"));
    return convert (anttab, a1, a2, baselineSelection, os);
  }

  casacore::Matrix<bool> BaselineSelect::convert (const Vector<String>& names,
                                              const vector<MPosition>& pos,
                                              const Vector<Int>& ant1,
                                              const Vector<Int>& ant2,
                                              const string& baselineSelection,
                                              std::ostream& os)
  {
    if (names.size() != pos.size())
      throw std::invalid_argument("Name and position arrays are of different size");
    // Create a temporary MSAntenna table in memory for parsing purposes.
    SetupNewTable antNew(String(), MSAntenna::requiredTableDesc(),
                         Table::New);
    Table anttab(antNew, Table::Memory, names.size());
    MSAntenna msant(anttab);
    MSAntennaColumns antcol(msant);
    antcol.name().putColumn (names);
    for (size_t i=0; i<pos.size(); ++i) {
      antcol.positionMeas().put (i, pos[i]);
    }
    // Create a temporary table holding the antenna numbers of the baselines.
    TableDesc td;
    td.addColumn (ScalarColumnDesc<Int>("ANTENNA1"));
    td.addColumn (ScalarColumnDesc<Int>("ANTENNA2"));
    SetupNewTable tabNew(String(), td, Table::New);
    Table tab(tabNew, Table::Memory, ant1.size());
    ScalarColumn<Int> ac1(tab, "ANTENNA1");
    ScalarColumn<Int> ac2(tab, "ANTENNA2");
    ac1.putColumn (ant1);
    ac2.putColumn (ant2);
    // Do the selection using the temporary tables.
    TableExprNode a1 (tab.col("ANTENNA1"));
    TableExprNode a2 (tab.col("ANTENNA2"));
    return convert (anttab, a1, a2, baselineSelection, os);
  }

  Matrix<bool> BaselineSelect::convert (Table& anttab,
                                        TableExprNode& a1,
                                        TableExprNode& a2,
                                        const string& baselineSelection,
                                        std::ostream& os)
  {
    // Overwrite the error handler to ignore errors for unknown antennas.
    // First construct MSSelection, because it resets the error handler.
    Vector<Int> selectedAnts1;
    Vector<Int> selectedAnts2;
    Matrix<Int> selectedBaselines;
    auto curHandler = MSAntennaParse::thisMSAErrorHandler;
    BaselineSelectErrorHandler errorHandler (os);
    MSAntennaParse::thisMSAErrorHandler = &errorHandler;
    try {
      // Create a table expression representing the selection.
      TableExprNode node = msAntennaGramParseCommand
        (anttab, a1, a2, baselineSelection, 
         selectedAnts1, selectedAnts2, selectedBaselines);
      // Get the antenna numbers.
      Table seltab = node.table()(node);
      Vector<Int> a1 = ROScalarColumn<Int>(seltab, "ANTENNA1").getColumn();
      Vector<Int> a2 = ROScalarColumn<Int>(seltab, "ANTENNA2").getColumn();
      int nant = anttab.nrow();
      Matrix<bool> bl(nant, nant, false);
      for (unsigned int i=0; i<a1.size(); ++i) {
        bl(a1[i], a2[i]) = true;
        bl(a2[i], a1[i]) = true;
      }
      MSAntennaParse::thisMSAErrorHandler = curHandler;
      return bl;
    } catch (const std::exception&) {
      MSAntennaParse::thisMSAErrorHandler = curHandler;
      throw;
    }
  }


  BaselineSelectErrorHandler::~BaselineSelectErrorHandler()
  {}

  void BaselineSelectErrorHandler::reportError (const char* token,
                                                const String message)
  {
    itsStream << message << token << endl;
  }

} // end namespace
