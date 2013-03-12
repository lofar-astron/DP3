//# BaselineSelection.h: Class to handle the baseline selection
//# Copyright (C) 2012
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

#ifndef DPPP_BASELINESELECTION_H
#define DPPP_BASELINESELECTION_H

// @file
// @brief Class to handle the baseline selection

#include <DPPP/DPInfo.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

namespace LOFAR {
  class ParameterSet;
  class ParameterValue;

  namespace DPPP {

    // Class containing a few static functions to parse a baseline selection
    // string.
    class BaselineSelection
    {
    public:
      // Default constructor has no selection.
      BaselineSelection();

      // Construct from the parset using the given prefix.
      // The keys used are:
      // <ul>
      //  <li> baseline: for a baseline selection
      //  <li> corrtype: for correlation selection (auto, cross, or empty)
      //  <li> blrange:  ranges of baseline lengths (in m)
      //  <li> minbl:    minimum baseline length (in m); only if minmax=true
      //  <li> maxbl:    maximum baseline length (in m); only if minmax=true
      // </ul>
      BaselineSelection (const ParameterSet&, const string& prefix,
                         bool minmax=false,
			 const string& defaultCorrType=string());

      // Is there any selection?
      bool hasSelection() const;

      // Show the parameters.
      void show (ostream& os) const;

      // Form the selection matrix telling for each baseline if it is selected.
      // An empty matrix is returned if no selection was made.
      casa::Matrix<bool> apply (const DPInfo& info) const;

    private:
      // Convert the baseline selection string.
      void handleBL (casa::Matrix<bool>& selectBL,
                     const DPInfo& info) const;

      // Handle a vector of baseline specifications.
      casa::Matrix<bool> handleBLVector (const ParameterValue& pvBL,
                                         const casa::Vector<casa::String>&) const;

      // Handle the correlation type selection.
      void handleCorrType (casa::Matrix<bool>& selectBL) const;

      // Handle the baseline length selection.
      void handleLength (casa::Matrix<bool>& selectBL,
                         const DPInfo& info) const;

      //# Data members
      string itsStrBL;
      string itsCorrType;
      vector<double> itsRangeBL;
    };

  } //# end namespace
}

#endif

