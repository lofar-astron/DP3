//# Copyright (C) 2006-8
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
//# @author Adriaan Renting

#ifndef __CS1_PP_RUN_DETAILS_H__
#define __CS1_PP_RUN_DETAILS_H__

/// @file
/// @brief Class to hold parameter settings for the steps in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

#include <casa/Arrays/Vector.h>
#include <vector>
#include <string>

namespace LOFAR
{
  namespace CS1
  {
    /// @ingroup DPPP
    ///
    /// This class is basically a fancy struct of the variables that need to be set
    /// in the parameterset for IDPPP for the various steps.
    /// It has some extra functionality to validate and print the values.
    /// It isn't using setter/getter functions because that seemed overkill for the required use.

    class RunDetails
    {
    public:
      RunDetails();
      ~RunDetails();

      unsigned int Fixed;        ///< BandpassCorrector
      unsigned int FreqWindow;   ///< FrequencyFlagger, MADFlagger
      unsigned int TimeWindow;   ///< ComplexMedianFlagger, MADFlagger
      double       Threshold;    ///< FrequencyFlagger, MADFlagger, ComplexFlagger2
      double       MinThreshold; ///< ComplexMedianFlagger
      double       MaxThreshold; ///< ComplexMedianFlagger, MADFlagger
      unsigned int Algorithm;    ///< FrequencyFlagger
      bool         Existing;     ///< all flaggers
      unsigned int NChan;        ///< DataSquasher
      unsigned int Start;        ///< DataSquasher
      unsigned int Step;         ///< DataSquasher
      bool         Skip;         ///< DataSquasher
      bool         AllColumns;   ///< DataSquasher
      unsigned int TimeStep;     ///< DataSquasher
      unsigned int TileNChan;    ///< Nr of channels per tile for the DATA columns
      unsigned int TileSize;     ///< Tile size (in kbytes) for the DATA columns
      std::string  FlagColumn;   ///< Data column to use in all flaggers
      std::vector<std::string>  DataColumns; ///< Data columns to handle
      casa::Vector<casa::String> AllParms;  ///< All parameters and their values

      bool CheckValues(void);    ///< Method to do some validity checks on the values
      void PrintInfo(void);      ///< Prints all values to cout, mainly for debugging purposes

    protected:
    private:
    }; // class RunDetails
  }; // namespace CS1
}; // namespace LOFAR

#endif // __CS1_PP_RUN_DETAILS_H__
