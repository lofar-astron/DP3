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

#ifndef LOFARDATASQUASHER_H
#define LOFARDATASQUASHER_H

/**
@author Adriaan Renting
*/
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <ms/MeasurementSets.h>
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>

/// @file
/// @brief Class to hold code for DataSquasher step in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    class DataBuffer;
    class TimeBuffer;
    class MsInfo;
    class RunDetails;

    class DataSquasher
    {
    public:
      DataSquasher(void);
      ~DataSquasher(void);

      /// All processing of one integration time happens in one go.
      void ProcessTimeslot(const DataBuffer& InData, DataBuffer& OutData,
                           MsInfo& Info, const RunDetails& Details,
                           const TimeBuffer& TimeData);
    private:
      /// This will squash every Step Channels starting from Start into one channel, until NChan have been processed.
      /// newWeights will reflect how many values were unflagged. Only unflagged values are retained
      void Squash(casa::Matrix<casa::Complex>& oldData, casa::Matrix<casa::Complex>& newData,
                  casa::Matrix<casa::Bool>& oldFlags, casa::Matrix<casa::Bool>& newFlags,
                  casa::Matrix<casa::Float>& oldWeights, casa::Matrix<casa::Float>& newWeights,
                  int itsNumPolarizations,
                  int Start, int Step, int NChan);

      // Add the input buffers to the sum buffers.
      // It takes flags and weights into account.
      void add (casa::Matrix<casa::Complex>& sumData,
                casa::Matrix<casa::Complex>& allData,
                casa::Matrix<casa::Int>& sumNPoint,
                casa::Matrix<casa::Float>& sumWeight,
                const casa::Matrix<casa::Complex>& inData,
                const casa::Matrix<casa::Bool>& inFlag,
                const casa::Matrix<casa::Float>& inWeight,
                int npol, int stchan, int step, int nchan);

      casa::Matrix<casa::Bool> average (casa::Matrix<casa::Complex>& sumData,
                                        casa::Matrix<casa::Float>& sumWeight,
                                        const casa::Matrix<casa::Complex>& allData,
                                        const casa::Matrix<casa::Int>& sumNPoint);

    }; //DataSquasher
  }; //namespace CS1
}; //namespace LOFAR
#endif
