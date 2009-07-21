/***************************************************************************
 *   Copyright (C) 2007-8 by Adriaan Renting, ASTRON                       *
 *   renting@astron.nl                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
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
      void ProcessTimeslot(DataBuffer& InData, DataBuffer& OutData,
                           MsInfo& Info, RunDetails& Details,
                           TimeBuffer& TimeData);
    private:
      /// This will squash every Step Channels starting from Start into one channel, until NChan have been processed.
      /// newWeights will reflect how many values were unflagged. Only unflagged values are retained
      void Squash(casa::Matrix<casa::Complex>& oldData, casa::Matrix<casa::Complex>& newData,
                  casa::Matrix<casa::Bool>& oldFlags, casa::Matrix<casa::Bool>& newFlags,
                  casa::Matrix<casa::Float>& newWeights, int itsNumPolarizations,
                  int Start, int Step, int NChan);

    }; //DataSquasher
  }; //namespace CS1
}; //namespace LOFAR
#endif
