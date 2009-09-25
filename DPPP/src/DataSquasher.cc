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

#include <lofar_config.h>
#include <casa/Arrays/ArrayMath.h>
#include <DPPP/DataSquasher.h>
#include <DPPP/MsInfo.h>
#include <DPPP/RunDetails.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/TimeBuffer.h>


using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  DataSquasher::DataSquasher  <<<===============

DataSquasher::DataSquasher(void)
{
}

//===============>>>  DataSquasher::~DataSquasher  <<<===============

DataSquasher::~DataSquasher(void)
{
}

//===============>>>  DataSquasher::Squash  <<<===============

void DataSquasher::Squash(Matrix<Complex>& oldData, Matrix<Complex>& newData,
                          Matrix<Bool>& oldFlags, Matrix<Bool>& newFlags,
                          Matrix<Float>& newWeights, int itsNumPolarizations,
                          int Start, int Step, int NChan)
{ //We only add to weight as it can have multiple timesteps integrated
  int incounter  = 0;
  int outcounter = 0;
  bool flagnew   = true;
  Vector<Complex> values(itsNumPolarizations, 0);
  Vector<Complex> allvalues(itsNumPolarizations, 0);
  Vector<Float>   weights(itsNumPolarizations, 0);
  while (incounter < NChan)
  {
    for (int i = 0; i < itsNumPolarizations; i++)
    {
      allvalues(i) += oldData(i, Start + incounter);
      if (!oldFlags(i, Start + incounter))
      { //existing weight <> 1 is not handled here, maybe sometime in the future?
        //On new data WEIGHT_SPECTRUM does not exist, only WEIGHT
        values(i) += oldData(i, Start + incounter);
        weights(i) += 1.0; //should be += old Weight?
        flagnew = false;
      }
    }
    incounter++;
    if ((incounter) % Step == 0)
    {
      for (int i = 0; i < itsNumPolarizations; i++)
      { if (flagnew) //Everything is flagged
        { values(i) = allvalues(i); //we take all values and put weight at 1.0
          newWeights(i, outcounter) += 1.0; //should be += old Weight * Step?
        }
        else //Not everything is flagged
        { values(i) = values(i) / weights(i); //We only take unflagged values and adjust weight betweem 1/step and 1.0
          newWeights(i, outcounter) += abs(weights(i)) / Step;
        }
      }
      newData.column(outcounter)  = newData.column(outcounter) + values;
      newFlags.column(outcounter) = flagnew || newFlags.column(outcounter);
      allvalues = 0;
      values    = 0;
      weights   = 0;
      outcounter++;
      flagnew = true;
    }
  }
}

//===============>>>  DataSquasher::ProcessTimeslot  <<<===============

void DataSquasher::ProcessTimeslot(const DataBuffer& InData, DataBuffer& OutData,
                                   MsInfo& Info, const RunDetails& Details,
                                   const TimeBuffer& TimeData)
{
  //Data.Position is the last filled timeslot, we need to process the one just in front of it.
  int inpos  = (InData.Position + 1) % InData.WindowSize;
  int outpos = 0;
  bool columns = InData.ModelData.size() > 0;
  Matrix<Complex> myOldData;
  Matrix<Complex> myNewData;
  Matrix<Bool>    myOldFlags;
  Matrix<Bool>    myNewFlags;
  Matrix<Float>   NewWeights;

  for (int i = 0; i < Info.NumBands; i++)
  {
    for(int j = 0; j < Info.NumAntennae; j++)
    {
      for(int k = j; k < Info.NumAntennae; k++)
      {
        int index = i * Info.NumPairs + Info.BaselineIndex[baseline_t(j, k)];

        myOldData.reference(InData.Data[index].xyPlane(inpos));
        myNewData.reference(OutData.Data[index].xyPlane(outpos));
        myOldFlags.reference(InData.Flags[index].xyPlane(inpos));
        myNewFlags.reference(OutData.Flags[index].xyPlane(outpos));
        NewWeights.reference(OutData.Weights[index].xyPlane(outpos));
        if (TimeData.Time[index].size() == 1)
        {
          myNewData  = 0.0;
          myNewFlags = false;
          NewWeights = 0.0;
        }
        Squash(myOldData, myNewData, myOldFlags, myNewFlags, NewWeights,
               Info.NumPolarizations, Details.Start, Details.Step, Details.NChan);
        if (columns)
        {
          myOldData.reference(InData.ModelData[index].xyPlane(inpos));
          myNewData.reference(OutData.ModelData[index].xyPlane(outpos));
          Squash(myOldData, myNewData, myOldFlags, myNewFlags, NewWeights,
                 Info.NumPolarizations, Details.Start, Details.Step, Details.NChan);

          myOldData.reference(InData.CorrectedData[index].xyPlane(inpos));
          myNewData.reference(OutData.CorrectedData[index].xyPlane(outpos));
          Squash(myOldData, myNewData, myOldFlags, myNewFlags, NewWeights,
                 Info.NumPolarizations, Details.Start, Details.Step, Details.NChan);
        }
      }
    }
  }
}

//===============>>>  DataSquasher  <<<===============
