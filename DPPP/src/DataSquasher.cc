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

/*
//===============>>>  DataSquasher::add  <<<===============

void DataSquasher::add (Matrix<Complex>& sumData,
                        Matrix<Complex>& allData,
                        Matrix<Int>& sumNPoint,
                        Matrix<Float>& sumWeight,
                        const Matrix<Complex>& inData,
                        const Matrix<Bool>& inFlag,
                        const Matrix<Float>& inWeight,
                        int npol, int stchan, int step, int nchan)
{
  const Complex* pinData   = inData.data()   + stchan*npol;
  const Bool*    pinFlag   = inFlag.data()   + stchan*npol;
  const Float*   pinWeight = inWeight.data() + stchan*npol;
  Complex* psumData   = sumData.data();
  Complex* pallData   = allData.data();
  Float*   psumWeight = sumWeight.data();
  Int*     psumNP     = sumNPoint.data();
  for (int i=0; i<nchan; ++i) {
    for (int j=0; j<npol; ++j) {
      pallData[j] += *pinData;
      if (! *pinFlag) {
        psumData[j]   += *pinData * *pinWeight;
        psumWeight[j] += *pinWeight;
        psumNP[j]     += 1;
      }
      ++pinData;
      ++pinFlag;
      ++pinWeight;
    }
    if (i%step == 0  ||  i == nchan-1) {
      pallData   += npol;
      psumData   += npol;
      psumWeight += npol;
      psumNP     += npol;
    }
  }
}

//===============>>>  DataSquasher::average  <<<===============

Matrix<Bool> DataSquasher::average (Matrix<Complex>& sumData,
                                    Matrix<Float>& sumWeight,
                                    const Matrix<Complex>& allData,
                                    const Matrix<Int>& sumNPoint)
{
  Matrix<Bool> flags(sumData.shape());
  Complex* pData   = sumData.data();
  Float*   pWeight = sumWeight.data();
  Bool*    pFlag   = flags.data();
  const Complex* pallData = allData.data();
  const Int*     pNP      = sumNPoint.data();
  int n = sumData.size();
  for (int i=0; i<n; ++i) {
    if (pNP[i] == 0  ||  pWeight[i] == 0) {
      pData[i]   = pallData[i];
      pWeight[i] = 0;
      pFlag[i]   = True;
    } else {
      pData[i]   /= pWeight[i];
      pWeight[i] /= pNP[i];
      pFlag[i]   = False;
    }
  }
  return flags;
}
*/

//===============>>>  DataSquasher::Squash  <<<===============

void DataSquasher::Squash(vector<Matrix<Complex> >& oldData,
                          vector<Matrix<Complex> >& newData,
                          Matrix<Bool>& oldFlags, Matrix<Bool>& newFlags,
                          Matrix<Float>& oldWeights, Matrix<Float>& newWeights,
                          int itsNumPolarizations,
                          int Start, int Step, int NChan)
{ //We only add to weight as it can have multiple timesteps integrated
  int nDataCol   = oldData.size();
  bool flagnew = true;
  Vector<Complex> values(nDataCol);
  Vector<Complex> allvalues(nDataCol);
  Vector<Complex> tmpvalues(nDataCol);
  // Use raw pointer as it is faster.
  Complex* tmpval = tmpvalues.data();
  Float weight = 0;
  for (int i = 0; i < itsNumPolarizations; i++)
  {
    int incounter  = 0;
    int outcounter = 0;
    int nstep = 0;
    while (incounter < NChan)
    {
      ++nstep;
      for (int j=0; j<nDataCol; ++j) {
        tmpval[j] = oldData[j](i, Start + incounter);
      }
      allvalues += tmpvalues;
      if (!oldFlags(i, Start + incounter)) {
        values += tmpvalues;
        weight += 1.0; //should be += old Weight?
        flagnew = false;
      }
      incounter++;
      if (incounter == NChan-1  ||  incounter%Step == 0)
      {
        if (flagnew) {
          // Everything is flagged
          // take all values and set weight to 1.0
          values = allvalues;
          newWeights(i, outcounter) += 1.0; //should be += old Weight * nstep?
        } else {
          //Not everything is flagged
          values /= weight;
          newWeights(i, outcounter) += abs(weight) / nstep;
        }
        for (int j=0; j<nDataCol; ++j) {
          newData[j](i, outcounter) += values[j];
        }
        newFlags(i, outcounter) = flagnew || newFlags(i, outcounter);
        allvalues = 0.;
        values    = 0.;
        weight    = 0.;
        outcounter++;
        flagnew = true;
        nstep = 0;
      }
    }
  }
}

//===============>>>  DataSquasher::ProcessTimeslot  <<<===============

void DataSquasher::ProcessTimeslot(const DataBuffer& InData,
                                   DataBuffer& OutData,
                                   MsInfo& Info, const RunDetails& Details,
                                   const TimeBuffer& TimeData)
{
  const int nDataCol = InData.dataColumns().size();
  //Data.Position is the last filled timeslot, we need to process the one just in front of it.
  int inpos  = (InData.Position + 1) % InData.WindowSize;
  int outpos = 0;
  vector<Matrix<Complex> > myOldData(nDataCol);
  vector<Matrix<Complex> > myNewData(nDataCol);
  Matrix<Bool>    myOldFlags;
  Matrix<Bool>    myNewFlags;
  Matrix<Float>   OldWeights;
  Matrix<Float>   NewWeights;

  for (int i = 0; i < Info.NumBands; i++)
  {
    for(int j = 0; j < Info.NumAntennae; j++)
    {
      for(int k = j; k < Info.NumAntennae; k++)
      {
        int inx = Info.getBaselineIndex(j, k);
        if (inx >= 0) {
          int index = i * Info.NumPairs + inx;
          for (int ic=0; ic<nDataCol; ++ic) {
            myOldData[ic].reference(InData.Data[ic][index].xyPlane(inpos));
            myNewData[ic].reference(OutData.Data[ic][index].xyPlane(outpos));
          }
          myOldFlags.reference(InData.Flags[index].xyPlane(inpos));
          myNewFlags.reference(OutData.Flags[index].xyPlane(outpos));
          OldWeights.reference(InData.Weights[index].xyPlane(inpos));
          NewWeights.reference(OutData.Weights[index].xyPlane(outpos));
          if (TimeData.Time[index].size() == 1) {
            for (int ic=0; ic<nDataCol; ++ic) {
              myNewData[ic] = 0.0;
            }
            myNewFlags = false;
            NewWeights = 0.0;
          }
          Squash(myOldData, myNewData, myOldFlags, myNewFlags, OldWeights, NewWeights,
                 Info.NumPolarizations, Details.Start, Details.Step, Details.NChan);
        }
      }
    }
  }
}

//===============>>>  DataSquasher  <<<===============
