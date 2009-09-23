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
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>
#include <DPPP/MADFlagger.h>
#include <casa/Quanta/MVEpoch.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Utilities/GenSort.h>

#include <DPPP/MsInfo.h>
#include <DPPP/RunDetails.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/FlaggerStatistics.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  MADFlagger::MADFlagger  <<<===============
/* initialize some meta data and get the datastorage the right size. */
MADFlagger::MADFlagger():
  NumChannels(0),
  NumPolarizations(0)
{
}

//===============>>>  MADFlagger::~MADFlagger  <<<===============

MADFlagger::~MADFlagger()
{
}

//===============>>> MADFlagger::ComputeThreshold  <<<============
/*Compute Thresholds */

Cube<Bool> MADFlagger::MakeMask(const Cube<Float>& Values, const Cube<Bool>& Flags,
                                double MaxLevel, int WindowSize,
                                int Position, bool Existing)
{
  Cube<Bool> Mask(Values.shape(), false);
  for (int k = WindowSize-1; k >= 0; k--)
  {
    for (int i = NumChannels-1; i >= 0; i--)
    {
      for (int j = NumPolarizations-1; j >= 0; j--)
      { Mask(j, i, k) = !isFinite(Values(j, i, k))
                      || MaxLevel && Values(j, i, Position) > MaxLevel
                      || (Existing && Flags(j, i, k));
      }
    }
  }
  return Mask;
}

//===============>>> MADFlagger::ComputeThreshold  <<<============
/*Compute Thresholds */
void MADFlagger::ComputeThreshold(const Cube<Float>& Values, const Cube<Bool>& Mask,
                                  int TWindowSize, int FWindowSize,
                                  int TimePos, int ChanPos, int PolPos,
                                  float& Z1, float& Z2, Vector<Float>& Medians)
{
  int tempF = 0;
  int tempT = 0;
  int F = FWindowSize/2;
  int T = TWindowSize/2;
  int index = 0;
  for (int i = -T; i <= T; i++)
  {
    for (int j = -F; j <= F; j++)
    {
      tempF = ((ChanPos + j < 0 || ChanPos + j >= NumChannels) ? ChanPos-j : ChanPos+j); //have the channels wrap back upon themselves.
      tempT = (TimePos + i + TWindowSize) % TWindowSize;
      if (!Mask(PolPos, tempF, tempT))
      { Medians(index++) = Values(PolPos, tempF, tempT); //Fill temp buffer
      }
    }
  }
  if (index > 0)
  {
    Float* mediansPointer = Medians.data();
    Z1 = GenSort<Float>::kthLargest(mediansPointer, index, index/2);      // Median Vt = Z
    Medians -= Z1;
    Medians = abs(Medians);
    Z2 = GenSort<Float>::kthLargest(mediansPointer, index, index/2); // Median abs(Vt - Z) = Z'
  }
  else //If there are NaN in the data, then what?
  { Z1 = -1.0;
    Z2 = 0.0;
  }
}

//===============>>> MADFlagger::FlagTimeslot  <<<===============
/* This function inspects each visibility in a cetain baseline-band
and flags on complexe distance, then determines to flag the entire baseline-band
based on the RMS of the points it didn't flag.*/
int MADFlagger::FlagBaselineBand(Cube<Bool>& Flags,
                                 const Cube<Float>& Data,
                                 const Cube<Bool>& Mask,
                                 int flagCounter,
                                 double Threshold,
                                 int Position, bool Existing,
                                 int TWindowSize, int FWindowSize)
{
  float Z1             = 0.0;
  float Z2             = 0.0;
  int    flagcount     = 0;
  double MAD           = 1.4826;
  Vector<Float> Medians(TWindowSize * FWindowSize);
  for (int i = NumChannels-1; i >= 0; i--)
  {
    bool FlagAllPolarizations = false;
    for (int j = NumPolarizations-1; j >= 0; j--)
    { //we need to loop twice, once to determine FlagAllPolarizations, then to set the flags
      if (!FlagAllPolarizations /*&& PolarizationsToCheck[j]*/)
      {
        if (Mask(j,i,Position))
        { FlagAllPolarizations = true;
        }
        else
        { ComputeThreshold(Data, Mask, TWindowSize, FWindowSize, Position, i, j, Z1, Z2, Medians);
          FlagAllPolarizations |= (Threshold * Z2 * MAD) < abs(Data(j, i, Position) - Z1);
        }
      }
    }
    for (int j = NumPolarizations-1; j >= 0; j--)
    { //the second loop we set the flags or calculate RMS
      if (FlagAllPolarizations)
      { Flags(j, i, Position) = true;
        flagcount++;
      }
      else
      { if (!Existing) { Flags(j, i, Position) = false;}
      }
    }
  }
  return flagCounter + flagcount;
}
//===============>>> MADFlagger::FlagBaseline  <<<===============
/* This function iterates over baseline and band and uses FlagBaselineBand() to determine
   for each one if it needs to be flagged. It treats the autocorrelations separately,
   to detect entire malfunctioning telescopes. Finally it writes the flags.
*/
void MADFlagger::ProcessTimeslot(DataBuffer& data,
                                           MsInfo& info,
                                           RunDetails& details,
                                           FlaggerStatistics& stats)
{
  //Data.Position is the last filled timeslot, the middle is 1/2 a window behind it
  int pos = (data.Position + (data.WindowSize+1)/2) % data.WindowSize;
  NumChannels      = info.NumChannels;
  NumPolarizations = info.NumPolarizations;
  int index        = 0;
  for (int i = 0; i < info.NumBands; i++)
  {
    for(int j = 0; j < info.NumAntennae; j++)
    {
      for(int k = j; k < info.NumAntennae; k++)
      {
        index    = i * info.NumPairs + info.BaselineIndex[baseline_t(j, k)];
//        if ((BaselineLengths[BaselineIndex[pairii(j, k)]] < 3000000))//radius of the Earth in meters? WSRT sometimes has fake telescopes at 3854243 m
        Cube<Float> RealData = amplitude(data.Data[index]);
        Cube<Bool> Mask      = MakeMask(RealData, data.Flags[index], details.MaxThreshold,
                                        data.WindowSize, pos, details.Existing);
        stats(i, j, k) = FlagBaselineBand(data.Flags[index],
                                          RealData, Mask,
                                          stats(i,j,k),
                                          details.Threshold,
                                          pos,
                                          details.Existing,
                                          details.TimeWindow,
                                          details.FreqWindow);
      }
    }
  }
}

//===============>>> MADFlagger  <<<===============
