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
#include <DPPP/FrequencyFlagger.h>
#include <casa/Quanta/MVEpoch.h>

#include <DPPP/MsInfo.h>
#include <DPPP/RunDetails.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/FlaggerStatistics.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  FrequencyFlagger::FrequencyFlagger  <<<===============
/* initialize some meta data and get the datastorage the right size. */
FrequencyFlagger::FrequencyFlagger():
  NumChannels(0),
  NumPolarizations(0)
{
}

//===============>>>  FrequencyFlagger::~FrequencyFlagger  <<<===============

FrequencyFlagger::~FrequencyFlagger()
{
}

//===============>>> FrequencyFlagger::FlagTimeslot  <<<===============
/* This function inspects each visibility in a cetain baseline-band
and flags on complexe distance, then determines to flag the entire baseline-band
based on the RMS of the points it didn't flag.*/
int FrequencyFlagger::FlagBaselineBand(casa::Matrix<casa::Bool>& Flags,
                                       const casa::Matrix<casa::Complex>& Data,
                                       int flagCounter,
                                       double FlagThreshold,
                                       bool Existing,
                                       int Algorithm)
{
  vector<double> MS(NumPolarizations, 0.0);
  vector<double> RMS(NumPolarizations);
  Matrix<double> MedianArray(NumPolarizations, NumChannels);
  vector<int>    RMSCounter(NumPolarizations, 0);
  for (int j = NumPolarizations-1; j >= 0; j--)
  {
    for (int i = NumChannels-1; i >= 0; i--) //calculata RMS of unflagged datapoints
    {
      if (Algorithm == 1)
      { if (!Existing || !Flags(j, i))
        { double temp = pow(abs(Data(j, i)), 2);
          if (!isNaN(temp))
          {
            MS[j] += temp;
            RMSCounter[j] += 1;
          }
        }
      }
      else
      { if (!Existing || !Flags(j, i))
        { double temp = abs(Data(i, j));
          if (!isNaN(temp))
          {
            MedianArray(j, RMSCounter[j]) = temp;
            RMSCounter[j] += 1;
          }
        }
      }
    }
    if (RMSCounter[j])
    { if (Algorithm == 1)
      { RMS[j] = sqrt(MS[j] /RMSCounter[j]);
      }
      else
      { RMS[j] = medianInPlace(MedianArray.row(j).operator()(Slice(0, RMSCounter[j]-1,1)));
      }
      for (int i = NumChannels-1; i >= 0; i--)
      {
        if (!Existing || !Flags(j, i))
        { double temp = abs(Data(j, i));
          bool flag   = isNaN(temp) || RMS[j] * FlagThreshold < temp;
          if (flag)
          { flagCounter++;
          }
          Flags(j, i) = flag || (Existing && Flags(j, i));
        }
      }
    }
    else
    {
      flagCounter += NumChannels;
      Flags.row(j) = true;
    }
  }
  return flagCounter;
}
//===============>>> FrequencyFlagger::FlagBaseline  <<<===============
/* This function iterates over baseline and band and uses FlagBaselineBand() to determine
  for each one if it needs to be flagged. It treats the autocorrelations separately,
  to detect entire malfunctioning telescopes. Finally it writes the flags.
*/
void FrequencyFlagger::ProcessTimeslot(DataBuffer& data,
                             MsInfo& info,
                             RunDetails& details,
                             FlaggerStatistics& stats)
{
  //Data.Position is the last filled timeslot, the middle is 1/2 a window behind it
  int pos = (data.Position + (data.WindowSize+1)/2) % data.WindowSize;
  NumChannels      = info.NumChannels;
  NumPolarizations = info.NumPolarizations;
  int index        = 0;
  Matrix<Bool>     flags;
  Matrix<Complex>  tempdata;

  for (int i = 0; i < info.NumBands; i++)
  {
    for(int j = 0; j < info.NumAntennae; j++)
    {
      for(int k = j; k < info.NumAntennae; k++)
      {
        index = i * info.NumPairs + info.BaselineIndex[baseline_t(j, k)];
        flags.reference(data.Flags[index].xyPlane(pos));
        tempdata.reference(data.Data[index].xyPlane(pos));
//        if ((BaselineLengths[BaselineIndex[pairii(j, k)]] < 3000000))//radius of the Earth in meters? WSRT sometimes has fake telescopes at 3854243 m
        stats(i, j, k) = FlagBaselineBand(flags,
                                          tempdata,
                                          stats(i,j,k),
                                          details.Threshold,
                                          details.Existing,
                                          details.Algorithm);
      }
    }
  }
}

//===============>>> FrequencyFlagger  <<<===============
