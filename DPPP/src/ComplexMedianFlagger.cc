/***************************************************************************
 *   Copyright (C) 2006-8 by ASTRON, Adriaan Renting                       *
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

#include <lofar_config.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>
#include <CS1_pp_lib/ComplexMedianFlagger.h>
#include <casa/Quanta/MVEpoch.h>

#include <CS1_pp_lib/MsInfo.h>
#include <CS1_pp_lib/RunDetails.h>
#include <CS1_pp_lib/DataBuffer.h>
#include <CS1_pp_lib/FlaggerStatistics.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  ComplexMedianFlagger::ComplexMedianFlagger  <<<===============
/* initialize some meta data and get the datastorage the right size. */
ComplexMedianFlagger::ComplexMedianFlagger(void):
  NumChannels(0),
  NumPolarizations(0)
{
}

//===============>>>  ComplexMedianFlagger::~ComplexMedianFlagger  <<<===============

ComplexMedianFlagger::~ComplexMedianFlagger()
{
}

//===============>>> ComplexMedianFlagger::FlagTimeslot  <<<===============
/* This function inspects each visibility in a cetain baseline-band
and flags on complexe distance, then determines to flag the entire baseline-band
based on the RMS of the points it didn't flag.*/
int ComplexMedianFlagger::FlagBaselineBand(Matrix<Bool>& Flags,
                                           Cube<Complex>& Data,
                                           int flagCounter,
                                           double FlagThreshold,
                                           int Position,
                                           bool Existing,
                                           int WindowSize)
{
  Vector<Float>  Reals(WindowSize);
  Vector<Float>  Imags(WindowSize);
  int            flagcount = 0;
  for (int i = NumChannels-1; i >= 0; i--)
  {
    bool FlagAllPolarizations = false;
    for (int j = NumPolarizations-1; j >= 0; j--)
    { //we need to loop twice, once to determine FlagAllCorrelations
      if (!FlagAllPolarizations /*&& PolarizationsToCheck[j]*/)
      {
        for (int k = 0; k < WindowSize; k++)
        { //This might be faster in some other way ?
          Reals[k] = Data(j, i, k).real();
          Imags[k] = Data(j, i, k).imag();
        }
        float real           = medianInPlace(Reals, false);
        float imag           = medianInPlace(Imags, false);
        FlagAllPolarizations |= FlagThreshold < abs(Data(j, i, Position)
                                                   - Complex(real, imag));
      }
    }
    for (int j = NumPolarizations-1; j >= 0; j--)
    { //the second loop we set the flags or calculate RMS
      if (FlagAllPolarizations)
      { Flags(j, i) = true;
        flagcount++;
      }
      else
      { if (!Existing) { Flags(j, i) = false;}
      }
    }
  }
  return flagCounter + flagcount;
}
//===============>>> ComplexMedianFlagger::FlagBaseline  <<<===============
/* This function iterates over baseline and band and uses FlagBaselineBand() to determine
   for each one if it needs to be flagged. It treats the autocorrelations separately,
   to detect entire malfunctioning telescopes. Finally it writes the flags.
*/
void ComplexMedianFlagger::ProcessTimeslot(DataBuffer& data,
                                           MsInfo& info,
                                           RunDetails& details,
                                           FlaggerStatistics& stats)
{
  //Data.Position is the last filled timeslot, the middle is 1/2 a window behind it
  int pos = (data.Position + (data.WindowSize+1)/2) % data.WindowSize;
  NumChannels      = info.NumChannels;
  NumPolarizations = info.NumPolarizations;
  int index        = 0;
  double treshold  = 0.0;
  Matrix<Bool> flags;

  for (int i = 0; i < info.NumBands; i++)
  {
    for(int j = 0; j < info.NumAntennae; j++)
    {
      for(int k = j; k < info.NumAntennae; k++)
      {
        index    = i * info.NumPairs + info.BaselineIndex[baseline_t(j, k)];
        treshold = (details.MinThreshold + (details.MaxThreshold - details.MinThreshold)
                    * info.BaselineLengths[info.BaselineIndex[baseline_t(j, k)]]
                    / info.MaxBaselineLength
                   ) * info.NoiseLevel;
        flags.reference(data.Flags[index].xyPlane(pos));
//        if ((BaselineLengths[BaselineIndex[pairii(j, k)]] < 3000000))//radius of the Earth in meters? WSRT sometimes has fake telescopes at 3854243 m
        stats(i, j, k) = FlagBaselineBand(flags,
                                          data.Data[index],
                                          stats(i,j,k),
                                          treshold, pos,
                                          details.Existing,
                                          data.WindowSize);
      }
    }
  }
}

//===============>>> ComplexMedianFlagger  <<<===============
