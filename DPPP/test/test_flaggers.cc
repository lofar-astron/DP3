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

#include <DPPP/MsInfo.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/RunDetails.h>
#include <DPPP/FlaggerStatistics.h>
#include <DPPP/MADFlagger.h>
#include <DPPP/ComplexMedianFlagger.h>
#include <DPPP/FrequencyFlagger.h>
#include <casa/BasicMath.h>

using namespace LOFAR::CS1;
using namespace casa;
#define SIZE 7;

int main(int argc, char *argv[])
{
  MsInfo  info("TEST.MS");
  info.NumSamples       = 1;
  info.NumAntennae      = 1;
  info.NumFields        = 1;
  info.NumBands         = 1;
  info.NumChannels      = 64;
  info.NumPolarizations = 4;
  info.NumPairs         = 1;
  info.NumTimeslots     = SIZE;
  info.NoiseLevel       = 0.001;
  info.AntennaNames.push_back("Test1");
  info.Polarizations.resize(4);
  info.Polarizations(0) = 9;  // XX, see DataBuffer.cc
  info.Polarizations(0) = 10; // XY, see DataBuffer.cc
  info.Polarizations(0) = 11; // YX, see DataBuffer.cc
  info.Polarizations(0) = 12; // YY, see DataBuffer.cc
  info.MaxBaselineLength = 30000.0;
  info.PairsIndex.push_back(baseline_t(0, 0));
  info.BaselineIndex[baseline_t(0,0)] = 0;
  info.BaselineLengths.push_back(30000.0);

  RunDetails details;
  details.FreqWindow = 9;
  details.TimeWindow = SIZE;
  details.Existing = false;
  details.Threshold = 4.5;
  DataBuffer data(&info, details.TimeWindow, false);

  ACG gen(11, 20);
  Normal norm(&gen, 0.0, 0.00001);
  for (int i= 0; i < info.NumChannels; i ++)
  {
    for (unsigned int j = 0; j < details.TimeWindow; j++)
    {
    for (int k=0; k < info.NumPolarizations;k++)
      { Complex c(norm() *(10 + j) / 10, norm());
        data.Data[0](k, i, j) = c;
      }
    }
  }

  data.Data[0](0, 30, 0) = Complex(1.0, 0.0);
  data.Data[0](0, 31, 1) = Complex(2.0, 0.0);
  data.Data[0](0, 32, 2) = Complex(4.0, 0.0);
  data.Data[0](0, 33, 3) = Complex(8.0, 0.0);
  data.Data[0](0, 34, 4) = Complex(16.0, 0.0);
  data.Data[0](0, 35, 5) = Complex(32.0, 0.0);
  data.Data[0](0, 36, 6) = Complex(64.0, 0.0);
  data.Data[0](0, 37, 5) = Complex(1.0, 0.0);
  data.Data[0](0, 38, 4) = Complex(0.5, 0.0);
  data.Data[0](0, 39, 3) = Complex(0.25, 0.0);
  data.Data[0](0, 40, 2) = Complex(0.125, 0.0);
  data.Data[0](0, 41, 1) = Complex(0.0625, 0.0);
  data.Data[0](0, 42, 0) = Complex(0.03125, 0.0);
  data.Data[0](0, 43, 1) = Complex(0.015625, 0.0);
  data.Data[0](0, 44, 2) = Complex(0.0078125, 0.0);
  data.Data[0](0, 45, 3) = Complex(0.00390625, 0.0);
  data.Data[0](0, 46, 4) = Complex(0.001953125, 0.0);
  data.Data[0](0, 47, 5) = Complex(0.0009765625, 0.0);
  data.Data[0](0, 48, 6) = Complex(0.00048828125, 0.0);

  for (int i= 0; i < info.NumChannels; i ++)
  { for (unsigned int j = 0; j < details.TimeWindow; j++)
    { std::cout << data.Data[0](0, i, j) << " ";
    }
    std::cout << std::endl;
  }
  FlaggerStatistics stats(info);
  MADFlagger flagger;
  for (unsigned int i = 0; i < details.TimeWindow; i++)
  {
    data.Position = i;
    flagger.ProcessTimeslot(data, info, details, stats);
  }
  for (int i= 0; i < info.NumChannels; i ++)
  { for (unsigned int j = 0; j < details.TimeWindow; j++)
    { std::cout << data.Flags[0](0, i, j) << " ";
    }
    std::cout << std::endl;
  }
  stats.PrintStatistics(std::cout);
}
