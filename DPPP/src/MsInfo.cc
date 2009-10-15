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
#include <casa/BasicMath/Math.h>
#include <casa/Arrays.h>

#include <iostream>

#include <DPPP/MsInfo.h>

#include <Common/LofarLogger.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  Ms_Info::Ms_Info  <<<===============

MsInfo::MsInfo():
  NumSamples(0),
  NumAntennae(0),
  NumFields(0),
  NumBands(0),
  NumChannels(0),
  NumPolarizations(0),
  NumPairs(0),
  NumTimeslots(0),
  NoiseLevel(0.0),
  MaxBaselineLength(0.0)
{}

MsInfo::MsInfo(const MeasurementSet& MS, const Table& orderedMainTable,
               bool checkRegularity):
  NumSamples(0),
  NumAntennae(0),
  NumFields(0),
  NumBands(0),
  NumChannels(0),
  NumPolarizations(0),
  NumPairs(0),
  NumTimeslots(0),
  NoiseLevel(0.0),
  MaxBaselineLength(0.0)
{
  //Number of samples
  NumSamples                       = orderedMainTable.nrow();
  //Number of Fields
  MSField fields                   = MS.field();
  NumFields                        = fields.nrow();

  //Number of Antennae
  MSAntenna antennae               = MS.antenna();
  NumAntennae                      = antennae.nrow();

  //Antenna Names
  ROScalarColumn<String>           ANT_NAME_col(antennae, "NAME");
  Vector<String>         ant_names = ANT_NAME_col.getColumn();
  ant_names.tovector(AntennaNames);

  //Number of channels in the Band
  MSSpectralWindow spectral_window = MS.spectralWindow();
  ROScalarColumn<Int>              NUM_CHAN_col(spectral_window, "NUM_CHAN");
  NumChannels                      = NUM_CHAN_col(0);

  //Number of polarizations
  MSPolarization      polarization = MS.polarization();
  ROScalarColumn<Int>              NUM_CORR_col(polarization, "NUM_CORR");
  NumPolarizations                 = NUM_CORR_col(0);
  ROArrayColumn<Int>               CORR_TYPE_col(polarization, "CORR_TYPE");
  Polarizations.resize(NumPolarizations);
  CORR_TYPE_col.get(0, Polarizations);

  //calculate theoretical noise level
  // MS might be empty (e.g. when called for output MS).
  Double exposure                  = 1;
  if (orderedMainTable.nrow() > 0)
  {
    ROScalarColumn<Double>         EXPOSURE_col(orderedMainTable, "EXPOSURE");
    exposure                       = EXPOSURE_col(0);
  }

  ROScalarColumn<Double>           TOTAL_BANDWIDTH_col(spectral_window, "TOTAL_BANDWIDTH");
  Double bandwidth                 = TOTAL_BANDWIDTH_col(0) / NumChannels;

  NoiseLevel                       = 1.0 / sqrt(bandwidth * exposure);

  //calculate number of timeslots
  if (orderedMainTable.nrow() > 0) {
    ROScalarColumn<Double>         INTERVAL_col(orderedMainTable, "INTERVAL");
    Double Interval                = INTERVAL_col(0);
    //Number of timeslots
    ROScalarColumn<Double>         TIME_col(orderedMainTable, "TIME");
    Double firstdate               = TIME_col(0);
    Double lastdate                = TIME_col(NumSamples-1);
    NumTimeslots                   = 1 + int((lastdate-firstdate)/Interval + 0.5);
  }

  //calculate number of baselines.
  // It assumes auto-correlations are present.
  NumPairs = (NumAntennae) * (NumAntennae + 1) / 2; //Triangular numbers formula

  //calculate number of Bands
  if (NumSamples)
  { NumBands                       = NumSamples / (NumPairs * NumTimeslots);
  }
  else
  { NumBands                       = spectral_window.nrow();
  }
  ///  if (checkRegularity) {
  ///    ASSERTSTR (NumBands*NumTimeslots*NumPairs == NumSamples,
  ///               "The MS cannot be handled by DPPP; it should contain:\n"
  ///               " - cross and auto-correlations for all antennae in ANTENNA table\n"
  ///               " - no missing time slots\n"
  ///               " - the same interval length for each time slot");
///  }
  PairsIndex.resize(NumPairs);

  int index = 0;
  for (int i = 0; i < NumAntennae; i++)
  { for(int j = i; j < NumAntennae; j++)
    { PairsIndex[index]               = baseline_t(i, j);
      BaselineIndex[baseline_t(i, j)] = index++;
    }
  }

  ComputeBaselineLengths(MS);
}

//===============>>>  Ms_Info::~Ms_Info  <<<===============

MsInfo::~MsInfo()
{
}

//===============>>> MS_Info::PrintInfo  <<<===============

void MsInfo::PrintInfo(void)
{
  std::cout << "NumSamples:       " << NumSamples << std::endl;
  std::cout << "NumFields:        " << NumFields << std::endl;
  std::cout << "NumAntennae:      " << NumAntennae << std::endl;
  std::cout << "NumChannels:      " << NumChannels << std::endl;
  std::cout << "NumPolarizations: " << NumPolarizations << std::endl;
  if(NoiseLevel < 1000)
  { std::cout << "NoiseLevel:       " << NoiseLevel << std::endl;
  }
  else
  { std::cout << "NoiseLevel:       " << "undetermined" << std::endl;
  }
  std::cout << "Numtimeslots:     " << NumTimeslots << std::endl;
  std::cout << "NumPairs:         " << NumPairs << std::endl;
  std::cout << "NumBands:         " << NumBands << std::endl;
}

//===============>>> MS_Info::ComputeBaselineLengths  <<<===============
/* compute baseline lengths, and determine the longest one.*/

void MsInfo::ComputeBaselineLengths(const casa::MeasurementSet& MS)
{
  BaselineLengths.resize(NumPairs);
  //Antenna positions
  MSAntenna antenna     = MS.antenna();
  ROArrayColumn<Double>  position(antenna, "POSITION");
  for (int i = 0; i < NumAntennae; i++ )
  {
    for (int j = i; j < NumAntennae; j++)
    {
      Vector<Double> p(position(i) - position(j));
      double temp   = sqrt(p(0)*p(0) + p(1)*p(1) + p(2)*p(2));
      BaselineLengths[BaselineIndex[baseline_t(i, j)]] = temp;
      if (temp > MaxBaselineLength && temp < 3000000) //radius of the Earth in meters? WSRT sometimes has fake telescopes at 3854243
      { MaxBaselineLength = temp;                     // non-existent antenna's can have position (0,0,0)
      }
    }
  }
}

//===============>>> MsInfo  <<<===============

