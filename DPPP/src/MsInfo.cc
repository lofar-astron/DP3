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
#include <tables/Tables/TableIter.h>

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

MsInfo::MsInfo(const MeasurementSet& MS, const Table& orderedMainTable):
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

  //Max number of Antennae
  MSAntenna antennae               = MS.antenna();

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
  //Get the first time slot and retrieve the unique antennae from it.
  Block<String> ms_iteration_variables(1);
  ms_iteration_variables[0] = "TIME";
  TableIterator iter(orderedMainTable, ms_iteration_variables,
                     TableIterator::Ascending, TableIterator::NoSort);
  Block<String> sort_cols(2);
  sort_cols[0] = "ANTENNA1";
  sort_cols[1] = "ANTENNA2";
  Table sortab = iter.table().sort (sort_cols, Sort::Ascending,
                                    Sort::QuickSort+Sort::NoDuplicates);
  NumPairs = sortab.nrow();

  // Create the index for all baselines.
  Vector<Int> ant1 = ROScalarColumn<Int>(sortab, "ANTENNA1").getColumn();
  Vector<Int> ant2 = ROScalarColumn<Int>(sortab, "ANTENNA2").getColumn();
  NumAntennae = std::max (max(ant1), max(ant2));
  int index = 0;
  BaselineIndex.resize (NumAntennae*NumAntennae);
  std::fill (BaselineIndex.begin(), BaselineIndex.end(), -1);
  for (uInt i=0; i<ant1.size(); ++i) {
    BaselineIndex[ant1[i] * NumAntennae + ant2[i]] = index++;
  }

  //calculate number of Bands
  // Take care of possible missing time slots.
  if (NumSamples)
  { NumBands  = int(double(NumSamples) / (NumPairs * NumTimeslots) + 0.9);
  }
  else
  { NumBands = spectral_window.nrow();
  }

  ComputeBaselineLengths(MS);

  cout << "Info: " <<NumSamples<<' '<<NumPairs<<' '<<NumAntennae<<' '<<index<<endl;
}

//===============>>>  Ms_Info::~Ms_Info  <<<===============

MsInfo::~MsInfo()
{
}


void MsInfo::update (const MeasurementSet& MS, int timestep)
{
  //Number of channels in the Band
  MSSpectralWindow spectral_window = MS.spectralWindow();
  ROScalarColumn<Int>              NUM_CHAN_col(spectral_window, "NUM_CHAN");
  NumChannels                      = NUM_CHAN_col(0);
  NumTimeslots = (NumTimeslots + timestep - 1) / timestep;
  NumSamples = NumTimeslots * NumPairs * NumBands;
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
      int inx = getBaselineIndex (i, j);
      if (inx >= 0) {
        Vector<Double> p(position(i) - position(j));
        double temp   = sqrt(p(0)*p(0) + p(1)*p(1) + p(2)*p(2));
        BaselineLengths[inx] = temp;
        if (temp > MaxBaselineLength && temp < 3000000) {
          //radius of the Earth in meters? WSRT sometimes has fake telescopes at 3854243
          MaxBaselineLength = temp;   // non-existent antenna's can have position (0,0,0)
        }
      }
    }
  }
}

//===============>>> MsInfo  <<<===============

