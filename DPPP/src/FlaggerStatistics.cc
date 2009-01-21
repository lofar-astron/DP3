/***************************************************************************
 *   Copyright (C) 2007-8 by ASTRON, Adriaan Renting                       *
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
#include <DPPP/FlaggerStatistics.h>

using namespace LOFAR::CS1;
using namespace casa;

//===============>>>  int2string  <<<===============
//ANSI C++ doesn't seem to have a decent function for this, or I'm not aware of it. Need to rename it to IntToStr(), to avoid confusion
std::string int2string(int value, int base=10)
{ // found on http://www.jb.man.ac.uk/~slowe/cpp/int2string.html
  //maybe copyright Robert Jan Schaper, Ray-Yuan Sheu, Rodrigo de Salvo Braz, Wes Garland and John Maloney
  enum { kMaxDigits = 35 };
  std::string buf;
  buf.reserve( kMaxDigits ); // Pre-allocate enough space.
        // check that the base if valid
  if (base < 2 || base > 16) return buf;
  int quotient = value;
        // Translating number to string with base:
  do {
    buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
    quotient /= base;
  } while ( quotient );
        // Append the negative sign for base 10
  if ( value < 0 && base == 10) buf += '-';
  std::reverse( buf.begin(), buf.end() );
  return buf;
}


//===============>>>  FlaggerStatistics::FlaggerStatistics  <<<===============
/* initialize some meta data and get the datastorage the right size. */
FlaggerStatistics::FlaggerStatistics(MsInfo& info)
{
  NumBands     = info.NumBands;
  NumAntennae  = info.NumAntennae;
  AntennaNames = info.AntennaNames;
  Normalizer   = info.NumChannels * info.NumTimeslots * info.NumPolarizations;
  Statistics   = Cube<int>(NumBands, NumAntennae, NumAntennae, 0);
}

//===============>>>  FlaggerStatistics::~FlaggerStatistics  <<<===============

FlaggerStatistics::~FlaggerStatistics()
{
}

//===============>>>  FlaggerStatistics::~operator[]  <<<===============

int& FlaggerStatistics::operator()(int x, int y, int z)
{
  return Statistics(x, y, z);
}


//===============>>> FlaggerStatistics::FlagDataOrBaselines  <<<===============
/*This function outputs the gathered statistics.*/
void FlaggerStatistics::PrintStatistics(ostream& output)
{
  vector<int>  bands(NumBands);
  vector<int>  antennae(NumAntennae);
  unsigned int namelength = 6;
  for(int i = 0; i < NumAntennae; i++)
  {
    if (namelength < AntennaNames[i].size())
    { namelength = AntennaNames[i].size();
    }
  }
  for (int i = 0; i < NumBands; i++)
  {
    output << "Band: " << i+1 << endl;
    output << string(namelength+1,' ');
    for(int j = 0; j < NumAntennae; j++)
    {
      string out = AntennaNames[j];
      out.resize(namelength+1,' ');
      output << out;
    }
    output << endl;
    for(int j = 0; j < NumAntennae; j++)
    {
      string out = AntennaNames[j];
      out.resize(namelength+1,' ');
      output << out;
      for(int k = 0; k < NumAntennae; k++)
      {
        if (k < j) //We print a complete array, but we have inspected only those where k >= j
        {
          int val = 100 * Statistics(i,k,j) / (Normalizer);
          bands[i]    += val;
          antennae[j] += val;
          antennae[k] += val;
          string out = int2string(val) + "%";
          out.resize(namelength+1,' ');
          output << out;
        }
        else
        {
          int val = 100 * Statistics(i,j,k) / (Normalizer);
          bands[i]    += val;
          antennae[j] += val;
          antennae[k] += val;
          string out = int2string(val) + "%";
          out.resize(namelength+1,' ');
          output << out;
        }
      }
      output << endl;
    }
  }
  output << "Bands (flagged %):    ";
  for (int i = 0; i < NumBands; i++)
  {
    string out = string("BND") + int2string(i);
    out.resize(namelength+1,' ');
    output << out;
  }
  output << endl << "                      ";
  for (int i = 0; i < NumBands; i++)
  {
    string out = int2string(bands[i] / (NumAntennae*NumAntennae)) + "%";
    out.resize(namelength+1,' ');
    output << out;
  }
  output << endl << "Antennae (flagged %): " ;
  for(int j = 0; j < NumAntennae; j++)
  {
    string out = AntennaNames[j];
    out.resize(namelength+1,' ');
    output << out;
  }
  output << endl << "                       ";
  for (int i = 0; i < NumAntennae; i++)
  {
    string out = int2string(antennae[i] / (NumBands*NumAntennae*2)) + "%";
    out.resize(namelength+1,' ');
    output << out;
  }
  output << endl;
}
//===============>>> FlaggerStatistics  <<<===============
