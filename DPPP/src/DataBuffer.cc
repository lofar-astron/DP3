/***************************************************************************
 *   Copyright (C) 2008 by ASTRON, Adriaan Renting                         *
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
#include <iostream>

#include <DPPP/DataBuffer.h>
#include <DPPP/MsInfo.h>

using namespace LOFAR::CS1;
using namespace casa;
using casa::Cube;
using casa::Matrix;
using casa::Complex;
using casa::Int;
using casa::Double;
using casa::Bool;
using casa::Table;
using std::vector;
using std::map;


enum CorrelationTypes {None=0,I=1,Q=2,U=3,V=4,RR=5,RL=6,LR=7,LL=8,XX=9,XY=10,YX=11,YY=12}; //found somewhere in AIPS++, don't remember where

//===============>>>  DataBuffer::DataBuffer  <<<===============

DataBuffer::DataBuffer(MsInfo* info, int TimeWindow, bool Columns):
  Position(-1),
  NumSlots(0),
  WindowSize(TimeWindow),
  myInfo(info)
{
  NumSlots   = myInfo->NumPairs * myInfo->NumBands;
  Init(Columns);
}

//===============>>>  DataBuffer::~DataBuffer  <<<===============

DataBuffer::~DataBuffer()
{
}

//===============>>> DataBufferr::DetermineCorrelationsToCheck <<<===============
/* create a list of polarizations we want to check, maybe we only want to to XY, YX */
void DataBuffer::DeterminePolarizationsToCheck(bool UseOnlyXpolarizations)
{
  if (UseOnlyXpolarizations)
  {
    bool noCorrError = true;
    for (int i = 0; i < myInfo->NumPolarizations; i++)
    {
      switch(myInfo->Polarizations[i])
      {
        case None:
        case I:
        case RR:
        case LL:
        case XX:
        case YY:
          if (UseOnlyXpolarizations)
            PolarizationsToCheck[i] = false;
          break;
        case Q:
        case U:
        case V:
        case RL:
        case LR:
        case XY:
        case YX:
          noCorrError = false;
          if (UseOnlyXpolarizations)
            PolarizationsToCheck[i] = true;
          break;
      }
    }
    if (noCorrError)
    {
      cout << "There are no crosspolarizations to flag!";
      exit(1);
    }
  }
  else
  {
    for (int i = 0; i < myInfo->NumPolarizations; i++)
    {  PolarizationsToCheck[i] = true;
    }
  }
}

//===============>>> DataBuffer::Init  <<<===============

void DataBuffer::Init(bool Columns)
{
  PolarizationsToCheck.resize(myInfo->NumPolarizations);
//  DeterminePolarizationsToCheck(UseOnlyXpolarizations);

  Data.resize(NumSlots);
  Flags.resize(NumSlots);
  Weights.resize(NumSlots);
  Uvw.resize(NumSlots);
  if (Columns)
  {
    ModelData.resize(NumSlots);
    CorrectedData.resize(NumSlots);
  }
  for (int i = 0; i < NumSlots; i++)
  {
    Data[i].resize(myInfo->NumPolarizations, myInfo->NumChannels, WindowSize);
    Flags[i].resize(myInfo->NumPolarizations, myInfo->NumChannels, WindowSize);
    Weights[i].resize(myInfo->NumPolarizations, myInfo->NumChannels, WindowSize);
    Uvw[i].resize(3, WindowSize);
    if (Columns)
    {
      ModelData[i].resize(myInfo->NumPolarizations, myInfo->NumChannels, WindowSize);
      CorrectedData[i].resize(myInfo->NumPolarizations, myInfo->NumChannels, WindowSize);
    }
  }
}

//===============>>> DataBuffer::PrintInfo  <<<===============

void DataBuffer::PrintInfo(void)
{
  std::cout << "Position        " << Position << std::endl;
  std::cout << "Baselines:      " << Data.size() << std::endl;
  std::cout << "Data:           " << Data[0].shape() << std::endl;
  std::cout << "Flags:          " << Flags[0].shape() << std::endl;
  std::cout << "Weights:        " << Weights[0].shape() << std::endl;
  if (ModelData.size())
  {
    std::cout << "ModelData:      " << ModelData[0].shape() << std::endl;
    std::cout << "CorrectedData:  " << CorrectedData[0].shape() << std::endl;
  }
}

//===============>>> DataBuffer  <<<===============
