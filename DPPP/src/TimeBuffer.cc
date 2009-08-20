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

#include <DPPP/TimeBuffer.h>
#include <casa/Arrays/ArrayMath.h>

using namespace LOFAR::CS1;
using namespace casa;
using casa::Vector;
using casa::Double;


//===============>>>  TimeBuffer::TimeBuffer  <<<===============

TimeBuffer::TimeBuffer(const int numslots):
  NumSlots(numslots)
{
  BufTime.resize(numslots);
  BufTimeCentroid.resize(numslots);
  BufInterval.resize(numslots);
  BufExposure.resize(numslots);
  BufUvw.resize(numslots);
  Time.resize(numslots);
  TimeCentroid.resize(numslots);
  Interval.resize(numslots);
  Exposure.resize(numslots);
  Uvw.resize(numslots);
}

//===============>>>  TimeBuffer::~TimeBuffer  <<<===============

TimeBuffer::~TimeBuffer()
{
}

//===============>>> TimeBuffer::Squash <<<===============
/* Clears the buffer */
void TimeBuffer::Clear(void)
{
  for (int i = 0; i < NumSlots; i++)
  {
    Time[i].resize(0);
    TimeCentroid[i].resize(0);
    Interval[i].resize(0);
    Exposure[i].resize(0);
    Uvw[i].resize(0);
  }
}

//===============>>> TimeBuffer::Squash <<<===============
/* Does the time compression, uses Time.size() to determine how. */
void TimeBuffer::Squash(void)
{
  for (int i = 0; i < NumSlots; i++)
  {
    unsigned int count   = Time[i].size();
    if (count > 1)
    {
      Double time          = 0.0;
      Double time_centroid = 0.0;
      Double interval      = 0.0;
      Double exposure      = 0.0;
      Vector<double> uvw(3, 0.0);
      Vector<double> dummy(3, count);
      for (unsigned int j = 0; j < count; j++)
      {
        time          += Time[i][j];
        time_centroid += TimeCentroid[i][j];
        interval      += Interval[i][j];
        exposure      += Exposure[i][j];
        uvw           += Uvw[i][j];
      }
      Time[i][0]         = time / count;
      TimeCentroid[i][0] = time_centroid / count;
      Interval[i][0]     = interval;
      Exposure[i][0]     = exposure;
      Uvw[i][0]          = uvw / dummy;
    }
  }
}

//===============>>> TimeBuffer::PrintInfo  <<<===============

void TimeBuffer::ShiftBuffer(void)
{
  for (int i = 0; i < NumSlots; i++)
  {
    Time[i].push_back(BufTime[i].back());
    BufTime[i].pop_back();
    TimeCentroid[i].push_back(BufTimeCentroid[i].back());
    BufTimeCentroid[i].pop_back();
    Interval[i].push_back(BufInterval[i].back());
    BufInterval[i].pop_back();
    Exposure[i].push_back(BufExposure[i].back());
    BufExposure[i].pop_back();
    Uvw[i].push_back(BufUvw[i].back());
    BufUvw[i].pop_back();
  }
}

//===============>>> TimeBuffer::PrintInfo  <<<===============

void TimeBuffer::PrintInfo(void)
{
  for (int i = 0; i < NumSlots; i++)
  {
    std::cout << "BufTime:        " << BufTime[i].size() << std::endl;
    std::cout << "BufTimeCentroid:" << BufTimeCentroid[i].size() << std::endl;
    std::cout << "BufInterval:    " << BufInterval[i].size() << std::endl;
    std::cout << "BufExposure:    " << BufExposure[i].size() << std::endl;
    std::cout << "BufUvw:         " << BufUvw[i].size() << std::endl;
    std::cout << "Time:           " << Time[i].size() << std::endl;
    std::cout << "TimeCentroid:   " << TimeCentroid[i].size() << std::endl;
    std::cout << "Interval:       " << Interval[i].size() << std::endl;
    std::cout << "Exposure:       " << Exposure[i].size() << std::endl;
    std::cout << "Uvw:            " << Uvw[i].size() << std::endl;
  }
}

//===============>>> TimeBuffer  <<<===============
