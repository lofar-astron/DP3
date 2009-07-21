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

TimeBuffer::TimeBuffer()
{
}

//===============>>>  TimeBuffer::~TimeBuffer  <<<===============

TimeBuffer::~TimeBuffer()
{
}

//===============>>> TimeBuffer::Squash <<<===============
/* Clears the buffer */
void TimeBuffer::Clear(void)
{
  Time.resize(0);
  TimeCentroid.resize(0);
  Interval.resize(0);
  Exposure.resize(0);
  Uvw.resize(0);
}

//===============>>> TimeBuffer::Squash <<<===============
/* Does the time compression, uses Time.size() to determine how. */
void TimeBuffer::Squash(void)
{
  if (Time.size() > 1)
  {
    unsigned int count   = Time.size();
    Double time          = 0.0;
    Double time_centroid = 0.0;
    Double interval      = 0.0;
    Double exposure      = 0.0;
    Vector<double> uvw;  //what a cludge, should be able to do this simpler.
    uvw.resize(3);
    uvw                  = 0.0;
    Vector<double> dummy;
    dummy = uvw;
    dummy = count;
    for (unsigned int i = 0; i < count; i++)
    {
      time          += Time[i];
      time_centroid += TimeCentroid[i];
      interval      += Interval[i];
      exposure      += Exposure[i];
      uvw           += Uvw[i];
    }
    Time[0]         = time / count;
    TimeCentroid[0] = time_centroid / count;
    Interval[0]     = interval;
    Exposure[0]     = exposure;
    Uvw[0]          = uvw / dummy;
  }
}

//===============>>> TimeBuffer::PrintInfo  <<<===============

void TimeBuffer::PrintInfo(void)
{
  std::cout << "Time:           " << Time.size() << std::endl;
  std::cout << "TimeCentroid:   " << TimeCentroid.size() << std::endl;
  std::cout << "Interval:       " << Interval.size() << std::endl;
  std::cout << "Exposure:       " << Exposure.size() << std::endl;
  std::cout << "Uvw:            " << Uvw.size() << std::endl;
}

//===============>>> TimeBuffer  <<<===============
