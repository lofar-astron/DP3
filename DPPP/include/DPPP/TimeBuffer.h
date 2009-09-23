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

#ifndef __CS1_PP_TIMEBUFFER_H__
#define __CS1_PP_TIMEBUFFER_H__

#include <casa/Arrays.h>
#include <vector>
#include <list>

/// @file
/// @brief Class to hold code for TimeBuffer in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

/// The class maintains a buffer with data that needs special treatment
/// when compressing in time.
namespace LOFAR
{
  namespace CS1
  {
    class TimeBuffer
    {
      public:
         TimeBuffer(const int numslots);
         ~TimeBuffer();

        /// We need a buffer for each variable before we can start processing them. This is a fifo buffer and
        /// thus implemented as a linked list.
        std::vector< std::list<casa::Double> >                   BufTime;
        std::vector< std::list<casa::Double> >                   BufTimeCentroid;
        std::vector< std::list<casa::Double> >                   BufInterval;
        std::vector< std::list<casa::Double> >                   BufExposure;
        std::vector< std::list< casa::Vector<casa::Double> > >   BufUvw;
        std::vector< std::vector<casa::Double> >                 Time;
        std::vector< std::vector<casa::Double> >                 TimeCentroid;
        std::vector< std::vector<casa::Double> >                 Interval;
        std::vector< std::vector<casa::Double> >                 Exposure;
        std::vector< std::vector< casa::Vector<casa::Double> > > Uvw;
        void Squash(void); ///Does the time compression, uses Time.size() to determine how.
        void Clear(void);  ///Clears after a Squash
        void ShiftBuffer(void);
        void PrintInfo(void);

      private:
        int NumSlots;
    }; // TimeBuffer
  }; //CS1
}; // namespace LOFAR

#endif // __CS1_PP_TIMEBUFFER_H__
