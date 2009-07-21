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
#ifndef __CS1_PP_TIMEBUFFER_H__
#define __CS1_PP_TIMEBUFFER_H__

#include <casa/Arrays.h>

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
         TimeBuffer();
         ~TimeBuffer();

        std::vector<casa::Double>                 Time;
        std::vector<casa::Double>                 TimeCentroid;
        std::vector<casa::Double>                 Interval;
        std::vector<casa::Double>                 Exposure;
        std::vector< casa::Vector<casa::Double> > Uvw;
        void Squash(void); ///Does the time compression, uses Time.size() to determine how.
        void Clear(void);  ///Clears after a Squash
        void PrintInfo(void);

      private:
    }; // TimeBuffer
  }; //CS1
}; // namespace LOFAR

#endif // __CS1_PP_TIMEBUFFER_H__
