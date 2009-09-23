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

#ifndef __CS1_PP_DATABUFFER_H__
#define __CS1_PP_DATABUFFER_H__

#include <casa/Arrays.h>

#include <DPPP/MsInfo.h>

/// @file
/// @brief Class to hold code for DataBuffer in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

/// The class maintains a buffer withe NumSlots Cubes of size:
/// [myInfo->NumPolarizations, myInfo->NumChannels, WindowSize]
/// ModelData and CorrectedData are only used if Colums==true
namespace LOFAR
{
  namespace CS1
  {
    class DataBuffer
    {
      public:
         DataBuffer(MsInfo* info, int TimeWindow, bool Columns);
         ~DataBuffer();
        int Position; ///< Position just updated in the buffer, -1 when it's uninitialised
        int NumSlots; ///< Number fo baselines X number of Spectral Windows
        int WindowSize;

        std::vector< bool >                       PolarizationsToCheck;
        std::vector< casa::Cube<casa::Complex> >  Data;
        std::vector< casa::Cube<casa::Complex> >  ModelData;
        std::vector< casa::Cube<casa::Complex> >  CorrectedData;
        std::vector< casa::Cube<casa::Bool> >     Flags;
        std::vector< casa::Cube<casa::Float> >    Weights;
        std::vector< casa::Matrix<casa::Double> > Uvw;
        void DeterminePolarizationsToCheck(bool UseOnlyXpolarizations); ///< Not used right now
        void PrintInfo(void);

      private:
        MsInfo* myInfo;
        void    Init(bool Columns);
    }; // DataBuffer
  }; //CS1
}; // namespace LOFAR

#endif // __CS1_PP_DATABUFFER_H__
