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
#ifndef __CS1_PP_MS_INFO_H__
#define __CS1_PP_MS_INFO_H__

#include <ms/MeasurementSets.h>
#include <tables/Tables.h>
#include <utility>
#include <vector>
#include <string>
#include <map>

/// @file
/// @brief Class to hold code for retrieving information from MeasurementSet in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    typedef std::pair<int, int> baseline_t;

    class MsInfo
    {
    public:
      MsInfo(const std::string& msname);
      ~MsInfo();

      int                       NumSamples;
      int                       NumAntennae;
      int                       NumFields;
      int                       NumBands;
      int                       NumChannels;
      int                       NumPolarizations;
      int                       NumPairs;
      int                       NumTimeslots;
      double                    NoiseLevel;
      std::vector<casa::String> AntennaNames;
      casa::Vector<casa::Int>   Polarizations;
      double                    MaxBaselineLength;
      std::vector<baseline_t>   PairsIndex;
      std::map<baseline_t, int> BaselineIndex;
      std::vector<double>       BaselineLengths;

      void                      Update(void); ///< reread the info from the MeasurementSet
      void                      PrintInfo(void); ///< prints info to cout, for debugging

    protected:
    private:
      std::string MsName;
      void        ComputeBaselineLengths(casa::MeasurementSet& MS); ///< for baseline dependent flagging or filtering
    }; // class MsInfo
  }; // CS1
}; // namespace LOFAR

#endif // __CS1_PP_MS_INFO_H__
