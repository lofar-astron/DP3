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

#ifndef __CS1_PP_PIPELINE_H__
#define __CS1_PP_PIPELINE_H__

#include <vector>

/// @file
/// @brief Class to hold code for Pipeline in IDPPP, which determines the order in which the
/// steps in IDPPP are run
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    ///Foreward declarations
    class MsInfo;
    class MsFile;
    class RunDetails;
    class BandpassCorrector;
    class Flagger;
    class DataSquasher;
    class DataBuffer;
    class TimeBuffer;
    class FlaggerStatistics;

    class Pipeline
    {
    public:
      /// Gets initialised after reading the ParameterSet and initalising the input and output
      /// MeasuremenSets and choosing a BandpassCorrector, Flagger and DataSquasher and initialising those
      Pipeline(MsInfo* info, MsFile* msfile, RunDetails* details,
               BandpassCorrector* bandpass, Flagger* flagger, DataSquasher* squasher);
      ~Pipeline();
      ///run the pipeline until there are no more timeslots in the input MeasurementSet.
      void Run(MsInfo* SquashedInfo, bool Columns);

    protected:
    private:
      MsInfo*             myInfo;
      MsFile*             myFile;
      RunDetails*         myDetails;
      BandpassCorrector*  myBandpass;
      Flagger*            myFlagger;
      DataSquasher*       mySquasher;
      DataBuffer*         BandpassData; ///< initial data read from MS
      DataBuffer*         FlaggerData; ///< data after bandpass correction to be flagged
      DataBuffer*         SquasherData; ///< output of the squasher
      FlaggerStatistics*  myStatistics; ///< Stores the statistics of the flaggers
      TimeBuffer*         myTime; ///< remmber what timeslots we are processing
      void MirrorBuffer(DataBuffer& buffer, MsInfo& info, int step); ///< for handling the start and stop edges of the data

    }; // class Pipeline
  }; // CS1
}; // namespace LOFAR

#endif // __CS1_PP_PIPELINE_H__
