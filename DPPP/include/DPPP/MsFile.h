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
#ifndef __CS1_PP_MS_FILE_H__
#define __CS1_PP_MS_FILE_H__

#include <ms/MeasurementSets.h>
#include <tables/Tables.h>
#include <tables/Tables/TableIter.h>
#include <DPPP/MsInfo.h>
#include <DPPP/DataBuffer.h>
#include <DPPP/TimeBuffer.h>

/// @file
/// @brief Class to hold code for access to MeasurementSet in IDPPP
/// @author Adriaan Renting (renting AT astron nl)

namespace LOFAR
{
  namespace CS1
  {
    ///Foreward declaration
    class MsInfo;
    class RunDetails;

    class MsFile
    {
    public:
      MsFile(const std::string& msin, const std::string& msout);
      ~MsFile();

      /// Iterator to process all data from one integration time at the same time
      casa::TableIterator TimeIterator();
      /// creates a new measurement set
      void Init(MsInfo& Info, RunDetails& Details, int Squashing);
      void PrintInfo(void); ///< prints some numbers for debug purposes
      /// processes the tables from TimeIterator to fill the next timeslot to be processed in DataBuffer
      void UpdateTimeslotData(casa::TableIterator& Data_iter,
                              MsInfo& Info,
                              DataBuffer& Buffer,
                              TimeBuffer& TimeData);
      /// Writes the data in DataBuffer->Position+1 to the file
      void WriteData(casa::TableIterator& Data_iter,
                     MsInfo& Info,
                     DataBuffer& Buffer,
                     TimeBuffer& TimeData);


    protected:
    private:
      /// Function for reshaping a table column
      void TableResize(casa::TableDesc tdesc,
                       casa::IPosition ipos,
                       std::string name,
                       casa::Table& table);
      casa::Block<casa::String> SELECTblock;
      std::string InName;
      std::string OutName;
      casa::MeasurementSet* InMS;
      casa::MeasurementSet* OutMS;
    }; // class MsFile
  }; // CS1
}; // namespace LOFAR

#endif // __CS1_PP_MS_FILE_H__
