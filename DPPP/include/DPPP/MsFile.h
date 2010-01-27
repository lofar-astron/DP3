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

#ifndef __CS1_PP_MS_FILE_H__
#define __CS1_PP_MS_FILE_H__

#include <Common/Exception.h>
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
    EXCEPTION_CLASS(PipelineException, LOFAR::Exception);

    //# Forward declaration
    class MsInfo;
    class RunDetails;

    class MsFile
    {
    public:
      MsFile(const std::string& msin, const std::string& msout);
      ~MsFile();

      // Get the input MS.
      const casa::MeasurementSet& getInMS() const
        { return *InMS; }
      // Get the ordered input main table.
      const casa::Table& getOrderedInTable() const
        { return itsOrderedTable; }
      // Get the output MS.
      const casa::MeasurementSet& getOutMS() const
        { return *OutMS; }

      /// Iterator to process all data from one integration time at the same time
      casa::TableIterator TimeIterator();
      /// Get nr of rows in input MS.
      unsigned int nrow() const
        { return InMS->nrow(); }
      /// creates a new measurement set, returns if the input has imaging columns
      void Init(MsInfo& Info, RunDetails& Details, int Squashing);
      void PrintInfo(void); ///< prints some numbers for debug purposes
      /// processes the tables from TimeIterator to fill the next timeslot to be processed in DataBuffer
      void UpdateTimeslotData(casa::TableIterator& Data_iter,
                              MsInfo& Info,
                              DataBuffer& Buffer,
                              TimeBuffer& TimeData,
			      bool missingTime,
			      double timeValue);
      /// Writes the data in DataBuffer->Position+1 to the file
      void WriteData(casa::TableIterator& Data_iter,
                     MsInfo& Info,
                     DataBuffer& Buffer,
                     TimeBuffer& TimeData);

      /// Flush the data.
      void flush()
        { OutMS->flush (true); }

    protected:
    private:
      /// Function for adding a table column.
      void TableResize(casa::ColumnDesc desc, const casa::IPosition& ipos,
                       casa::TiledColumnStMan* tsm, casa::Table& table);

      // Check for time gaps.
      void checkGaps(const MsInfo&, bool updateMS) const;

      casa::IPosition DetermineDATAshape(const casa::Table& MS);
      casa::Block<casa::String> SELECTblock;
      std::string InName;
      std::string OutName;
      casa::MeasurementSet* InMS;
      casa::MeasurementSet* OutMS;
      casa::Table           itsOrderedTable;
      casa::TableIterator   itsIterator;
      bool itsHasWeightSpectrum;
      bool itsIsOrdered;
    }; // class MsFile
  }; // CS1
}; // namespace LOFAR

#endif // __CS1_PP_MS_FILE_H__
