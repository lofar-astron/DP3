//# MSWriter.h: DPPP step writing to an MS
//# Copyright (C) 2010
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
//# @author Ger van Diepen

#ifndef DPPP_MSWRITER_H
#define DPPP_MSWRITER_H

// @file
// @brief DPPP step writing to an MS

#include <DPPP/DPStep.h>
#include <DPPP/MSReader.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ColumnDesc.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/TiledColumnStMan.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {
    class AverageInfo;

    // @ingroup DPPP

    class MSWriter: public DPStep
    {
    public:
      explicit MSWriter (MSReader* reader, const std::string& outName,
                         const AverageInfo&,
                         const ParameterSet&, const string& prefix);

      virtual ~MSWriter();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Write the parset info into the HISTORY table of the MS.
      static void writeHistory (casa::Table& ms,
                                const ParameterSet& parset);

    private:
      // Create an array column description and add to table with given
      // stoage manager (if given).
      void makeArrayColumn (casa::ColumnDesc desc, const casa::IPosition& shape,
                            casa::TiledColumnStMan* tsm, casa::Table& table);

      // Create the MS by cloning all subtables from the input MS.
      // All output columns in the main table are using normal storage managers.
      // The SPECTRAL_WINDOW table is adapated as needed.
      void createMS (const std::string& outName, const AverageInfo& avgInfo,
                     uint tileSize, uint tileNChan);

      // Update the SPECTRAL_WINDOW table for averaged channels.
      void updateSpw (const string& outName, const AverageInfo& avgInfo);

      // Write the data, flags, etc.
      void writeData (casa::Table& out, const DPBuffer& buf);

      // Write the full resolution flags (flags before any averaging).
      void writePreAvgFlags (casa::Table& out, const DPBuffer& buf);

      // Write the time info (TIME, TIME_CENTROID, INTERVAL, EXPOSURE).
      void writeTimeInfo (casa::Table& out, double time,
                          const casa::Matrix<double>& uvws);

      // Copy meta data columns for a time slot (ANTENNA1, etc.)
      // It also copies all time info if possible.
      void copyMeta (const casa::Table& in, casa::Table& out,
                     bool copyTimeInfo);

      // Copy the contents of a scalar column.
      template<typename T> void copySca (const casa::Table& in,
                                         casa::Table& out,
                                         const casa::String& columnName)
      {
        casa::ROScalarColumn<T> inCol(in, columnName);
        casa::ScalarColumn<T>  outCol(out, columnName);
        outCol.putColumn (inCol.getColumn());
      }

      // Copy the contents of an array column.
      template<typename T> void copyArr (const casa::Table& in,
                                         casa::Table& out,
                                         const casa::String& columnName)
      {
        casa::ROArrayColumn<T> inCol(in, columnName);
        casa::ArrayColumn<T>  outCol(out, columnName);
        outCol.putColumn (inCol.getColumn());
      }

      //# Data items.
      MSReader*       itsReader;
      casa::Table     itsMS;
      casa::String    itsDataColName;
      double          itsInterval;
      bool            itsCopyTimeInfo;
      bool            itsCopyCorrData;
      bool            itsCopyModelData;
      bool            itsWritePreAvgFlags;
      uint            itsNrCorr;
      uint            itsNrChan;
      uint            itsNrBl;
      uint            itsNrTimes;
      uint            itsOrigNrChan;    //# original nr of chan before averaging
      uint            itsNTimeAvg;      //# nr of times averaged to 1
      bool            itsCountFlags;
      FlagCounter     itsFlagCounter;
    };

  } //# end namespace
}

#endif
