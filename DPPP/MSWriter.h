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

#include "DPStep.h"
#include "MSReader.h"
#include "StManParsetKeys.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/DataMan/TiledColumnStMan.h>

namespace DP3 {
  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep creating a new MeasurementSet and writing
    // all data in it.
    // Most meta information (subtables and meta columns in main table) is
    // copied from the input MeasurementSet given by the MSReader object.
    // <br>
    // In principle the new MS uses the same storage managers as used in the
    // input MS, but in case of an MS stored with LofarStMan it will use the
    // optimal storage managers (ISM for slowly varying meta data, TSM for
    // bulk data, SSM for others).
    //
    // The SPECTRAL_WINDOW table will be changed to reflect the channels
    // being used or averaged.
    // The OBSERVATION table will be updated for the correct start and end time.
    // The HISTORY table gets an entry containing the parset values and the
    // DPPP version.

    class MSWriter: public DPStep
    {
    public:
      explicit MSWriter (MSReader* reader, const std::string& outName,
                         const ParameterSet&, const string& prefix);

      virtual ~MSWriter();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Add some data to the MeasurementSet written/updated.
      // Calls addToMS from the previous step, with the current output msname.
      virtual void addToMS (const string&);

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

      // Write the parset info into the HISTORY table of the MS.
      static void writeHistory (casacore::Table& ms,
                                const ParameterSet& parset);

      static void updateBeam (const std::string& outName, const std::string& outColName, const DPInfo& info);

    private:
      // Create an array column description and add to table with given
      // stoage manager (if given).
      void makeArrayColumn (casacore::ColumnDesc desc, const casacore::IPosition& shape,
                            casacore::DataManager* dm, casacore::Table& table, bool makeDirectColumn = false);

      // Create the MS by cloning all subtables from the input MS.
      // All output columns in the main table are using normal storage managers.
      // The SPECTRAL_WINDOW table is adapted as needed.
      void createMS (const std::string& outName, const DPInfo& info,
                     uint tileSize, uint tileNChan);

      // Update the SPECTRAL_WINDOW table for averaged channels.
      void updateSpw (const string& outName, const DPInfo& info);

      // Update the OBSERVATION table with the correct start and end time.
      void updateObs (const string& outName);

      // Update the FIELD table with the new phase center.
      void updatePhaseCentre (const string& outName, const DPInfo& info);
      
      // Write the data, flags, etc.
      void writeData (casacore::Table& out, const DPBuffer& buf);

      // Write the full resolution flags (flags before any averaging).
      void writeFullResFlags (casacore::Table& out, const DPBuffer& buf);

      // Write all meta data columns for a time slot (ANTENNA1, etc.)
      void writeMeta (casacore::Table& out, const DPBuffer& buf);

      // Copy meta data columns for a time slot (ANTENNA1, etc.)
      // It also copies all time info if possible.
      void copyMeta (const casacore::Table& in, casacore::Table& out,
                     bool copyTimeInfo);

      // Copy the contents of a scalar column.
      template<typename T> void fillSca (const T& value,
                                         casacore::Table& out,
                                         const casacore::String& columnName)
      {
        casacore::ScalarColumn<T>  outCol(out, columnName);
        outCol.fillColumn (value);
      }

      // Copy the contents of an array column.
      template<typename T> void fillArr (const casacore::Array<T>& value,
                                         casacore::Table& out,
                                         const casacore::String& columnName)
      {
        casacore::ArrayColumn<T> outCol(out, columnName);
        outCol.fillColumn (value);
      }

      // Copy the contents of a scalar column.
      template<typename T> void copySca (const casacore::Table& in,
                                         casacore::Table& out,
                                         const casacore::String& columnName)
      {
        casacore::ROScalarColumn<T> inCol(in, columnName);
        casacore::ScalarColumn<T>  outCol(out, columnName);
        outCol.putColumn (inCol.getColumn());
      }

      // Copy the contents of an array column.
      template<typename T> void copyArr (const casacore::Table& in,
                                         casacore::Table& out,
                                         const casacore::String& columnName)
      {
        casacore::ROArrayColumn<T> inCol(in, columnName);
        casacore::ArrayColumn<T>  outCol(out, columnName);
        outCol.putColumn (inCol.getColumn());
      }

      //# Data items.
      MSReader*       itsReader;
      string          itsName;
      string          itsOutName;
      DPBuffer        itsBuffer;
      casacore::Table     itsMS;
      ParameterSet    itsParset; //# parset for writing history
      casacore::String    itsDataColName;
      casacore::String    itsWeightColName;
      double          itsInterval;
      bool            itsOverwrite;   //# Overwrite an existing output MS?
      bool            itsCopyCorrData;
      bool            itsCopyModelData;
      bool            itsWriteFullResFlags;
      uint            itsTileSize;
      uint            itsTileNChan;
      uint            itsNrCorr;
      uint            itsNrChan;
      uint            itsNrBl;
      uint            itsNrTimes;
      uint            itsNChanAvg;    //# nr of channels in input averaged to 1
      uint            itsNTimeAvg;    //# nr of times in input averaged to 1
      uint            itsNrTimesFlush;//# flush every N time slots (0=no flush)
      uint            itsNrDone;      //# nr of time slots written
      std::string     itsVdsDir;      //# directory where to put VDS file
      std::string     itsClusterDesc; //# name of clusterdesc file
      NSTimer         itsTimer;
      StManParsetKeys itsStManKeys;
    };

  } //# end namespace
}

#endif
