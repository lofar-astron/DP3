//# MSUpdater.cc: DPPP step updating an MS
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

#include <lofar_config.h>
#include <DPPP/MSUpdater.h>
#include <DPPP/MSReader.h>
#include <DPPP/MSWriter.h>
#include <DPPP/DPBuffer.h>
#include <Common/ParameterSet.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/ColumnDesc.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/LinearSearch.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSUpdater::MSUpdater (MSReader* reader, String msName,
                          const ParameterSet& parset,
                          const string& prefix, bool writeHistory)
      : itsReader         (reader),
        itsName           (prefix),
        itsMSName         (msName),
        itsParset         (parset),
        itsWriteData      (false),
        itsWriteWeights   (false),
        itsWriteFlags     (false),
        itsNrDone         (0),
        itsDataColAdded   (false),
        itsWeightColAdded (false),
        itsWriteHistory   (writeHistory)
    {
      itsDataColName   = parset.getString (prefix+"datacolumn",  "");
      itsWeightColName = parset.getString (prefix+"weightcolumn","");
      itsNrTimesFlush  = parset.getUint (prefix+"flush", 0);
    }

    MSUpdater::~MSUpdater()
    {}

    bool MSUpdater::addColumn (const string& colName, const casa::DataType
                               dataType, const ColumnDesc& cd)
    {
      if (itsMS.tableDesc().isColumn(colName)) {
        const ColumnDesc& cd = itsMS.tableDesc().columnDesc(colName);
        ASSERTSTR (cd.dataType() == dataType  &&  cd.isArray(),
                   "Column " << itsDataColName
                   << " already exists, but is not of the right type");
        return false;
      }

      TableDesc td;
      td.addColumn (cd, colName);

      // Use the same data manager as the DATA column.
      // Get the data manager info and find the DATA column in it.
      Record dminfo = itsMS.dataManagerInfo();
      Record colinfo;
      for (uInt i=0; i<dminfo.nfields(); ++i) {
        const Record& subrec = dminfo.subRecord(i);
        if (linearSearch1 (Vector<String>(subrec.asArrayString("COLUMNS")),
                           "DATA") >= 0) {
          colinfo = subrec;
          break;
        }
      }
      ASSERT(colinfo.nfields()>0);
      colinfo.define ("NAME", colName + "_dm");
      itsMS.addColumn (td, colinfo);
      return true;
    }

    bool MSUpdater::process (const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      if (itsWriteFlags) {
        putFlags (buf.getRowNrs(), buf.getFlags());
      }
      if (itsWriteData) {
        putData (buf.getRowNrs(), buf.getData());
      }
      if (itsWriteWeights) {
        if (!buf.getWeights().empty()) {
          // Use weights from buffer
          putWeights (buf.getRowNrs(), buf.getWeights());
        } else {
          itsBuffer.referenceFilled (buf);
          putWeights (buf.getRowNrs(),
                      itsReader->fetchWeights(buf, itsBuffer, itsTimer));
        }
      }
      itsNrDone++;
      if (itsNrTimesFlush > 0  &&  itsNrDone%itsNrTimesFlush == 0) {
        itsMS.flush();
      }
      getNextStep()->process(buf);
      return true;
    }

    void MSUpdater::finish()
    {}

    void MSUpdater::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      itsWriteFlags  = getInfo().writeFlags();

      String origDataColName = getInfo().getDataColName();
      if (itsDataColName.empty()) {
        itsDataColName = origDataColName;
      } else if (itsDataColName != origDataColName) {
        info().setNeedVisData();
        info().setWriteData();
      }
      itsWriteData = getInfo().writeData();

      String origWeightColName = getInfo().getWeightColName();
      if (itsWeightColName.empty()) {
        if (origWeightColName == "WEIGHT") {
          itsWeightColName = "WEIGHT_SPECTRUM";
        } else {
          itsWeightColName = origWeightColName;
        }
      }
      ASSERT(itsWeightColName != "WEIGHT");
      if (itsWeightColName != origWeightColName) {
        info().setWriteWeights();
      }
      itsWriteWeights = getInfo().writeWeights();

      if (getInfo().metaChanged()) {
        THROW(Exception, "Update step " + itsName + 
              " is not possible because meta data changes"
              " (by averaging, adding/removing stations, etc.)");
      }

      if (itsWriteData || itsWriteFlags || itsWriteWeights) {
        NSTimer::StartStop sstime(itsTimer);
        itsMS = MeasurementSet (itsMSName, TableLock::AutoNoReadLocking,
                                Table::Update);
        // Add the data + weight column if needed and if it does not exist yet.
        if (itsWriteData) {
          // use same layout as DATA column
          ColumnDesc cd = itsMS.tableDesc().columnDesc("DATA");
          itsDataColAdded = addColumn(itsDataColName, TpComplex, cd);
        }
        if (itsWriteWeights) {
          IPosition dataShape =
            itsMS.tableDesc().columnDesc("DATA").shape();
          ArrayColumnDesc<float> cd("WEIGHT_SPECTRUM", "weight per corr/chan",
                                    dataShape, ColumnDesc::FixedShape);
          itsWeightColAdded = addColumn(itsWeightColName, TpFloat, cd);
        }
      }
      // Subsequent steps have to set again if writes need to be done.
      info().clearWrites();
      info().clearMetaChanged();
      // Tell the reader if visibility data needs to be read.
      itsReader->setReadVisData (info().needVisData());
    }
      
    void MSUpdater::addToMS (const string&)
    {
      getPrevStep()->addToMS (itsMSName);
      if (itsWriteHistory) {
        MSWriter::writeHistory (itsMS, itsParset);
      }
    }

    void MSUpdater::show (std::ostream& os) const
    {
      os << "MSUpdater " << itsName << std::endl;
      os << "  MS:             " << itsMSName << std::endl;
      os << "  datacolumn:     " << itsDataColName;
      if (itsDataColAdded) {
        os << "  (has been added to the MS)";
      }
      os << std::endl;
      os << "  weightcolumn    " << itsWeightColName;
      if (itsWeightColAdded) {
        os << "  (has been added to the MS)";
      }
      os << std::endl;
      if (itsWriteData || itsWriteFlags || itsWriteWeights) {
        os << "  writing:       ";
        if (itsWriteData)    os << " data";
        if (itsWriteFlags)   os << " flags";
        if (itsWriteWeights) os << " weights";
        os << std::endl;
      }
      os << std::endl;
      os << "  flush:          " << itsNrTimesFlush << std::endl;
    }

    void MSUpdater::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " MSUpdater " << itsName << endl;
    }

    void MSUpdater::putFlags (const RefRows& rowNrs,
                              const Cube<bool>& flags)
    {
      // Only put if rownrs are filled, thus if data were not inserted.
      if (! rowNrs.rowVector().empty()) {
        Slicer colSlicer(IPosition(2, 0, info().startchan()),
                         IPosition(2, info().ncorr(), info().nchan()) );
        ArrayColumn<bool> flagCol(itsMS, "FLAG");
        ScalarColumn<bool> flagRowCol(itsMS, "FLAG_ROW");
        // Loop over all rows of this subset.
        // (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<bool> flagIter (flags, 2);
        for (uint i=0; i<rows.size(); ++i) {
          flagCol.putSlice (rows[i], colSlicer, flagIter.array());
          // If a new flag in a row is clear, the ROW_FLAG should not be set.
          // If all new flags are set, we leave it because we might have a
          // subset of the channels, so other flags might still be clear.
          if (anyEQ (flagIter.array(), False)) {
            flagRowCol.put (rows[i], False);
          }
          flagIter.next();
        }
      }
    }

    void MSUpdater::putWeights (const RefRows& rowNrs,
                                const Cube<float>& weights)
    {
      // Only put if rownrs are filled, thus if data were not inserted.
      if (! rowNrs.rowVector().empty()) {
        Slicer colSlicer(IPosition(2, 0, info().startchan()),
                         IPosition(2, info().ncorr(), info().nchan()) );
        ArrayColumn<float> weightCol(itsMS, itsWeightColName);
        // Loop over all rows of this subset.
        // (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<float> weightIter (weights, 2);
        for (uint i=0; i<rows.size(); ++i) {
          weightCol.putSlice (rows[i], colSlicer, weightIter.array());
          weightIter.next();
        }
      }
    }


    void MSUpdater::putData (const RefRows& rowNrs,
                             const Cube<Complex>& data)
    {
      // Only put if rownrs are filled, thus if data were not inserted.
      if (! rowNrs.rowVector().empty()) {
        Slicer colSlicer(IPosition(2, 0, info().startchan()),
                         IPosition(2, info().ncorr(), info().nchan()) );
        ArrayColumn<Complex> dataCol(itsMS, itsDataColName);
        // Loop over all rows of this subset.
        // (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<Complex> dataIter (data, 2);
        for (uint i=0; i<rows.size(); ++i) {
          dataCol.putSlice (rows[i], colSlicer, dataIter.array());
          dataIter.next();
        }
      }
    }
  } //# end namespace
}
