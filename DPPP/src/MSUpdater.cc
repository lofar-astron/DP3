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
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/ColumnDesc.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/LinearSearch.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSUpdater::MSUpdater (MSReader* reader, const ParameterSet& parset,
                          const string& prefix, int needWrite)
      : itsReader       (reader),
        itsWriteData    ((needWrite & DPInfo::NeedWriteData) != 0),
        itsWriteWeight  ((needWrite & DPInfo::NeedWriteWeight) != 0),
        itsNrDone       (0),
        itsDataColAdded (false),
        itsWeightColAdded (false)
    {
      itsDataColName  = parset.getString (prefix+"datacolumn",
                                          itsReader->dataColumnName());
      itsWeightColName  = parset.getString (prefix+"weightcolumn",
        itsReader->weightColumnName()=="WEIGHT" ? "WEIGHT_SPECTRUM":
        itsReader->weightColumnName());
      ASSERT(itsWeightColName!="WEIGHT");
      if (itsWeightColName != itsReader->weightColumnName()) {
        itsWriteWeight = true;
      }
      itsNrTimesFlush = parset.getUint (prefix+"flush", 0);
      NSTimer::StartStop sstime(itsTimer);
      // Reopen the MS for read/write.
      itsReader->table().reopenRW();
      // Add the data + weight column if needed and if it does not exist yet.
      if (itsWriteData) {
        // use same layout as DATA column
        ColumnDesc cd = itsReader->table().tableDesc().columnDesc("DATA");
        itsDataColAdded = addColumn(itsDataColName, TpComplex, cd);
      }
      if (itsWriteWeight) {
        IPosition dataShape =
            itsReader->table().tableDesc().columnDesc("DATA").shape();
        ArrayColumnDesc<float> cd("WEIGHT_SPECTRUM", "weight per corr/chan",
            dataShape, ColumnDesc::FixedShape);
        itsWeightColAdded = addColumn(itsWeightColName, TpFloat, cd);
      }
      MSWriter::writeHistory (reader->table(), parset);
    }

    MSUpdater::~MSUpdater()
    {}

    bool MSUpdater::isNewDataColumn (MSReader* reader,
                                     const ParameterSet& parset,
                                     const string& prefix)
    {
      // Only test if the output column name is given.
      String colName = parset.getString (prefix+"datacolumn",
                                         reader->dataColumnName());
      return colName != reader->dataColumnName();
    }

    bool MSUpdater::updateAllowed (DPInfo& info, MSReader* reader,
                                    bool throwError) {
      if (info.nchanAvg() != 1 || info.ntimeAvg() != 1) {
        if (throwError) {
          THROW(Exception, "A new MS has to be given in msout if averaging is done");
        }
        return false;
      }
      if (!info.phaseCenterIsOriginal()) {
        if (throwError) {
          THROW(Exception, "A new MS has to be given in msout if a phase shift is done");
        }
        return false;
      }
      if (info.antennaNames().size() != reader->getInfo().antennaNames().size() ||
          ! allEQ(info.antennaNames(), reader->getInfo().antennaNames())) {
        if (throwError) {
          THROW(Exception, "A new MS has to be given if antennas are added or removed");
        }
        return false;
      }
      return true;
    }

    bool MSUpdater::addColumn(const string& colName, const casa::DataType
        dataType, const ColumnDesc& cd) {
      Table& tab = itsReader->table();
      if (tab.tableDesc().isColumn(colName)) {
        const ColumnDesc& cd = tab.tableDesc().columnDesc(colName);
        ASSERTSTR (cd.dataType() == dataType  &&  cd.isArray(),
                   "Column " << itsDataColName
                   << " already exists, but is not of the right type");
        return false;
      }

      TableDesc td;
      td.addColumn (cd, colName);

      // Use the same data manager as the DATA column.
      // Get the data manager info and find the DATA column in it.
      Record dminfo = tab.dataManagerInfo();
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
      tab.addColumn (td, colinfo);
      return true;
    }

    bool MSUpdater::process (const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      putFlags (buf.getRowNrs(), buf.getFlags());
      if (itsWriteData) {
        putData (buf.getRowNrs(), buf.getData());
      }
      if (itsWriteWeight) {
        putWeights (buf.getRowNrs(),
            itsReader->fetchWeights(buf, buf.getRowNrs(), itsTimer));
      }
      itsNrDone++;
      if (itsNrTimesFlush > 0  &&  itsNrDone%itsNrTimesFlush == 0) {
        itsReader->table().flush();
      }
      return true;
    }

    void MSUpdater::finish()
    {}

    void MSUpdater::show (std::ostream& os) const
    {
      os << "MSUpdater" << std::endl;
      os << "  MS:             " << itsReader->msName() << std::endl;
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
      os << "  flush:          " << itsNrTimesFlush << std::endl;
    }

    void MSUpdater::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " MSUpdater" << endl;
    }

    void MSUpdater::putFlags (const RefRows& rowNrs,
                              const Cube<bool>& flags)
    {
      // Only put if rownrs are filled, thus if data were not inserted.
      if (! rowNrs.rowVector().empty()) {
        Slicer colSlicer(IPosition(2, 0, info().startchan()),
                         IPosition(2, info().ncorr(), info().nchan()) );
        ArrayColumn<bool> flagCol(itsReader->table(), "FLAG");
        ScalarColumn<bool> flagRowCol(itsReader->table(), "FLAG_ROW");
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
        ArrayColumn<float> weightCol(itsReader->table(), itsWeightColName);
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
      //cout << "Data shape"<< data.shape()<<endl;
      if (! rowNrs.rowVector().empty()) {
        Slicer colSlicer(IPosition(2, 0, info().startchan()),
                         IPosition(2, info().ncorr(), info().nchan()) );
        ArrayColumn<Complex> dataCol(itsReader->table(), itsDataColName);
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
