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

#include "MSUpdater.h"
#include "MSReader.h"
#include "MSWriter.h"
#include "DPBuffer.h"

#include "../Common/ParameterSet.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/LinearSearch.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <iostream>
#include <limits>

using namespace casacore;

namespace DP3 {
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
      itsTileSize      = parset.getUint (prefix+"tilesize", 1024);
      itsStManKeys.Set(parset, prefix);
    }

    MSUpdater::~MSUpdater()
    {}

    bool MSUpdater::addColumn (const string& colName, const casacore::DataType
                               dataType, const ColumnDesc& cd)
    {
      if (itsMS.tableDesc().isColumn(colName)) {
        const ColumnDesc& cd = itsMS.tableDesc().columnDesc(colName);
        if (cd.dataType() != dataType || !cd.isArray())
                   "Column " + itsDataColName +
                   " already exists, but is not of the right type";
        return false;
      }

      if (itsStManKeys.stManName == "dysco" && itsStManKeys.dyscoDataBitRate != 0) {
        casacore::Record dyscoSpec = itsStManKeys.GetDyscoSpec();
        DataManagerCtor dyscoConstructor = DataManager::getCtor("DyscoStMan");
        std::unique_ptr<DataManager> dyscoStMan(dyscoConstructor(colName + "_dm", dyscoSpec));
        ColumnDesc directColumnDesc(cd);
        directColumnDesc.setOptions(casacore::ColumnDesc::Direct | casacore::ColumnDesc::FixedShape);
        TableDesc td;
        td.addColumn (directColumnDesc, colName);
        itsMS.addColumn (td, *dyscoStMan);
      }
      else {
        // When no specific storage manager is requested, use the same
        // as for the DATA column.
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
        assert(colinfo.nfields()>0);
        // When the storage manager is compressed, do not implicitly (re)compress it. Use TiledStMan instead.
        std::string dmType = colinfo.asString("TYPE");
        TableDesc td;
        td.addColumn (cd, colName);
        if(dmType == "DyscoStMan")
        {
          IPosition tileShape(3, info().ncorr(), info().nchan(), 1);
          tileShape[2] = itsTileSize * 1024 / (8 * tileShape[0] * tileShape[1]);
          if (tileShape[2] < 1) {
            tileShape[2] = 1;
          }
          TiledColumnStMan tsm(colName + "_dm", tileShape);
          itsMS.addColumn (td, tsm);
        }
        else {
          colinfo.define ("NAME", colName + "_dm");
          itsMS.addColumn (td, colinfo);
        }
      }
      return true;
    }

    bool MSUpdater::process (const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      if (itsWriteFlags) {
        putFlags (buf.getRowNrs(), buf.getFlags());
      }
      if (itsWriteData) {
        // If compressing, flagged values need to be set to NaN to decrease the dynamic range
        if (itsStManKeys.stManName == "dysco") {
          Cube<Complex> dataCopy = buf.getData().copy();
          Cube<Complex>::iterator dataIter = dataCopy.begin();
          for (Cube<bool>::const_iterator flagIter = buf.getFlags().begin();
               flagIter != buf.getFlags().end(); ++flagIter) {
            if(*flagIter) {
              *dataIter = Complex(std::numeric_limits<float>::quiet_NaN(),
                                  std::numeric_limits<float>::quiet_NaN());
            }
            ++dataIter;
          }
          putData (buf.getRowNrs(), dataCopy);
        }
        else {
          putData (buf.getRowNrs(), buf.getData());
        }
      }
      if (itsWriteWeights) {
        Cube<float> weights;
        if (!buf.getWeights().empty()) {
          // Use weights from buffer
          weights = buf.getWeights();
        } else {
          itsBuffer.referenceFilled (buf);
          weights = itsReader->fetchWeights(buf, itsBuffer, itsTimer);
        }

        // If compressing, set weights of flagged points to zero to decrease the
        // dynamic range
        if (itsStManKeys.stManName == "dysco") {
          Cube<float> weightsCopy = weights.copy();
          Cube<float>::iterator weightsIter = weightsCopy.begin();
          for (Cube<bool>::const_iterator flagIter = buf.getFlags().begin();
               flagIter != buf.getFlags().end(); ++flagIter) {
            if(*flagIter) {
              *weightsIter = 0.;
            }
            ++weightsIter;
          }
          putWeights (buf.getRowNrs(), weightsCopy);
        } else {
          putWeights (buf.getRowNrs(), weights);
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
      assert(itsWeightColName != "WEIGHT");
      if (itsWeightColName != origWeightColName) {
        info().setWriteWeights();
      }
      itsWriteWeights = getInfo().writeWeights();

      if (getInfo().metaChanged()) {
        throw Exception("Update step " + itsName + 
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
      MSWriter::updateBeam(itsMSName, itsDataColName, info());
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
      if(itsStManKeys.stManName == "dysco") {
        os
          << "  Compressed:     yes\n"
          << "  Data bitrate:   " << itsStManKeys.dyscoDataBitRate << '\n'
          << "  Weight bitrate: " << itsStManKeys.dyscoWeightBitRate << '\n'
          << "  Dysco mode:     " << itsStManKeys.dyscoNormalization << ' '
          << itsStManKeys.dyscoDistribution << '(' << itsStManKeys.dyscoDistTruncation << ")\n";
      }
      else {
        os << "  Compressed:     no\n";
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
