// MSUpdater.cc: DPPP step updating an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MSUpdater.h"
#include "MSReader.h"
#include "MSWriter.h"

#include "../base/DPBuffer.h"

#include "../common/ParameterSet.h"

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledColumnStMan.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/LinearSearch.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <iostream>
#include <limits>

using casacore::ArrayColumn;
using casacore::ColumnDesc;
using casacore::Cube;
using casacore::DataManager;
using casacore::DataManagerCtor;
using casacore::IPosition;
using casacore::MeasurementSet;
using casacore::Record;
using casacore::RefRows;
using casacore::ScalarColumn;
using casacore::Slicer;
using casacore::Table;
using casacore::TableDesc;
using casacore::TableLock;
using casacore::TiledColumnStMan;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

MSUpdater::MSUpdater(InputStep* reader, casacore::String msName,
                     const common::ParameterSet& parset,
                     const std::string& prefix, bool writeHistory)
    : itsReader(reader),
      itsName(prefix),
      itsMSName(msName),
      itsParset(parset),
      itsWriteData(false),
      itsWriteWeights(false),
      itsWriteFlags(false),
      itsNrDone(0),
      itsDataColAdded(false),
      itsFlagColAdded(false),
      itsWeightColAdded(false),
      itsWriteHistory(writeHistory) {
  itsDataColName = parset.getString(prefix + "datacolumn", "");
  itsFlagColName = parset.getString(prefix + "flagcolumn", "");
  itsWeightColName = parset.getString(prefix + "weightcolumn", "");
  itsNrTimesFlush = parset.getUint(prefix + "flush", 0);
  itsTileSize = parset.getUint(prefix + "tilesize", 1024);
  itsStManKeys.Set(parset, prefix);
}

MSUpdater::~MSUpdater() {}

bool MSUpdater::addColumn(const string& colName,
                          const casacore::DataType dataType,
                          const ColumnDesc& cd) {
  if (itsMS.tableDesc().isColumn(colName)) {
    const ColumnDesc& cd = itsMS.tableDesc().columnDesc(colName);
    if (cd.dataType() != dataType || !cd.isArray())
      throw std::runtime_error("Column " + colName +
                               " already exists, but is not of the right type");
    return false;
  }

  if (dataType == casacore::TpBool) {
    // Dysco should never be used for the FLAG column. Use the same storage
    // manager used for the FLAG column. To do so, get the data manager info and
    // find the FLAG column in it.
    Record dminfo = itsMS.dataManagerInfo();
    Record colinfo;
    for (size_t i = 0; i < dminfo.nfields(); ++i) {
      const Record& subrec = dminfo.subRecord(i);
      if (linearSearch1(casacore::Vector<casacore::String>(
                            subrec.asArrayString("COLUMNS")),
                        "FLAG") >= 0) {
        colinfo = subrec;
        break;
      }
    }
    if (colinfo.nfields() == 0)
      throw std::runtime_error("Could not obtain column info");
    TableDesc td;
    td.addColumn(cd, colName);
    colinfo.define("NAME", colName + "_dm");
    itsMS.addColumn(td, colinfo);
  } else if (itsStManKeys.stManName == "dysco" &&
             itsStManKeys.dyscoDataBitRate != 0) {
    casacore::Record dyscoSpec = itsStManKeys.GetDyscoSpec();
    DataManagerCtor dyscoConstructor = DataManager::getCtor("DyscoStMan");
    std::unique_ptr<DataManager> dyscoStMan(
        dyscoConstructor(colName + "_dm", dyscoSpec));
    ColumnDesc directColumnDesc(cd);
    directColumnDesc.setOptions(casacore::ColumnDesc::Direct |
                                casacore::ColumnDesc::FixedShape);
    TableDesc td;
    td.addColumn(directColumnDesc, colName);
    itsMS.addColumn(td, *dyscoStMan);
  } else {
    // When no specific storage manager is requested, use the same
    // as for the DATA column.
    // Get the data manager info and find the DATA column in it.
    Record dminfo = itsMS.dataManagerInfo();
    Record colinfo;
    for (size_t i = 0; i < dminfo.nfields(); ++i) {
      const Record& subrec = dminfo.subRecord(i);
      if (linearSearch1(casacore::Vector<casacore::String>(
                            subrec.asArrayString("COLUMNS")),
                        "DATA") >= 0) {
        colinfo = subrec;
        break;
      }
    }
    if (colinfo.nfields() == 0)
      throw std::runtime_error("Could not obtain column info");
    // When the storage manager is compressed, do not implicitly (re)compress
    // it. Use TiledStMan instead.
    std::string dmType = colinfo.asString("TYPE");
    TableDesc td;
    td.addColumn(cd, colName);
    if (dmType == "DyscoStMan") {
      IPosition tileShape(3, info().ncorr(), info().nchan(), 1);
      tileShape[2] = itsTileSize * 1024 / (8 * tileShape[0] * tileShape[1]);
      if (tileShape[2] < 1) {
        tileShape[2] = 1;
      }
      TiledColumnStMan tsm(colName + "_dm", tileShape);
      itsMS.addColumn(td, tsm);
    } else {
      colinfo.define("NAME", colName + "_dm");
      itsMS.addColumn(td, colinfo);
    }
  }
  return true;
}

bool MSUpdater::process(const DPBuffer& buf) {
  common::NSTimer::StartStop sstime(itsTimer);
  if (itsWriteFlags) {
    putFlags(buf.getRowNrs(), buf.getFlags());
  }
  if (itsWriteData) {
    // If compressing, flagged values need to be set to NaN to decrease the
    // dynamic range
    if (itsStManKeys.stManName == "dysco") {
      Cube<casacore::Complex> dataCopy = buf.getData().copy();
      Cube<casacore::Complex>::iterator dataIter = dataCopy.begin();
      for (Cube<bool>::const_iterator flagIter = buf.getFlags().begin();
           flagIter != buf.getFlags().end(); ++flagIter) {
        if (*flagIter) {
          *dataIter =
              casacore::Complex(std::numeric_limits<float>::quiet_NaN(),
                                std::numeric_limits<float>::quiet_NaN());
        }
        ++dataIter;
      }
      putData(buf.getRowNrs(), dataCopy);
    } else {
      putData(buf.getRowNrs(), buf.getData());
    }
  }
  if (itsWriteWeights) {
    Cube<float> weights;
    if (!buf.getWeights().empty()) {
      // Use weights from buffer
      weights = buf.getWeights();
    } else {
      itsBuffer.referenceFilled(buf);
      weights = itsReader->fetchWeights(buf, itsBuffer, itsTimer);
    }

    // If compressing, set weights of flagged points to zero to decrease the
    // dynamic range
    if (itsStManKeys.stManName == "dysco") {
      Cube<float> weightsCopy = weights.copy();
      Cube<float>::iterator weightsIter = weightsCopy.begin();
      for (Cube<bool>::const_iterator flagIter = buf.getFlags().begin();
           flagIter != buf.getFlags().end(); ++flagIter) {
        if (*flagIter) {
          *weightsIter = 0.;
        }
        ++weightsIter;
      }
      putWeights(buf.getRowNrs(), weightsCopy);
    } else {
      putWeights(buf.getRowNrs(), weights);
    }
  }
  itsNrDone++;
  if (itsNrTimesFlush > 0 && itsNrDone % itsNrTimesFlush == 0) {
    itsMS.flush();
  }
  getNextStep()->process(buf);
  return true;
}

void MSUpdater::finish() {}

void MSUpdater::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  if (itsReader->outputs() != this->outputs()) {
    throw Exception(
        "Update step is not possible because input/output types are "
        "incompatible (BDA buffer - Regular buffer).\nSpecify a name in the "
        "parset for \"msout\"");
  }

  std::string origDataColName = getInfo().getDataColName();
  if (itsDataColName.empty()) {
    itsDataColName = origDataColName;
  } else if (itsDataColName != origDataColName) {
    info().setNeedVisData();
    info().setWriteData();
  }
  itsWriteData = getInfo().writeData();

  std::string origWeightColName = getInfo().getWeightColName();
  if (itsWeightColName.empty()) {
    if (origWeightColName == "WEIGHT") {
      itsWeightColName = "WEIGHT_SPECTRUM";
    } else {
      itsWeightColName = origWeightColName;
    }
  }
  if (itsWeightColName == "WEIGHT")
    throw std::runtime_error(
        "Can't use WEIGHT column as spectral weights column");
  if (itsWeightColName != origWeightColName) {
    info().setWriteWeights();
  }
  itsWriteWeights = getInfo().writeWeights();

  std::string origFlagColName = getInfo().getFlagColName();
  if (itsFlagColName.empty()) {
    itsFlagColName = origFlagColName;
  } else if (itsFlagColName != origFlagColName) {
    info().setWriteFlags();
  }

  itsWriteFlags = getInfo().writeFlags();

  if (getInfo().metaChanged()) {
    throw Exception("Update step " + itsName +
                    " is not possible because meta data changes"
                    " (by averaging, adding/removing stations, etc.)");
  }

  if (itsWriteData || itsWriteFlags || itsWriteWeights) {
    common::NSTimer::StartStop sstime(itsTimer);
    itsMS =
        MeasurementSet(itsMSName, TableLock::AutoNoReadLocking, Table::Update);
    // Add the data + flag + weight column if needed and if it does not exist
    // yet.
    if (itsWriteData) {
      // use same layout as DATA column
      ColumnDesc cd = itsMS.tableDesc().columnDesc("DATA");
      itsDataColAdded = addColumn(itsDataColName, casacore::TpComplex, cd);
    }
    if (itsWriteFlags) {
      // use same layout as FLAG column
      ColumnDesc cd = itsMS.tableDesc().columnDesc("FLAG");
      itsFlagColAdded = addColumn(itsFlagColName, casacore::TpBool, cd);
    }
    if (itsWriteWeights) {
      IPosition dataShape = itsMS.tableDesc().columnDesc("DATA").shape();
      casacore::ArrayColumnDesc<float> cd("WEIGHT_SPECTRUM",
                                          "weight per corr/chan", dataShape,
                                          ColumnDesc::FixedShape);
      itsWeightColAdded = addColumn(itsWeightColName, casacore::TpFloat, cd);
    }
  }
  MSWriter::updateBeam(itsMSName, itsDataColName, info());
  // Subsequent steps have to set again if writes need to be done.
  info().clearWrites();
  info().clearMetaChanged();
}

void MSUpdater::addToMS(const string&) {
  getPrevStep()->addToMS(itsMSName);
  if (itsWriteHistory) {
    MSWriter::writeHistory(itsMS, itsParset);
  }
}

void MSUpdater::show(std::ostream& os) const {
  os << "MSUpdater " << itsName << '\n';
  os << "  MS:             " << itsMSName << '\n';
  os << "  datacolumn:     " << itsDataColName;
  if (itsDataColAdded) {
    os << "  (has been added to the MS)";
  }
  os << '\n';
  os << "  flagcolumn:     " << itsFlagColName;
  if (itsFlagColAdded) {
    os << "  (has been added to the MS)";
  }
  os << '\n';
  os << "  weightcolumn    " << itsWeightColName;
  if (itsWeightColAdded) {
    os << "  (has been added to the MS)";
  }
  os << '\n';
  if (itsWriteData || itsWriteFlags || itsWriteWeights) {
    os << "  writing:       ";
    if (itsWriteData) os << " data";
    if (itsWriteFlags) os << " flags";
    if (itsWriteWeights) os << " weights";
    os << '\n';
  }
  if (itsStManKeys.stManName == "dysco") {
    os << "  Compressed:     yes\n"
       << "  Data bitrate:   " << itsStManKeys.dyscoDataBitRate << '\n'
       << "  Weight bitrate: " << itsStManKeys.dyscoWeightBitRate << '\n'
       << "  Dysco mode:     " << itsStManKeys.dyscoNormalization << ' '
       << itsStManKeys.dyscoDistribution << '('
       << itsStManKeys.dyscoDistTruncation << ")\n";
  } else {
    os << "  Compressed:     no\n";
  }
  os << '\n';
  os << "  flush:          " << itsNrTimesFlush << '\n';
}

void MSUpdater::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " MSUpdater " << itsName << '\n';
}

void MSUpdater::putFlags(const RefRows& rowNrs, const Cube<bool>& flags) {
  // Only put if rownrs are filled, thus if data were not inserted.
  if (!rowNrs.rowVector().empty()) {
    Slicer colSlicer(IPosition(2, 0, info().startchan()),
                     IPosition(2, info().ncorr(), info().nchan()));
    ArrayColumn<bool> flagCol(itsMS, itsFlagColName);
    ScalarColumn<bool> flagRowCol(itsMS, "FLAG_ROW");
    // Loop over all rows of this subset.
    // (it also avoids StandardStMan putCol with RefRows problem).
    casacore::Vector<common::rownr_t> rows = rowNrs.convert();
    casacore::ReadOnlyArrayIterator<bool> flagIter(flags, 2);
    for (common::rownr_t i = 0; i < rows.size(); ++i) {
      flagCol.putSlice(rows[i], colSlicer, flagIter.array());
      // If a new flag in a row is clear, the ROW_FLAG should not be set.
      // If all new flags are set, we leave it because we might have a
      // subset of the channels, so other flags might still be clear.
      if (anyEQ(flagIter.array(), false)) {
        flagRowCol.put(rows[i], false);
      }
      flagIter.next();
    }
  }
}

void MSUpdater::putWeights(const RefRows& rowNrs, const Cube<float>& weights) {
  // Only put if rownrs are filled, thus if data were not inserted.
  if (!rowNrs.rowVector().empty()) {
    Slicer colSlicer(IPosition(2, 0, info().startchan()),
                     IPosition(2, info().ncorr(), info().nchan()));
    ArrayColumn<float> weightCol(itsMS, itsWeightColName);
    // Loop over all rows of this subset.
    // (it also avoids StandardStMan putCol with RefRows problem).
    casacore::Vector<common::rownr_t> rows = rowNrs.convert();
    casacore::ReadOnlyArrayIterator<float> weightIter(weights, 2);
    for (size_t i = 0; i < rows.size(); ++i) {
      weightCol.putSlice(rows[i], colSlicer, weightIter.array());
      weightIter.next();
    }
  }
}

void MSUpdater::putData(const RefRows& rowNrs,
                        const Cube<casacore::Complex>& data) {
  // Only put if rownrs are filled, thus if data were not inserted.
  if (!rowNrs.rowVector().empty()) {
    Slicer colSlicer(IPosition(2, 0, info().startchan()),
                     IPosition(2, info().ncorr(), info().nchan()));
    ArrayColumn<casacore::Complex> dataCol(itsMS, itsDataColName);
    // Loop over all rows of this subset.
    // (it also avoids StandardStMan putCol with RefRows problem).
    casacore::Vector<common::rownr_t> rows = rowNrs.convert();
    casacore::ReadOnlyArrayIterator<casacore::Complex> dataIter(data, 2);
    for (size_t i = 0; i < rows.size(); ++i) {
      dataCol.putSlice(rows[i], colSlicer, dataIter.array());
      dataIter.next();
    }
  }
}
}  // namespace steps
}  // namespace dp3
