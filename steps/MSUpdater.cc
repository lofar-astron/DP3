// MSUpdater.cc: DPPP step updating an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MSUpdater.h"

#include <cassert>
#include <iostream>
#include <limits>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledColumnStMan.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/LinearSearch.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <dp3/base/DPBuffer.h>

#include "../common/ParameterSet.h"

#include "InputStep.h"
#include "MSWriter.h"

using casacore::ArrayColumn;
using casacore::ColumnDesc;
using casacore::Cube;
using casacore::DataManager;
using casacore::DataManagerCtor;
using casacore::IPosition;
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

MSUpdater::MSUpdater(std::string msName, const common::ParameterSet& parset,
                     const std::string& prefix, bool writeHistory)
    : itsName(prefix),
      itsMSName(std::move(msName)),
      itsParset(parset),
      itsDataColName(parset.getString(prefix + "datacolumn", "")),
      itsFlagColName(parset.getString(prefix + "flagcolumn", "")),
      itsWeightColName(parset.getString(prefix + "weightcolumn", "")),
      itsNrTimesFlush(parset.getUint(prefix + "flush", 0)),
      itsNrDone(0),
      itsDataColAdded(false),
      itsFlagColAdded(false),
      itsWeightColAdded(false),
      itsWriteHistory(writeHistory),
      itsTileSize(parset.getUint(prefix + "tilesize", 1024)),
      itsStManKeys(parset, prefix) {
  // Call SetFieldsToWrite, since getRequiredFields() uses GetFieldsToWrite().
  SetFieldsToWrite(common::Fields());
}

bool MSUpdater::addColumn(const std::string& colName,
                          const casacore::DataType dataType,
                          const ColumnDesc& cd) {
  if (itsMS.tableDesc().isColumn(colName)) {
    const ColumnDesc& cd = itsMS.tableDesc().columnDesc(colName);
    if (cd.dataType() != dataType || !cd.isArray())
      throw std::runtime_error("Column " + colName +
                               " already exists, but is not of the right type");
    return false;
  }

  const bool use_sisco = itsStManKeys.storage_manager_name == "sisco" &&
                         dataType == casacore::TpComplex;
  const bool use_stokes_i = itsStManKeys.storage_manager_name == "stokes_i";
  const bool use_dysco = itsStManKeys.storage_manager_name == "dysco" &&
                         itsStManKeys.dysco_data_bit_rate != 0 &&
                         dataType != casacore::TpBool;
  if (use_sisco || use_stokes_i || use_dysco) {
    if (use_sisco && getInfoOut().origNChan() != getInfoOut().nchan())
      throw std::runtime_error(
          "Sisco can only write the entire bandwidth, subselecting channels is "
          "not possible");
    const std::string stdman_class_name(
        itsStManKeys.GetStorageManagerClassName());
    const casacore::DataManagerCtor data_constructor =
        DataManager::getCtor(stdman_class_name);
    if (!data_constructor)
      throw std::runtime_error(
          stdman_class_name +
          " storage manager requested, but it is not available in "
          "casacore");
    ColumnDesc directColumnDesc(cd);
    directColumnDesc.setOptions(casacore::ColumnDesc::Direct |
                                casacore::ColumnDesc::FixedShape);
    TableDesc td;
    td.addColumn(directColumnDesc, colName);
    std::unique_ptr<DataManager> data_manager(
        data_constructor(colName + "_dm", itsStManKeys.GetSpecification()));
    itsMS.addColumn(td, *data_manager);
  } else if (dataType == casacore::TpBool) {
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
      IPosition tileShape(3, getInfoOut().ncorr(), getInfoOut().nchan(), 1);
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

bool MSUpdater::process(std::unique_ptr<DPBuffer> buffer) {
  {
    common::NSTimer::StartStop sstime(itsTimer);
    const casacore::IPosition casacore_shape(3, getInfoOut().ncorr(),
                                             getInfoOut().nchan(),
                                             getInfoOut().nbaselines());
    if (GetFieldsToWrite().Flags()) {
      const Cube<bool> flags(casacore_shape, buffer->GetFlags().data(),
                             casacore::SHARE);
      putFlags(buffer->GetRowNumbers(), flags);
    }
    if (GetFieldsToWrite().Data()) {
      const Cube<casacore::Complex> data(
          casacore_shape, buffer->GetData().data(), casacore::SHARE);

      // If compressing, flagged values need to be set to NaN to decrease the
      // dynamic range
      if (itsStManKeys.storage_manager_name == "dysco") {
        Cube<casacore::Complex> data_copy = data.copy();
        Cube<casacore::Complex>::iterator data_iterator = data_copy.begin();
        for (const bool flag : buffer->GetFlags()) {
          if (flag) {
            *data_iterator =
                casacore::Complex(std::numeric_limits<float>::quiet_NaN(),
                                  std::numeric_limits<float>::quiet_NaN());
          }
          ++data_iterator;
        }
        putData(buffer->GetRowNumbers(), data_copy);
      } else {
        putData(buffer->GetRowNumbers(), data);
      }
    }
    if (GetFieldsToWrite().Weights()) {
      const Cube<float> weights(casacore_shape, buffer->GetWeights().data(),
                                casacore::SHARE);

      // If compressing, set weights of flagged points to zero to decrease the
      // dynamic range
      if (itsStManKeys.storage_manager_name == "dysco") {
        Cube<float> weights_copy = weights.copy();
        Cube<float>::iterator weights_iterator = weights_copy.begin();
        for (const bool flag : buffer->GetFlags()) {
          if (flag) {
            *weights_iterator = 0.0f;
          }
          ++weights_iterator;
        }
        putWeights(buffer->GetRowNumbers(), weights_copy);
      } else {
        putWeights(buffer->GetRowNumbers(), weights);
      }
    }
    itsNrDone++;
    if (itsNrTimesFlush > 0 && itsNrDone % itsNrTimesFlush == 0) {
      itsMS.flush();
    }
  }
  getNextStep()->process(std::move(buffer));
  return true;
}

void MSUpdater::finish() {
  addToMS(itsMSName);
  if (itsWriteHistory) {
    MSWriter::WriteHistory(itsMS, itsParset);
  }
  if (getNextStep()) getNextStep()->finish();
}

common::Fields MSUpdater::getRequiredFields() const {
  common::Fields fields;
  if (GetFieldsToWrite().Data()) fields |= kDataField;
  if (GetFieldsToWrite().Flags()) fields |= kFlagsField;
  if (GetFieldsToWrite().Weights()) fields |= kWeightsField;
  return fields;
}

void MSUpdater::SetFieldsToWrite(const common::Fields& base_fields) {
  // A non-empty column name indicates it should be written.
  common::Fields fields = base_fields;
  if (!itsDataColName.empty()) fields |= kDataField;
  if (!itsFlagColName.empty()) fields |= kFlagsField;
  if (!itsWeightColName.empty()) fields |= kWeightsField;
  OutputStep::SetFieldsToWrite(fields);
}

void MSUpdater::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  if (getInfoOut().metaChanged()) {
    throw std::runtime_error("Update step " + itsName +
                             " is not possible because meta data changes"
                             " (by averaging, adding/removing stations, etc.)");
  }

  // Since an MSUpdater supports updating an MS created by MSWriter, and
  // MSWriter creates its MS in updateInfo(), MSUpdater should delay opening the
  // MS until updateInfo().
  itsMS = casacore::MeasurementSet(itsMSName, TableLock::AutoNoReadLocking,
                                   Table::Update);
  if (InputStep::HasBda(itsMS)) {
    throw std::runtime_error(
        R"(Update step is not possible because the input/output types are incompatible
(BDA buffer - Regular buffer).
Specify a name in the parset for "msout".)");
  }

  if (itsDataColName.empty()) {
    itsDataColName = infoIn.dataColumnName();
  }

  if (itsWeightColName.empty()) {
    if (infoIn.weightColumnName() == "WEIGHT") {
      itsWeightColName = "WEIGHT_SPECTRUM";
      SetFieldsToWrite(GetFieldsToWrite() | kWeightsField);
    } else {
      itsWeightColName = infoIn.weightColumnName();
    }
  }

  if (itsWeightColName == "WEIGHT")
    throw std::runtime_error(
        "Can't use WEIGHT column as spectral weights column");

  if (itsFlagColName.empty()) {
    itsFlagColName = infoIn.flagColumnName();
  }

  if (GetFieldsToWrite().Data() || GetFieldsToWrite().Flags() ||
      GetFieldsToWrite().Weights()) {
    common::NSTimer::StartStop sstime(itsTimer);
    // Add the data + flag + weight column if needed and if it does not exist
    // yet.
    if (GetFieldsToWrite().Data()) {
      // use same layout as DATA column
      ColumnDesc cd = itsMS.tableDesc().columnDesc("DATA");
      itsDataColAdded = addColumn(itsDataColName, casacore::TpComplex, cd);
    }
    if (GetFieldsToWrite().Flags()) {
      // use same layout as FLAG column
      ColumnDesc cd = itsMS.tableDesc().columnDesc("FLAG");
      itsFlagColAdded = addColumn(itsFlagColName, casacore::TpBool, cd);
    }
    if (GetFieldsToWrite().Weights()) {
      IPosition dataShape = itsMS.tableDesc().columnDesc("DATA").shape();
      casacore::ArrayColumnDesc<float> cd("WEIGHT_SPECTRUM",
                                          "weight per corr/chan", dataShape,
                                          ColumnDesc::FixedShape);
      itsWeightColAdded = addColumn(itsWeightColName, casacore::TpFloat, cd);
    }
  }
  MSWriter::UpdateBeam(itsMS, itsDataColName, getInfoOut());
  // Subsequent steps have to set again if writes need to be done.
  GetWritableInfoOut().clearMetaChanged();
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
  if (GetFieldsToWrite().Data() || GetFieldsToWrite().Flags() ||
      GetFieldsToWrite().Weights()) {
    os << "  writing:       ";
    if (GetFieldsToWrite().Data()) os << " data";
    if (GetFieldsToWrite().Flags()) os << " flags";
    if (GetFieldsToWrite().Weights()) os << " weights";
    os << '\n';
  }
  os << GetCompressionString(itsStManKeys);
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
    Slicer colSlicer(IPosition(2, 0, getInfoOut().startchan()),
                     IPosition(2, getInfoOut().ncorr(), getInfoOut().nchan()));
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
    Slicer colSlicer(IPosition(2, 0, getInfoOut().startchan()),
                     IPosition(2, getInfoOut().ncorr(), getInfoOut().nchan()));
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

void MSUpdater::putData(const RefRows& row_nrs,
                        const Cube<casacore::Complex>& data) {
  // Only put if rownrs are filled, thus if data were not inserted.
  if (!row_nrs.rowVector().empty()) {
    ArrayColumn<casacore::Complex> data_col(itsMS, itsDataColName);
    // Loop over all rows of this subset.
    // (it also avoids StandardStMan putCol with RefRows problem).
    casacore::Vector<common::rownr_t> rows = row_nrs.convert();
    casacore::ReadOnlyArrayIterator<casacore::Complex> dataIter(data, 2);
    if (getInfoOut().origNChan() == getInfoOut().nchan()) {
      // If all channels are written, using putSlice() is unnecessary. This is
      // important for Sisco, which doesn't support putSlice().
      for (size_t i = 0; i < rows.size(); ++i) {
        data_col.put(rows[i], dataIter.array());
        dataIter.next();
      }
    } else {
      Slicer col_slicer(
          IPosition(2, 0, getInfoOut().startchan()),
          IPosition(2, getInfoOut().ncorr(), getInfoOut().nchan()));
      for (size_t i = 0; i < rows.size(); ++i) {
        data_col.putSlice(rows[i], col_slicer, dataIter.array());
        dataIter.next();
      }
    }
  }
}
}  // namespace steps
}  // namespace dp3
