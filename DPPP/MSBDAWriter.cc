// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/DataMan/IncrementalStMan.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableCopy.h>
#include <casacore/tables/Tables/TableDesc.h>

#include "../Common/ParameterSet.h"
#include "BDABuffer.h"
#include "MSBDAWriter.h"
#include "MSWriter.h"

using casacore::ArrayColumn;
using casacore::ArrayColumnDesc;
using casacore::Bool;
using casacore::ColumnDesc;
using casacore::Double;
using casacore::Int;
using casacore::MeasurementSet;
using casacore::ObjectID;
using casacore::ScalarColumn;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::Table;
using casacore::TableCopy;
using casacore::TableDesc;
using casacore::TableLock;
using casacore::True;
using casacore::Vector;

namespace {
/// Measurement Set column names (main table). Uses the same order as
/// https://casacore.github.io/casacore-notes/229.html#x1-630005.1
/// @{
const std::string kTime = "TIME";
const std::string kAntenna1 = "ANTENNA1";
const std::string kAntenna2 = "ANTENNA2";
const std::string kFeed1 = "FEED1";
const std::string kFeed2 = "FEED2";
const std::string kDataDescId = "DATA_DESC_ID";
const std::string kProcessorId = "PROCESSOR_ID";
const std::string kFieldId = "FIELD_ID";
const std::string kInterval = "INTERVAL";
const std::string kExposure = "EXPOSURE";
const std::string kTimeCentroid = "TIME_CENTROID";
const std::string kScanNumber = "SCAN_NUMBER";
const std::string kArrayId = "ARRAY_ID";
const std::string kObservationId = "OBSERVATION_ID";
const std::string kStateId = "STATE_ID";
const std::string kUVW = "UVW";
const std::string kData = "DATA";
const std::string kSigma = "SIGMA";
const std::string kWeight = "WEIGHT";
const std::string kWeightSpectrum = "WEIGHT_SPECTRUM";
const std::string kFlag = "FLAG";
const std::string kFlagCategory = "FLAG_CATEGORY";  // Is ignored when writing.
const std::string kFlagRow = "FLAG_ROW";
/// @}
}  // namespace

namespace {
/// BDA_TIME_AXIS table column names.
const std::string kBDATimeAxisTable = "BDA_TIME_AXIS";

const std::string kTimeAxisId = "BDA_TIME_AXIS_ID";
const std::string kIsBdaApplied = "IS_BDA_APPLIED";
const std::string kSingleFactorPerBL = "SINGLE_FACTOR_PER_BASELINE";
const std::string kMaxTimeInterval = "MAX_TIME_INTERVAL";
const std::string kMinTimeInterval = "MIN_TIME_INTERVAL";
const std::string kUnitTimeInterval = "UNIT_TIME_INTERVAL";
const std::string kIntervalFactors = "INTEGER_INTERVAL_FACTORS";
const std::string kHasBDAOrdering = "HAS_BDA_ORDERING";
/// @}
}  // namespace

namespace {
/// BDA_TIME_FACTOR table column names.
const std::string kBDATimeFactorTable = "BDA_TIME_FACTOR";

const std::string kFactor = "FACTOR";
/// @}
}  // namespace

namespace {
/// BDA metadata table column names for SPECTRAL_WINDOW.
const std::string kSpectralWindowTable = "SPECTRAL_WINDOW";

const std::string kBDAFreqAxisId = "BDA_FREQ_AXIS_ID";
const std::string kBDASetId = "BDA_SET_ID";
/// @}
}  // namespace

namespace DP3 {
namespace DPPP {

MSBDAWriter::MSBDAWriter(MSReader* reader, const std::string& out_name,
                         const ParameterSet& parset, const std::string& prefix)
    : reader_(reader),
      outName_(out_name),
      parset_(parset),
      prefix_(prefix),
      overwrite_(parset.getBool(prefix + "overwrite", false)) {}

MSBDAWriter::~MSBDAWriter() {}

void MSBDAWriter::updateInfo(const DPInfo& info_in) {
  DPStep::updateInfo(info_in);
  CreateMS();

  WriteMetaData();

  MSWriter::writeHistory(ms_, parset_);
  ms_.flush(true, true);
}

bool MSBDAWriter::process(std::unique_ptr<BDABuffer> buffer) {
  buffer->SetBaseRowNr(ms_.nrow());

  const std::vector<BDABuffer::Row>& rows = buffer->GetRows();

  ms_.addRow(rows.size());

  ScalarColumn<casacore::Double> time(ms_, kTime);
  ScalarColumn<casacore::Double> time_centroid(ms_, kTimeCentroid);
  ScalarColumn<casacore::Double> exposure(ms_, kExposure);
  ScalarColumn<casacore::Int> ant1(ms_, kAntenna1);
  ScalarColumn<casacore::Int> ant2(ms_, kAntenna2);
  ArrayColumn<casacore::Complex> data(ms_, kData);
  ArrayColumn<casacore::Float> weights(ms_, kWeightSpectrum);
  ArrayColumn<casacore::Bool> flags(ms_, kFlag);
  ScalarColumn<casacore::Bool> flags_row(ms_, kFlagRow);
  ArrayColumn<casacore::Double> uvw(ms_, kUVW);
  std::vector<DP3::rownr_t> row_nrs;
  row_nrs.reserve(rows.size());
  for (const BDABuffer::Row& row : rows) {
    time.put(row.row_nr_, row.time_);
    time_centroid.put(row.row_nr_, row.time_);
    exposure.put(row.row_nr_, row.interval_);

    ant1.put(row.row_nr_, info().getAnt1()[row.baseline_nr_]);
    ant2.put(row.row_nr_, info().getAnt2()[row.baseline_nr_]);

    // TODO: When DPInfo is updated, remove the assert and the comment below.
    assert(row.baseline_nr_ == 0);
    const std::size_t n_chan =
        row.n_channels_;  // row.chanFreqs(row.baseline_nr_).size();
    const casacore::IPosition dim(2, info().ncorr(), n_chan);
    data.put(row.row_nr_, casacore::Array<casacore::Complex>(dim, row.data_,
                                                             casacore::SHARE));
    weights.put(row.row_nr_, casacore::Array<casacore::Float>(dim, row.weights_,
                                                              casacore::SHARE));
    flags.put(row.row_nr_, casacore::Array<casacore::Bool>(dim, row.flags_,
                                                           casacore::SHARE));
    // Set the row_flag if all flags in the row are true / none are false.
    const bool row_flag =
        std::count(row.flags_, row.flags_ + row.GetDataSize(), false) == 0;
    flags_row.put(row.row_nr_, row_flag);

    const casacore::IPosition uvw_dim(1, 3);
    uvw.put(row.row_nr_, casacore::Array<casacore::Double>(uvw_dim, row.uvw_));

    row_nrs.push_back(row.row_nr_);
  }

  // Create a table view containing only the added rows.
  casacore::Table tbl_added(ms_(row_nrs));

  // Fill various mandatory columns with default values.
  ScalarColumn<casacore::Int>(tbl_added, kFeed1).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kFeed2).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kDataDescId).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kProcessorId).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kFieldId).fillColumn(0);
  ScalarColumn<casacore::Double>(tbl_added, kInterval)
      .fillColumn(info().timeInterval());
  ScalarColumn<casacore::Int>(tbl_added, kScanNumber).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kArrayId).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kObservationId).fillColumn(0);
  ScalarColumn<casacore::Int>(tbl_added, kStateId).fillColumn(0);
  casacore::Vector<casacore::Float> sigma_weight(info().ncorr(), 1);
  ArrayColumn<casacore::Float>(tbl_added, kSigma).fillColumn(sigma_weight);
  ArrayColumn<casacore::Float>(tbl_added, kWeight).fillColumn(sigma_weight);

  return true;
}

void MSBDAWriter::finish() {}

void MSBDAWriter::show(std::ostream& os) const {}

void MSBDAWriter::CreateMS() {
  CreateMainTable();
  // TODO fill the main table

  CreateMetaDataFrequencyColumns();
  CreateBDATimeAxis();
  CreateBDATimeFactor();
}

void MSBDAWriter::CreateMainTable() {
  // Create an empty table.
  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable setup(outName_, TableDesc(), opt);
  ms_ = Table(setup);

  TableDesc desc_inc;
  TableDesc desc_std;

  // Add fixed fields with the incremental storage manager.
  for (const std::string& name :
       {kFeed1, kFeed2, kDataDescId, kProcessorId, kFieldId, kScanNumber,
        kArrayId, kObservationId, kStateId}) {
    ScalarColumnDesc<casacore::Int> column(name, ColumnDesc::FixedShape);
    desc_inc.addColumn(column);
  }
  for (const std::string& name : {kInterval}) {
    ScalarColumnDesc<casacore::Double> column(name, ColumnDesc::FixedShape);
    desc_inc.addColumn(column);
  }
  for (const std::string& name : {kSigma, kWeight}) {
    const casacore::IPosition dim(1, info().ncorr());
    ArrayColumnDesc<casacore::Float> column(name, dim, ColumnDesc::FixedShape);
    desc_inc.addColumn(column);
  }

  // Add fixed-shape fields with the standard storange manager.
  for (const std::string& name : {kTime, kExposure, kTimeCentroid}) {
    ScalarColumnDesc<casacore::Double> column(name, ColumnDesc::FixedShape);
    desc_std.addColumn(column);
  }
  for (const std::string& name : {kAntenna1, kAntenna2}) {
    ScalarColumnDesc<casacore::Int> column(name, ColumnDesc::FixedShape);
    desc_std.addColumn(column);
  }
  {
    ScalarColumnDesc<Bool> column(kFlagRow, ColumnDesc::FixedShape);
    desc_std.addColumn(column);
  }
  {
    const casacore::IPosition dim(1, 3);
    ArrayColumnDesc<casacore::Double> column(kUVW, dim, ColumnDesc::FixedShape);
    desc_std.addColumn(column);
  }

  // Add variable-shaped fields.
  ArrayColumnDesc<casacore::Complex> col_data(kData);
  ArrayColumnDesc<casacore::Float> col_weights(kWeightSpectrum);
  ArrayColumnDesc<casacore::Bool> col_flags(kFlag);
  desc_std.addColumn(col_data);
  desc_std.addColumn(col_weights);
  desc_std.addColumn(col_flags);

  casacore::IncrementalStMan man_inc;
  casacore::StandardStMan man_std(32768);
  ms_.addColumn(desc_inc, man_inc);
  ms_.addColumn(desc_std, man_std);

  if (reader_) {
    // Copy the info and subtables.
    TableCopy::copyInfo(ms_, reader_->table());

    casacore::Block<casacore::String> omitted_subtables(1);
    omitted_subtables[0] = "BDA_TIME_AXIS";
    TableCopy::copySubTables(ms_, reader_->table(), false, omitted_subtables);
  } else {
    // TODO: Create info and subtables.
    assert(false);
  }
}  // namespace DPPP

void MSBDAWriter::CreateBDATimeAxis() {
  // Build the table description for BDA_TIME_AXIS.
  TableDesc td(kBDATimeAxisTable, TableDesc::Scratch);
  td.comment() = "Meta information that specify the regularity of the MS.";
  td.addColumn(ScalarColumnDesc<Int>(kTimeAxisId));
  td.addColumn(ScalarColumnDesc<Bool>(kIsBdaApplied));
  td.addColumn(ScalarColumnDesc<Bool>(kSingleFactorPerBL));
  td.addColumn(ScalarColumnDesc<Double>(kMaxTimeInterval));
  td.addColumn(ScalarColumnDesc<Double>(kMinTimeInterval));
  td.addColumn(ScalarColumnDesc<Double>(kUnitTimeInterval));
  td.addColumn(ScalarColumnDesc<Bool>(kIntervalFactors));
  td.addColumn(ScalarColumnDesc<Bool>(kHasBDAOrdering));

  // Add the BDA_TIME_AXIS as a subtable to the output measurementset.
  SetupNewTable new_table(outName_ + '/' + kBDATimeAxisTable, td, Table::New);
  Table bda_time_axis_table(new_table);
  ms_.rwKeywordSet().defineTable(kBDATimeAxisTable, bda_time_axis_table);
}

void MSBDAWriter::CreateBDATimeFactor() {
  // Build the table description for BDA_TIME_FACTOR.
  TableDesc td(kBDATimeFactorTable, TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<Int>(kTimeAxisId));
  td.addColumn(ScalarColumnDesc<Int>(kAntenna1));
  td.addColumn(ScalarColumnDesc<Int>(kAntenna2));
  td.addColumn(ScalarColumnDesc<Double>(kFactor));

  // Add the BDA_TIME_FACTOR as a subtable to the output measurementset.
  SetupNewTable new_table(outName_ + '/' + kBDATimeFactorTable, td, Table::New);
  Table bda_time_factor_table(new_table);
  ms_.rwKeywordSet().defineTable(kBDATimeFactorTable, bda_time_factor_table);
}

void MSBDAWriter::CreateMetaDataFrequencyColumns() {
  Table out_spw = Table(outName_ + '/' + kSpectralWindowTable, Table::Update);

  // Add column BDA_FREQ_AXIS_ID if not exists
  if (!out_spw.tableDesc().isColumn(kBDAFreqAxisId)) {
    ScalarColumnDesc<Int> bdaFreqAxisIdColumn(kBDAFreqAxisId);
    bdaFreqAxisIdColumn.setDefault(-1);
    out_spw.addColumn(bdaFreqAxisIdColumn);
  }

  // Add column BDA_SET_ID if not exists
  if (!out_spw.tableDesc().isColumn(kBDASetId)) {
    ScalarColumnDesc<Int> bdaFreqAxisIdColumn(kBDASetId);
    bdaFreqAxisIdColumn.setDefault(-1);
    out_spw.addColumn(bdaFreqAxisIdColumn);
  }
}

void MSBDAWriter::WriteMetaData() {
  Table bda_time_axis(outName_ + '/' + kBDATimeAxisTable, Table::Update);
  if (bda_time_axis.nrow() > 0) {
    // BDA metadata already exists, do nothing.
    // TODO test that the SPECTRAL_WINDOW has not been overwritten in this case.
    return;
  }

  const Int pid = ObjectID().pid();
  std::size_t min_factor_time = 65535;
  std::size_t max_factor_time = 1;

  WriteTimeFactorRows(pid, min_factor_time, max_factor_time);
  WriteTimeAxisRow(bda_time_axis, pid, min_factor_time, max_factor_time);
  FillSpectralWindowColumns(pid);
}

void MSBDAWriter::WriteTimeFactorRows(const Int& pid, size_t& min_factor_time,
                                      size_t& max_factor_time) {
  Table bda_time_factor(outName_ + '/' + kBDATimeFactorTable, Table::Update);
  int row = bda_time_factor.nrow();
  const Vector<Int>& ant1 = info().getAnt1();
  const Vector<Int>& ant2 = info().getAnt2();
  const std::vector<size_t>& factors = info().getBDAFactors();
  for (std::size_t i = 0; i < info().nbaselines(); ++i) {
    bda_time_factor.addRow();
    const size_t factor = factors[i];
    min_factor_time = std::min(min_factor_time, factor);
    max_factor_time = std::max(max_factor_time, factor);

    ScalarColumn<casacore::Int>(bda_time_factor, kTimeAxisId).put(row, pid);
    ScalarColumn<casacore::Int>(bda_time_factor, kAntenna1).put(row, ant1[i]);
    ScalarColumn<casacore::Int>(bda_time_factor, kAntenna2).put(row, ant2[i]);
    ScalarColumn<casacore::Int>(bda_time_factor, kFactor).put(row, factor);
    ++row;
  }
}

void MSBDAWriter::WriteTimeAxisRow(Table& bda_time_axis, const Int& pid,
                                   const double& min_factor_time,
                                   const double& max_factor_time) {
  const double interval = info().timeInterval();
  int row = bda_time_axis.nrow();
  bda_time_axis.addRow();
  ScalarColumn<casacore::Int>(bda_time_axis, kTimeAxisId).put(row, pid);
  ScalarColumn<casacore::Bool>(bda_time_axis, kIsBdaApplied).put(row, True);
  ScalarColumn<casacore::Bool>(bda_time_axis, kSingleFactorPerBL)
      .put(row, True);
  ScalarColumn<casacore::Double>(bda_time_axis, kMaxTimeInterval)
      .put(row, max_factor_time * interval);
  ScalarColumn<casacore::Double>(bda_time_axis, kMinTimeInterval)
      .put(row, min_factor_time * interval);
  ScalarColumn<casacore::Double>(bda_time_axis, kUnitTimeInterval)
      .put(row, interval);
  ScalarColumn<casacore::Bool>(bda_time_axis, kIntervalFactors).put(row, True);
  ScalarColumn<casacore::Bool>(bda_time_axis, kHasBDAOrdering).put(row, True);
}

void MSBDAWriter::FillSpectralWindowColumns(const Int& pid) {
  Table spectral_window(outName_ + '/' + kSpectralWindowTable, Table::Update);
  ScalarColumn<casacore::Int>(spectral_window, kBDAFreqAxisId).fillColumn(pid);
  ScalarColumn<casacore::Int>(spectral_window, kBDASetId).fillColumn(0);
}

}  // namespace DPPP
}  // namespace DP3
