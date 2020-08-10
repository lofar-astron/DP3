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
using casacore::Bool;
using casacore::Double;
using casacore::Int;
using casacore::MeasurementSet;
using casacore::ScalarColumn;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::StandardStMan;
using casacore::Table;
using casacore::TableCopy;
using casacore::TableDesc;
using casacore::TableLock;

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
const std::string kFlagCategory = "FLAG_CATEGORY";
const std::string kFlagRow = "FLAG_ROW";
const std::string kBdaFreqAxisId = "BDA_FREQ_AXIS_ID";
/// @}
}  // namespace

namespace DP3 {
namespace DPPP {

MSBDAWriter::MSBDAWriter(MSReader* reader, const std::string& outName,
                         const ParameterSet& parset, const std::string& prefix)
    : reader_(reader), outName_(outName), parset_(parset), prefix_(prefix) {
  overwrite_ = parset.getBool(prefix + "overwrite", false);
}

MSBDAWriter::~MSBDAWriter() {}

void MSBDAWriter::updateInfo(const DPInfo& infoIn) {
  DPStep::updateInfo(infoIn);
  CreateMS();

  MSWriter::writeHistory(ms_, parset_);
  ms_.flush(true, true);
}

bool MSBDAWriter::process(std::unique_ptr<BDABuffer> buffer) {
  // TODO: buffer->SetBaseRowNr(ms_.nrow());

  const std::vector<BDABuffer::Row>& rows = buffer->GetRows();

  ms_.addRow(rows.size());

  ScalarColumn<casacore::Int> ant1(ms_, kAntenna1);
  ScalarColumn<casacore::Int> ant2(ms_, kAntenna2);
  ArrayColumn<casacore::Complex> data(ms_, kData);
  ArrayColumn<casacore::Float> weights(ms_, kWeightSpectrum);
  ArrayColumn<casacore::Bool> flags(ms_, kFlag);
  std::vector<DP3::rownr_t> row_nrs;
  row_nrs.reserve(rows.size());
  for (const BDABuffer::Row& row : rows) {
    ant1.put(row.row_nr_, info().getAnt1()[row.baseline_nr_]);
    ant2.put(row.row_nr_, info().getAnt2()[row.baseline_nr_]);

    // TODO: When DPInfo is updated, remove the assert and the comment below.
    assert(row.baseline_nr_ == 0);
    const std::size_t n_chan = info().chanFreqs(/*row.baseline_nr_*/).size();
    const casacore::IPosition dim(2, info().ncorr(), n_chan);
    data.put(row.row_nr_, casacore::Array<casacore::Complex>(dim, row.data_,
                                                             casacore::SHARE));
    weights.put(row.row_nr_, casacore::Array<casacore::Float>(dim, row.weights_,
                                                              casacore::SHARE));
    flags.put(row.row_nr_, casacore::Array<casacore::Bool>(dim, row.flags_,
                                                           casacore::SHARE));

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
  bool create_bda_time_axis = true;

  CreateMainTable();
  // TODO fill the main table

  AddMetaDataFrequency();

  if (create_bda_time_axis) {
    CreateBDATimeAxis();
  }
}

void MSBDAWriter::CreateMainTable() {
  // Build the table description.
  casacore::Block<casacore::String> fixedColumns(16);
  fixedColumns[0] = kTime;
  fixedColumns[1] = kAntenna1;
  fixedColumns[2] = kAntenna2;
  fixedColumns[3] = kFeed1;
  fixedColumns[4] = kFeed2;
  fixedColumns[5] = kDataDescId;
  fixedColumns[6] = kProcessorId;
  fixedColumns[7] = kFieldId;
  fixedColumns[8] = kInterval;
  fixedColumns[9] = kExposure;
  fixedColumns[10] = kTimeCentroid;
  fixedColumns[11] = kScanNumber;
  fixedColumns[12] = kArrayId;
  fixedColumns[13] = kObservationId;
  fixedColumns[14] = kStateId;
  fixedColumns[15] = kUVW;
  fixedColumns[16] = kSigma;
  fixedColumns[17] = kWeight;
  // fixedColumns[1] = "FLAG_CATEGORY";
  // fixedColumns[12] = "FLAG_ROW";
  // fixedColumns[20] = "BDA_FREQ_AXIS_ID";
  const TableDesc& td = reader_->table().project(fixedColumns).tableDesc();

  // Add DATA, WEIGHT_SPECTRUM and FLAG columns.

  // Setup a new table from the description
  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable table(outName_, td, opt);
  table.bindAll(StandardStMan(32768));

  // Create the table
  ms_ = Table(table);

  // Copy the info and subtables.
  TableCopy::copyInfo(ms_, reader_->table());

  casacore::Block<casacore::String> omitted_subtables(1);
  omitted_subtables[0] = "BDA_TIME_AXIS";
  TableCopy::copySubTables(ms_, reader_->table(), false, omitted_subtables);
}

void MSBDAWriter::CreateBDATimeAxis() {
  const std::string tableName = "BDA_TIME_AXIS";
  const std::string version_bda_time_axis = "1.0";

  // Build the table description for BDA_TIME_AXIS.
  TableDesc td(tableName, version_bda_time_axis, TableDesc::Scratch);
  td.comment() = "Meta information that specify the regularity of the MS.";
  td.addColumn(ScalarColumnDesc<Int>("BDA_TIME_AXIS_ID", "unique id"));
  td.addColumn(ScalarColumnDesc<Bool>("IS_BDA_APPLIED",
                                      "BDA has been applied to the time axis"));
  td.addColumn(ScalarColumnDesc<Bool>(
      "SINGLE_FACTOR_PER_BASELINE",
      "If for every baseline, the averaging factor is constant in time."));
  td.addColumn(ScalarColumnDesc<Double>(
      "MAX_TIME_INTERVAL", "maximum TIME_INTERVAL over this subset"));
  td.addColumn(ScalarColumnDesc<Double>(
      "MIN_TIME_INTERVAL", "minimum TIME_INTERVAL over this subset"));
  td.addColumn(ScalarColumnDesc<Double>(
      "UNIT_TIME_INTERVAL",
      "An integer multiple of the UNIT_TIME_INTERVAL value"));
  td.addColumn(ScalarColumnDesc<Double>(
      "INTEGER_INTERVAL_FACTORS",
      "The TIME_INTERVAL between two consecutive timesteps"));
  td.addColumn(ScalarColumnDesc<Bool>(
      "HAS_BDA_ORDERING",
      "if a row starts at T_0 (where T_0 = TIME - 0.5 * TIME_INTERVAL) then "
      "all visibilities that end before T_0 are before this row"));

  // Add the BDA_TIME_AXIS as a subtable to the output measurementset.
  SetupNewTable newTable(outName_ + '/' + tableName, td, Table::New);
  Table bdaTimeAxisTable(newTable);
  ms_.rwKeywordSet().defineTable(tableName, bdaTimeAxisTable);
}

void MSBDAWriter::CreateBDATimeFactor() {
  const std::string tableName = "BDA_TIME_FACTOR";
  const std::string version_bda_time_axis = "1.0";

  // Build the table description for BDA_TIME_AXIS.
  TableDesc td(tableName, version_bda_time_axis, TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<Int>(
      "BDA_TIME_AXIS_ID", "refers back to a row in the table BDA_TIME_AXIS"));
  td.addColumn(ScalarColumnDesc<Int>("ANTENNA1"));
  td.addColumn(ScalarColumnDesc<Int>("ANTENNA2"));
  td.addColumn(ScalarColumnDesc<Double>(
      "FACTOR",
      "Averaging factor for this baseline relative to UNIT_TIME_INTERVAL in "
      "the table TIME_AXIS_BDA"));

  // Add the BDA_TIME_FACTOR as a subtable to the output measurementset.
  SetupNewTable newTable(outName_ + '/' + tableName, td, Table::New);
  Table bdaTimeFactorTable(newTable);
  ms_.rwKeywordSet().defineTable(tableName, bdaTimeFactorTable);
}

void MSBDAWriter::AddMetaDataFrequency() {
  Table outSPW = Table(outName_ + "/SPECTRAL_WINDOW", Table::Update);

  // Add table if not exists
  if (outSPW.tableDesc().isColumn("BDA_FREQ_AXIS_ID")) {
    ScalarColumnDesc<Int> bdaFreqAxisIdColumn("BDA_FREQ_AXIS_ID");
    bdaFreqAxisIdColumn.setDefault(-1);
    outSPW.addColumn(bdaFreqAxisIdColumn);
  }

  // TODO fill / overwrite column
}

}  // namespace DPPP
}  // namespace DP3
