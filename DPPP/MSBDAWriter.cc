// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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

#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/DataMan/StandardStMan.h>

#include "../Common/ParameterSet.h"
#include "MSBDAWriter.h"

using casacore::Block;
using casacore::Bool;
using casacore::Double;
using casacore::Int;
using casacore::MeasurementSet;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::StandardStMan;
using casacore::String;
using casacore::Table;
using casacore::TableDesc;
using casacore::TableLock;

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
  TableDesc td(reader_->table().tableDesc());
  // Block<String> fixedColumns(20);
  // fixedColumns[0] = "UVW";
  // fixedColumns[1] = "FLAG_CATEGORY";
  // fixedColumns[2] = "WEIGHT";
  // fixedColumns[3] = "SIGMA";
  // fixedColumns[4] = "ANTENNA1";
  // fixedColumns[5] = "ANTENNA2";
  // fixedColumns[6] = "ARRAY_ID";
  // fixedColumns[7] = "DATA_DESC_ID";
  // fixedColumns[8] = "EXPOSURE";
  // fixedColumns[9] = "FEED1";
  // fixedColumns[10] = "FEED2";
  // fixedColumns[11] = "FIELD_ID";
  // fixedColumns[12] = "FLAG_ROW";
  // fixedColumns[13] = "INTERVAL";
  // fixedColumns[14] = "OBSERVATION_ID";
  // fixedColumns[15] = "PROCESSOR_ID";
  // fixedColumns[16] = "SCAN_NUMBER";
  // fixedColumns[17] = "STATE_ID";
  // fixedColumns[18] = "TIME";
  // fixedColumns[19] = "TIME_CENTROID";
  // // TODO, check if here or SPECTRAL_WINDOW table
  // fixedColumns[20] = "BDA_FREQ_AXIS_ID";
  // Table tempTable = reader_->table().project(fixedColumns);
  // TableDesc newdesc = tempTable.tableDesc();

  // TODO monday: this gives a table that is not valid

  // Setup a new table from the description
  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable newtab(outName_, td, opt);

  // Create storage managers
  StandardStMan stmanStand;
  newtab.bindAll(stmanStand);

  // Create the table
  ms_ = Table(newtab);
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
