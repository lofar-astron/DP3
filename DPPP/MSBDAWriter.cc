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

using casacore::MeasurementSet;
using casacore::Table;
using casacore::TableLock;
using casacore::Block;
using casacore::String;
using casacore::Int;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::TableDesc;
using casacore::StandardStMan;

namespace DP3 {
namespace DPPP {

MSBDAWriter::MSBDAWriter(MSReader* reader, const std::string& outName,
                         const ParameterSet& parset, const std::string& prefix)
    : reader_(reader), outName_(outName), parset_(parset), prefix_(prefix) {
      overwrite_ = parset.getBool(prefix + "overwrite", false);
    }

MSBDAWriter::~MSBDAWriter() {}

void MSBDAWriter::updateInfo(const DPInfo& infoIn) { createMS(infoIn); }

void MSBDAWriter::finish() {}

void MSBDAWriter::show(std::ostream& os) const {}

void MSBDAWriter::createMS(const DPInfo& info) {
  bool create_bda_time_axis = true;

  createMainTable(info);
  // TODO fill the main table

  if (create_bda_time_axis) {
    createBDATimeAxis(info);
  }
}

void MSBDAWriter::createMainTable(const DPInfo& info) {
  // Build the table description.
  TableDesc td(reader_->table().tableDesc());

  // Setup a new table from the description
  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable newtab(outName_, td, opt);

  // Create storage managers
  StandardStMan stmanStand;
  newtab.bindAll (stmanStand);

  // Create the table
  ms_ = Table(newtab);
}

void MSBDAWriter::createBDATimeAxis(const DPInfo& info) {
  const std::string tableName = "BDA_TIME_AXIS";
  const std::string version_bda_time_axis = "1.0";

  // Build the table description for BDA_TIME_AXIS.
  TableDesc td(tableName, version_bda_time_axis, TableDesc::Scratch);
  td.comment() = "This subtable specifies regularity only for the time axis.";
  td.addColumn(ScalarColumnDesc<Int>("BDA_TIME_AXIS_ID", "unique id"));
  td.addColumn(ScalarColumnDesc<bool>("IS_BDA_APPLIED", "unique id"));

  // Add the BDA_TIME_AXIS as a subtable to the output measurementset.
  SetupNewTable newTable(outName_ + '/' + tableName, td, Table::New);
  Table bdaTimeAxisTable(newTable);
  ms_.rwKeywordSet().defineTable(tableName, bdaTimeAxisTable);
}

}  // namespace DPPP
}  // namespace DP3
