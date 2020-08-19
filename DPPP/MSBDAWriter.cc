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
#include <map>

#include "../Common/ParameterSet.h"
#include "BDABuffer.h"
#include "MSBDAWriter.h"
#include "MSWriter.h"
#include "DPLogger.h"

#include <casacore/tables/TaQL/TableParse.h>

using casacore::Array;
using casacore::ArrayColumn;
using casacore::ArrayColumnDesc;
using casacore::Block;
using casacore::Bool;
using casacore::ColumnDesc;
using casacore::Complex;
using casacore::Double;
using casacore::Float;
using casacore::IncrementalStMan;
using casacore::Int;
using casacore::IPosition;
using casacore::MeasurementSet;
using casacore::MS;
using casacore::ScalarColumn;
using casacore::ScalarColumnDesc;
using casacore::SetupNewTable;
using casacore::StandardStMan;
using casacore::String;
using casacore::Table;
using casacore::TableCopy;
using casacore::TableDesc;
using casacore::TableLock;
using casacore::TiledColumnStMan;
using casacore::True;
using casacore::Vector;

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
const std::string kDataDescTable = "DATA_DESCRIPTION";

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
  if (info_in.ntimeAvgs().size() != info_in.nbaselines()) {
    throw std::invalid_argument("Invalid time averaging factors");
  }

  DPStep::updateInfo(info_in);
  CreateMS();

  WriteMetaData();

  MSWriter::writeHistory(ms_, parset_);
  ms_.flush(true, true);
  DPLOG_INFO("Finished preparing output MS", false);
}

bool MSBDAWriter::process(std::unique_ptr<BDABuffer> buffer) {
  buffer->SetBaseRowNr(ms_.nrow());

  const std::vector<BDABuffer::Row>& rows = buffer->GetRows();

  ms_.addRow(rows.size());

  Vector<Float> sigma_weight(info().ncorr(), 1);

  ScalarColumn<Double> time(ms_, MS::columnName(MS::TIME));
  ScalarColumn<Double> time_centroid(ms_, MS::columnName(MS::TIME_CENTROID));
  ScalarColumn<Double> exposure(ms_, MS::columnName(MS::EXPOSURE));
  ScalarColumn<Int> ant1(ms_, MS::columnName(MS::ANTENNA1));
  ScalarColumn<Int> ant2(ms_, MS::columnName(MS::ANTENNA2));
  ArrayColumn<Complex> data(ms_, MS::columnName(MS::DATA));
  ArrayColumn<Float> weights(ms_, MS::columnName(MS::WEIGHT_SPECTRUM));
  ArrayColumn<Bool> flags(ms_, MS::columnName(MS::FLAG));
  ScalarColumn<Bool> flags_row(ms_, MS::columnName(MS::FLAG_ROW));
  ArrayColumn<Double> uvw(ms_, MS::columnName(MS::UVW));
  ScalarColumn<Double> interval(ms_, MS::columnName(MS::INTERVAL));
  ArrayColumn<Float> sigma(ms_, MS::columnName(MS::SIGMA));
  ArrayColumn<Float> weight(ms_, MS::columnName(MS::WEIGHT));
  ScalarColumn<Int> dataDescId(ms_, MS::columnName(MS::DATA_DESC_ID));

  std::vector<DP3::rownr_t> row_nrs;
  row_nrs.reserve(rows.size());
  for (const BDABuffer::Row& row : rows) {
    time.put(row.row_nr, row.time);
    time_centroid.put(row.row_nr, row.time);
    interval.put(row.row_nr, row.interval);
    exposure.put(row.row_nr, row.exposure);

    ant1.put(row.row_nr, info().getAnt1()[row.baseline_nr]);
    ant2.put(row.row_nr, info().getAnt2()[row.baseline_nr]);

    const std::size_t n_chan = info().chanFreqs(row.baseline_nr).size();
    const IPosition dim(2, info().ncorr(), n_chan);
    data.put(row.row_nr, Array<Complex>(dim, row.data, casacore::SHARE));
    weights.put(row.row_nr, Array<Float>(dim, row.weights, casacore::SHARE));
    flags.put(row.row_nr, Array<Bool>(dim, row.flags, casacore::SHARE));
    // Set the row_flag if all flags in the row are true / none are false.
    const bool row_flag =
        std::count(row.flags, row.flags + row.GetDataSize(), false) == 0;
    flags_row.put(row.row_nr, row_flag);

    const IPosition uvw_dim(1, 3);
    uvw.put(row.row_nr, Array<Double>(uvw_dim, row.uvw));

    // Fill values in all the cells of various columns.
    sigma.put(row.row_nr, sigma_weight);
    weight.put(row.row_nr, sigma_weight);
    dataDescId.put(row.row_nr, nchanToDescId[n_chan]);

    row_nrs.push_back(row.row_nr);
  }

  // Create a table view containing only the added rows.
  Table tbl_added(ms_(row_nrs));

  return true;
}

void MSBDAWriter::finish() {}

void MSBDAWriter::show(std::ostream& os) const {}

void MSBDAWriter::CreateMS() {
  CreateMainTable();

  CreateMetaDataFrequencyColumns();
  CreateBDATimeAxis();
  CreateBDATimeFactor();
}

/// Based on the example on
/// https://casacore.github.io/casacore/group__Tables__module.html#Tables:performance
void MSBDAWriter::CreateMainTable() {
  // Create an empty table.
  TableDesc td = MS::requiredTableDesc();
  MS::addColumnToDesc(td, MS::DATA);
  MS::addColumnToDesc(td, MS::WEIGHT_SPECTRUM);

  IncrementalStMan man_inc;
  StandardStMan man_std(32768);

  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable setup(outName_, td, opt);

  // Set fixed shapes for columns that allow it.
  setup.setShapeColumn(MS::columnName(MS::UVW), IPosition(1, 3));

  // Set fixed fields with the incremental storage manager.
  setup.bindAll(man_inc);
  // Set fixed-shape fields with the standard storage manager.
  setup.bindColumn(MS::columnName(MS::TIME), man_std);
  setup.bindColumn(MS::columnName(MS::ANTENNA1), man_std);
  setup.bindColumn(MS::columnName(MS::ANTENNA2), man_std);
  setup.bindColumn(MS::columnName(MS::EXPOSURE), man_std);
  setup.bindColumn(MS::columnName(MS::TIME_CENTROID), man_std);
  setup.bindColumn(MS::columnName(MS::UVW), man_std);
  setup.bindColumn(MS::columnName(MS::DATA), man_std);
  setup.bindColumn(MS::columnName(MS::WEIGHT_SPECTRUM), man_std);
  setup.bindColumn(MS::columnName(MS::FLAG), man_std);
  setup.bindColumn(MS::columnName(MS::FLAG_CATEGORY), man_std);
  setup.bindColumn(MS::columnName(MS::FLAG_ROW), man_std);

  ms_ = Table(setup);

  // Create empty subtables.
  MeasurementSet(ms_).createDefaultSubtables(Table::New);
  if (reader_) {
    DPLOG_INFO("Copying info and subtables ...", false);
    // Copy the info and subtables.
    TableCopy::copyInfo(ms_, reader_->table());

    Block<String> omitted_subtables(4);
    omitted_subtables[0] = kBDATimeAxisTable;
    omitted_subtables[1] = kBDATimeFactorTable;
    omitted_subtables[2] = kSpectralWindowTable;
    omitted_subtables[3] = kDataDescTable;
    TableCopy::copySubTables(ms_, reader_->table(), false, omitted_subtables);
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
  td.addColumn(ScalarColumnDesc<Int>(MS::columnName(MS::ANTENNA1)));
  td.addColumn(ScalarColumnDesc<Int>(MS::columnName(MS::ANTENNA2)));
  td.addColumn(ScalarColumnDesc<Int>(kFactor));

  // Add the BDA_TIME_FACTOR as a subtable to the output measurementset.
  SetupNewTable new_table(outName_ + '/' + kBDATimeFactorTable, td, Table::New);
  Table bda_time_factor_table(new_table);
  ms_.rwKeywordSet().defineTable(kBDATimeFactorTable, bda_time_factor_table);
}

void MSBDAWriter::CreateMetaDataFrequencyColumns() {
  Table out_spw = Table(outName_ + '/' + kSpectralWindowTable, Table::Update);

  // Add column BDA_FREQ_AXIS_ID
  ScalarColumnDesc<Int> bdaFreqAxisIdColumn(kBDAFreqAxisId);
  bdaFreqAxisIdColumn.setDefault(-1);
  out_spw.addColumn(bdaFreqAxisIdColumn);

  // Add column BDA_SET_ID
  ScalarColumnDesc<Int> bdaSetIdColumn(kBDASetId);
  bdaSetIdColumn.setDefault(-1);
  out_spw.addColumn(bdaSetIdColumn);

  // Remove fixed size options from columns
  TableDesc td = out_spw.tableDesc();
  td.rwColumnDesc("CHAN_FREQ").setOptions(ColumnDesc::Direct);
  td.rwColumnDesc("CHAN_WIDTH").setOptions(ColumnDesc::Direct);
  td.rwColumnDesc("EFFECTIVE_BW").setOptions(ColumnDesc::Direct);
  td.rwColumnDesc("RESOLUTION").setOptions(ColumnDesc::Direct);
}

void MSBDAWriter::WriteMetaData() {
  const Int pid = 0;
  unsigned int min_factor_time = 65535;
  unsigned int max_factor_time = 1;

  WriteTimeFactorRows(pid, min_factor_time, max_factor_time);
  WriteTimeAxisRow(pid, min_factor_time, max_factor_time);

  OverwriteSubTables(pid);
}

void MSBDAWriter::WriteTimeFactorRows(const Int& pid,
                                      unsigned int& min_factor_time,
                                      unsigned int& max_factor_time) {
  Table bda_time_factor(outName_ + '/' + kBDATimeFactorTable, Table::Update);
  int row = bda_time_factor.nrow();
  const Vector<Int>& ant1 = info().getAnt1();
  const Vector<Int>& ant2 = info().getAnt2();
  for (std::size_t i = 0; i < info().nbaselines(); ++i) {
    bda_time_factor.addRow();
    const unsigned int factor = info().ntimeAvg(i);
    min_factor_time = std::min(min_factor_time, factor);
    max_factor_time = std::max(max_factor_time, factor);

    ScalarColumn<Int>(bda_time_factor, kTimeAxisId).put(row, pid);
    ScalarColumn<Int>(bda_time_factor, MS::columnName(MS::ANTENNA1))
        .put(row, ant1[i]);
    ScalarColumn<Int>(bda_time_factor, MS::columnName(MS::ANTENNA2))
        .put(row, ant2[i]);
    ScalarColumn<Int>(bda_time_factor, kFactor).put(row, factor);
    ++row;
  }
}

void MSBDAWriter::WriteTimeAxisRow(const Int& pid,
                                   const unsigned int& min_factor_time,
                                   const unsigned int& max_factor_time) {
  Table bda_time_axis(outName_ + '/' + kBDATimeAxisTable, Table::Update);
  const double interval = info().timeInterval();
  int row = bda_time_axis.nrow();
  bda_time_axis.addRow();
  ScalarColumn<Int>(bda_time_axis, kTimeAxisId).put(row, pid);
  ScalarColumn<Bool>(bda_time_axis, kIsBdaApplied).put(row, True);
  ScalarColumn<Bool>(bda_time_axis, kSingleFactorPerBL).put(row, True);
  ScalarColumn<Double>(bda_time_axis, kMaxTimeInterval)
      .put(row, max_factor_time * interval);
  ScalarColumn<Double>(bda_time_axis, kMinTimeInterval)
      .put(row, min_factor_time * interval);
  ScalarColumn<Double>(bda_time_axis, kUnitTimeInterval).put(row, interval);
  ScalarColumn<Bool>(bda_time_axis, kIntervalFactors).put(row, True);
  ScalarColumn<Bool>(bda_time_axis, kHasBDAOrdering).put(row, True);
}

void MSBDAWriter::OverwriteSubTables(const Int& pid) {
  unsigned int id = 0;
  std::string name;
  int measFreqRef;

  Table outDD(outName_ + '/' + kDataDescTable, Table::Update);
  Table outSPW(outName_ + '/' + kSpectralWindowTable, Table::Update);

  if (outSPW.nrow() != outDD.nrow())
    throw std::runtime_error(
        "nrow in SPECTRAL_WINDOW table is not the same as nrow in "
        "DATA_DESCRIPTION table");

  // Remove all rows before and after the selected band.
  // Do it from the end, otherwise row numbers change.
  for (unsigned int i = outSPW.nrow(); i > 0;) {
    if (--i == reader_->spectralWindow()) {
      measFreqRef = outSPW.col("MEAS_FREQ_REF").getInt(i);
      name = outSPW.col("NAME").getString(i);
    }
    outSPW.removeRow(i);
    outDD.removeRow(i);
  }

  for (std::size_t i = 0; i < info().nbaselines(); ++i) {
    std::size_t nchanFreqs = info().chanFreqs(i).size();
    if (nchanToDescId.count(nchanFreqs)) {
      continue;
    }

    outDD.addRow();
    ScalarColumn<Int>(outDD, "SPECTRAL_WINDOW_ID").put(id, id);

    outSPW.addRow();
    ScalarColumn<Int>(outSPW, "NUM_CHAN").put(id, nchanFreqs);
    ScalarColumn<Int>(outSPW, "MEAS_FREQ_REF").put(id, measFreqRef);
    ScalarColumn<String>(outSPW, "NAME").put(id, name);
    ArrayColumn<Double>(outSPW, "CHAN_FREQ")
        .put(id, Vector<double>(info().chanFreqs(i)));
    ArrayColumn<Double>(outSPW, "CHAN_WIDTH")
        .put(id, Vector<double>(info().chanWidths(i)));
    ArrayColumn<Double>(outSPW, "EFFECTIVE_BW")
        .put(id, Vector<double>(info().effectiveBW(i)));
    ArrayColumn<Double>(outSPW, "RESOLUTION")
        .put(id, Vector<double>(info().resolutions(i)));
    ScalarColumn<Double>(outSPW, "TOTAL_BANDWIDTH").put(id, info().totalBW());
    ScalarColumn<Double>(outSPW, "REF_FREQUENCY").put(id, info().refFreq());
    ScalarColumn<Int>(outSPW, kBDAFreqAxisId).put(id, pid);
    // When we support multiple SPECTRAL_WINDOW entries this value will be used
    // to keep track of the old windows
    ScalarColumn<Int>(outSPW, kBDASetId).put(id, 0);

    nchanToDescId[nchanFreqs] = id;
    ++id;
  }
}

}  // namespace DPPP
}  // namespace DP3
