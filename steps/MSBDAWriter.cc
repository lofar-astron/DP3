// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MSBDAWriter.h"
#include "MSWriter.h"

#include "../common/ParameterSet.h"
#include "../base/BDABuffer.h"
#include "../base/DPLogger.h"
#include "../base/MS.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/DataMan/IncrementalStMan.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableCopy.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <map>

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
using casacore::Table;
using casacore::TableCopy;
using casacore::TableDesc;
using casacore::TableLock;
using casacore::True;
using casacore::Vector;

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace {
// Keywords
const std::string kBDATimeAxisVersionKW = "BDA_TIME_AXIS_VERSION";
const std::string kBDATimeAxisVersion = "1.0";
/// @}
}  // namespace

namespace dp3 {
namespace steps {

MSBDAWriter::MSBDAWriter(InputStep* reader, const std::string& out_name,
                         const common::ParameterSet& parset,
                         const std::string& prefix)
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

  Step::updateInfo(info_in);
  nbl_ = info().nbaselines();
  ncorr_ = info().ncorr();
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

  std::vector<common::rownr_t> row_nrs;
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

void MSBDAWriter::addToMS(const std::string&) {
  getPrevStep()->addToMS(outName_);
}

void MSBDAWriter::show(std::ostream& os) const {
  os << "MSBDAWriter " << prefix_ << '\n';
  os << "  output MS:      " << ms_.tableName() << '\n';
  os << "  ncorrelations:  " << ncorr_ << '\n';
  os << "  nbaselines:     " << nbl_ << '\n';
  os << "  DATA column:    DATA" << '\n';
  os << "  Compressed:     no\n";
}

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

    Block<casacore::String> omitted_subtables(4);
    omitted_subtables[0] = base::DP3MS::kBDATimeAxisTable;
    omitted_subtables[1] = base::DP3MS::kBDAFactorsTable;
    omitted_subtables[2] = base::DP3MS::kSpectralWindowTable;
    omitted_subtables[3] = base::DP3MS::kDataDescTable;
    TableCopy::copySubTables(ms_, reader_->table(), false, omitted_subtables);
  }
}  // namespace DPPP

void MSBDAWriter::CreateBDATimeAxis() {
  // Build the table description for BDA_TIME_AXIS.
  TableDesc td(base::DP3MS::kBDATimeAxisTable, TableDesc::Scratch);
  td.comment() = "Meta information that specify the regularity of the MS.";
  td.rwKeywordSet().define(kBDATimeAxisVersionKW, kBDATimeAxisVersion);
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kTimeAxisId));
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kFieldId));
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kBDAFreqAxisId));
  td.addColumn(ScalarColumnDesc<Bool>(base::DP3MS::kIsBdaApplied));
  td.addColumn(ScalarColumnDesc<Bool>(base::DP3MS::kSingleFactorPerBL));
  td.addColumn(ScalarColumnDesc<Double>(base::DP3MS::kMaxTimeInterval));
  td.addColumn(ScalarColumnDesc<Double>(base::DP3MS::kMinTimeInterval));
  td.addColumn(ScalarColumnDesc<Double>(base::DP3MS::kUnitTimeInterval));
  td.addColumn(ScalarColumnDesc<Bool>(base::DP3MS::kIntervalFactors));
  td.addColumn(ScalarColumnDesc<Bool>(base::DP3MS::kHasBDAOrdering));

  // Add the BDA_TIME_AXIS as a subtable to the output measurementset.
  SetupNewTable new_table(outName_ + '/' + base::DP3MS::kBDATimeAxisTable, td,
                          Table::New);
  Table bda_time_axis_table(new_table);
  ms_.rwKeywordSet().defineTable(base::DP3MS::kBDATimeAxisTable,
                                 bda_time_axis_table);
}

void MSBDAWriter::CreateBDATimeFactor() {
  // Build the table description for BDA_FACTORS.
  TableDesc td(base::DP3MS::kBDAFactorsTable, TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kTimeAxisId));
  td.addColumn(ScalarColumnDesc<Int>(MS::columnName(MS::ANTENNA1)));
  td.addColumn(ScalarColumnDesc<Int>(MS::columnName(MS::ANTENNA2)));
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kFactor));
  td.addColumn(ScalarColumnDesc<Int>(base::DP3MS::kSpectralWindowId));

  // Add the BDA_FACTORS as a subtable to the output measurementset.
  SetupNewTable new_table(outName_ + '/' + base::DP3MS::kBDAFactorsTable, td,
                          Table::New);
  Table bda_time_factor_table(new_table);
  ms_.rwKeywordSet().defineTable(base::DP3MS::kBDAFactorsTable,
                                 bda_time_factor_table);
}

void MSBDAWriter::CreateMetaDataFrequencyColumns() {
  using MS_SPW = casacore::MSSpectralWindow;

  Table out_spw =
      Table(outName_ + '/' + base::DP3MS::kSpectralWindowTable, Table::Update);

  // Don't add BDA_FREQ_AXIS_ID

  // Add column BDA_SET_ID
  ScalarColumnDesc<Int> bdaSetIdColumn(base::DP3MS::kBDASetId);
  // When we support multiple SPECTRAL_WINDOW entries this value will be used
  // to keep track of the old windows
  bdaSetIdColumn.setDefault(0);
  out_spw.addColumn(bdaSetIdColumn);

  // Remove fixed size options from columns
  TableDesc td = out_spw.tableDesc();
  td.rwColumnDesc(MS_SPW::columnName(MS_SPW::CHAN_FREQ))
      .setOptions(ColumnDesc::Direct);
  td.rwColumnDesc(MS_SPW::columnName(MS_SPW::CHAN_WIDTH))
      .setOptions(ColumnDesc::Direct);
  td.rwColumnDesc(MS_SPW::columnName(MS_SPW::EFFECTIVE_BW))
      .setOptions(ColumnDesc::Direct);
  td.rwColumnDesc(MS_SPW::columnName(MS_SPW::RESOLUTION))
      .setOptions(ColumnDesc::Direct);
}

void MSBDAWriter::WriteMetaData() {
  const Int bda_set_id = 0;
  unsigned int min_factor_time = 65535;
  unsigned int max_factor_time = 1;

  OverwriteSubTables(bda_set_id);
  WriteTimeFactorRows(bda_set_id, min_factor_time, max_factor_time);
  WriteTimeAxisRow(bda_set_id, min_factor_time, max_factor_time);
}

void MSBDAWriter::WriteTimeFactorRows(Int bda_set_id,
                                      unsigned int& min_factor_time,
                                      unsigned int& max_factor_time) {
  Table bda_time_factor(outName_ + '/' + base::DP3MS::kBDAFactorsTable,
                        Table::Update);

  ScalarColumn<Int> col_time_axis_id(bda_time_factor, base::DP3MS::kTimeAxisId);
  ScalarColumn<Int> col_ant1(bda_time_factor, MS::columnName(MS::ANTENNA1));
  ScalarColumn<Int> col_ant2(bda_time_factor, MS::columnName(MS::ANTENNA2));
  ScalarColumn<Int> col_factor(bda_time_factor, base::DP3MS::kFactor);
  ScalarColumn<Int> col_spw_id(bda_time_factor, base::DP3MS::kSpectralWindowId);

  const Vector<Int>& ant1 = info().getAnt1();
  const Vector<Int>& ant2 = info().getAnt2();
  for (std::size_t i = 0; i < info().nbaselines(); ++i) {
    std::size_t nchan = info().chanFreqs(i).size();
    const int row = bda_time_factor.nrow();
    bda_time_factor.addRow();
    const unsigned int factor = info().ntimeAvg(i);
    min_factor_time = std::min(min_factor_time, factor);
    max_factor_time = std::max(max_factor_time, factor);

    col_time_axis_id.put(row, bda_set_id);
    col_ant1.put(row, ant1[i]);
    col_ant2.put(row, ant2[i]);
    col_factor.put(row, factor);
    col_spw_id.put(row, nchanToDescId[nchan]);
  }
}

void MSBDAWriter::WriteTimeAxisRow(Int bda_set_id, unsigned int min_factor_time,
                                   unsigned int max_factor_time) {
  Table bda_time_axis(outName_ + '/' + base::DP3MS::kBDATimeAxisTable,
                      Table::Update);
  const double interval = info().timeInterval();
  int row = bda_time_axis.nrow();
  bda_time_axis.addRow();
  ScalarColumn<Int>(bda_time_axis, base::DP3MS::kTimeAxisId)
      .put(row, bda_set_id);
  ScalarColumn<Bool>(bda_time_axis, base::DP3MS::kIsBdaApplied).put(row, True);
  ScalarColumn<Bool>(bda_time_axis, base::DP3MS::kSingleFactorPerBL)
      .put(row, True);
  ScalarColumn<Double>(bda_time_axis, base::DP3MS::kMaxTimeInterval)
      .put(row, max_factor_time * interval);
  ScalarColumn<Double>(bda_time_axis, base::DP3MS::kMinTimeInterval)
      .put(row, min_factor_time * interval);
  ScalarColumn<Double>(bda_time_axis, base::DP3MS::kUnitTimeInterval)
      .put(row, interval);
  ScalarColumn<Bool>(bda_time_axis, base::DP3MS::kIntervalFactors)
      .put(row, True);
  ScalarColumn<Bool>(bda_time_axis, base::DP3MS::kHasBDAOrdering)
      .put(row, True);

  ScalarColumn<Int>(bda_time_axis, base::DP3MS::kFieldId).put(row, -1);
  ScalarColumn<Int>(bda_time_axis, base::DP3MS::kBDAFreqAxisId).put(row, -1);
}

void MSBDAWriter::OverwriteSubTables(const Int bda_set_id) {
  using MS_DD = casacore::MSDataDescription;
  using MS_SPW = casacore::MSSpectralWindow;

  std::string name;
  int measFreqRef;

  Table outDD(outName_ + '/' + base::DP3MS::kDataDescTable, Table::Update);
  Table outSPW(outName_ + '/' + base::DP3MS::kSpectralWindowTable,
               Table::Update);

  if (outSPW.nrow() != outDD.nrow())
    throw std::runtime_error(
        "nrow in SPECTRAL_WINDOW table is not the same as nrow in "
        "DATA_DESCRIPTION table");

  // Remove all rows before and after the selected band.
  // Do it from the end, otherwise row numbers change.
  for (unsigned int i = outSPW.nrow(); i > 0;) {
    --i;
    if (i == reader_->spectralWindow()) {
      measFreqRef =
          outSPW.col(MS_SPW::columnName(MS_SPW::MEAS_FREQ_REF)).getInt(i);
      name = outSPW.col(MS_SPW::columnName(MS_SPW::NAME)).getString(i);
    }
    outSPW.removeRow(i);
    outDD.removeRow(i);
  }

  ScalarColumn<Int> col_spw_id(outDD,
                               MS_DD::columnName(MS_DD::SPECTRAL_WINDOW_ID));

  ArrayColumn<Double> col_chan_freq(outSPW,
                                    MS_SPW::columnName(MS_SPW::CHAN_FREQ));
  ArrayColumn<Double> col_chan_widths(outSPW,
                                      MS_SPW::columnName(MS_SPW::CHAN_WIDTH));
  ArrayColumn<Double> col_effective_bw(
      outSPW, MS_SPW::columnName(MS_SPW::EFFECTIVE_BW));
  ScalarColumn<Int> col_meas_freq_ref(
      outSPW, MS_SPW::columnName(MS_SPW::MEAS_FREQ_REF));
  ScalarColumn<casacore::String> col_name(outSPW,
                                          MS_SPW::columnName(MS_SPW::NAME));
  ScalarColumn<Int> col_num_chan(outSPW, MS_SPW::columnName(MS_SPW::NUM_CHAN));
  ScalarColumn<Double> col_ref_freq(outSPW,
                                    MS_SPW::columnName(MS_SPW::REF_FREQUENCY));
  ArrayColumn<Double> col_resolution(outSPW,
                                     MS_SPW::columnName(MS_SPW::RESOLUTION));
  ScalarColumn<Double> col_total_bw(
      outSPW, MS_SPW::columnName(MS_SPW::TOTAL_BANDWIDTH));
  ScalarColumn<Int> col_bda_set_id(outSPW, base::DP3MS::kBDASetId);

  for (std::size_t i = 0; i < info().nbaselines(); ++i) {
    std::size_t nchanFreqs = info().chanFreqs(i).size();
    if (nchanToDescId.count(nchanFreqs)) {
      continue;
    }

    const common::rownr_t row = outDD.nrow();
    assert(row == outSPW.nrow());

    outDD.addRow();
    col_spw_id.put(row, row);

    outSPW.addRow();
    col_chan_freq.put(row, Vector<double>(info().chanFreqs(i)));
    col_chan_widths.put(row, Vector<double>(info().chanWidths(i)));
    col_effective_bw.put(row, Vector<double>(info().effectiveBW(i)));
    col_meas_freq_ref.put(row, measFreqRef);
    col_name.put(row, name);
    col_num_chan.put(row, nchanFreqs);
    col_ref_freq.put(row, info().refFreq());
    col_resolution.put(row, Vector<double>(info().resolutions(i)));
    col_total_bw.put(row, info().totalBW());
    col_bda_set_id.put(row, bda_set_id);

    nchanToDescId[nchanFreqs] = row;
  }
}

}  // namespace steps
}  // namespace dp3
