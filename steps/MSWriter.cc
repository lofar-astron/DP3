// MSWriter.cc: DPPP step writing to an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MSWriter.h"

#include <algorithm>
#include <iostream>
#include <limits>

#include <EveryBeam/correctionmode.h>

#include "InputStep.h"
#include "NullStep.h"

#include <Version.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../base/MS.h"

#include "../common/VdsMaker.h"
#include "../common/ParameterSet.h"

#include <casacore/tables/Tables/TableCopy.h>
#include <casacore/tables/Tables/TableLocker.h>
#include <casacore/tables/DataMan/DataManInfo.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledColumnStMan.h>
#include <casacore/tables/DataMan/TiledStManAccessor.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/TableMeasures/TableMeasDesc.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/casa/version.h>

#include <aocommon/logger.h>

using casacore::Array;
using casacore::ArrayColumn;
using casacore::ArrayColumnDesc;
using casacore::ArrayMeasColumn;
using casacore::Block;
using casacore::ColumnDesc;
using casacore::Cube;
using casacore::DataManager;
using casacore::DataManInfo;
using casacore::IPosition;
using casacore::Matrix;
using casacore::MDirection;
using casacore::MeasureHolder;
using casacore::Record;
using casacore::ScalarColumn;
using casacore::SetupNewTable;
using casacore::Table;
using casacore::TableCopy;
using casacore::TableDesc;
using casacore::TableRecord;
using casacore::TiledColumnStMan;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {
namespace {

/// Create an array column description and add to table with given
/// storage manager (if given).
void MakeArrayColumn(ColumnDesc desc, const IPosition& shape, DataManager* dm,
                     Table& table, bool make_direct_column = false) {
  desc.setOptions(0);
  int options = 0;
  if (!shape.empty()) {
    desc.setShape(shape);
    options |= ColumnDesc::FixedShape;
  }
  if (make_direct_column) {
    options |= ColumnDesc::Direct;
  }
  desc.setOptions(options);
  if (table.tableDesc().isColumn(desc.name())) {
    table.removeColumn(desc.name());
  }
  // Use storage manager if given.
  if (dm == nullptr) {
    table.addColumn(desc);
  } else {
    table.addColumn(desc, *dm);
  }
}

}  // namespace

MSWriter::MSWriter(const std::string& out_name,
                   const common::ParameterSet& parset,
                   const std::string& prefix)
    : name_(prefix),
      out_name_(out_name),
      parset_(parset),
      data_col_name_(parset.getString(prefix + "datacolumn", "DATA")),
      flag_col_name_(parset.getString(prefix + "flagcolumn", "FLAG")),
      weight_col_name_(
          parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM")),
      overwrite_(parset.getBool(prefix + "overwrite", false)),
      copy_corr_data_(parset.getBool(prefix + "copycorrecteddata", false)),
      copy_model_data_(parset.getBool(prefix + "copymodeldata", false)),
      // Get tile size (default 1024 KBytes).
      tile_size_(parset.getUint(prefix + "tilesize", 1024)),
      tile_n_chan_(parset.getUint(prefix + "tilenchan", 64)),
      nr_times_flush_(parset.getUint(prefix + "flush", 60)),
      nr_done_(0),
      chunk_duration_(parset.getDouble(prefix + "chunkduration", 0.0)),
      vds_dir_(parset.getString(prefix + "vdsdir", std::string())),
      cluster_desc_(parset.getString(prefix + "clusterdesc", std::string())),
      st_man_keys_(parset, prefix),
      scalar_flags_(
          parset.getBool(prefix + "scalarflags", METADATA_COMPRESSION_DEFAULT)),
      uvw_compression_(parset.getBool(prefix + "uvwcompression",
                                      METADATA_COMPRESSION_DEFAULT)),
      antenna_compression_(parset.getBool(prefix + "antennacompression",
                                          METADATA_COMPRESSION_DEFAULT)) {
  if (data_col_name_ != "DATA")
    throw std::runtime_error(
        "Currently only the DATA column"
        " can be used as output when writing a new MS");
  if (flag_col_name_ != "FLAG")
    throw std::runtime_error(
        "Currently only the FLAG column can be used as output for flags when "
        "writing a new MS");
  if (weight_col_name_ != "WEIGHT_SPECTRUM")
    throw std::runtime_error(
        "Currently only the "
        "WEIGHT_SPECTRUM column can be used as output when writing a new MS");
}

MSWriter::~MSWriter() { StopWriteThread(); }

std::string MSWriter::InsertNumberInFilename(const std::string& name,
                                             size_t number) {
  size_t dot = name.find_last_of('.');
  // If the name doesn't contain an extension, just add the number at the end
  if (dot == std::string::npos) {
    dot = name.size();
  }
  std::string n = std::to_string(number);
  n = number < 10 ? "00" + n : number < 100 ? "0" + n : n;
  return name.substr(0, dot) + '-' + n + name.substr(dot);
}

bool MSWriter::process(std::unique_ptr<DPBuffer> buffer) {
  if (chunk_start_time_ == 0.0) chunk_start_time_ = buffer->GetTime();

  if (chunk_duration_ != 0.0 &&
      buffer->GetTime() - chunk_start_time_ >= chunk_duration_) {
    FinishMs();
    ++current_chunk_index_;
    chunk_start_time_ = buffer->GetTime();
    StartNewMs();
  }

  common::NSTimer::StartStop sstime(timer_);

  if (use_write_thread_) {
    CreateTask(std::move(buffer));
  } else {
    ProcessBuffer(*buffer);
    getNextStep()->process(std::move(buffer));
  }

  return true;
}

void MSWriter::ProcessBuffer(DPBuffer& buffer) {
  const common::NSTimer::StartStop timer(writer_timer_);

  // Form the vector of the output table containing new rows.
  casacore::Vector<common::rownr_t> rownrs(getInfoOut().nbaselines());
  indgen(rownrs, ms_.nrow());
  // Add the necessary rows to the table.
  ms_.addRow(getInfoOut().nbaselines());
  // Form the subset of the tables containing the rows.
  // It can happen that a missing slot was inserted. In that case
  // the rownr vector is empty and we use the first itsNrBl input rows.
  // Time related info can only be copied if not averaging and if the
  // the time slot was not missing.
  Table out(ms_(rownrs));
  // Copy the input columns that do not change.
  WriteMeta(out, buffer);
  // Now write the data and flags.
  WriteData(out, buffer);
  // Flush if sufficient time slots are written.
  nr_done_++;
  if (nr_times_flush_ > 0 && nr_done_ % nr_times_flush_ == 0) {
    ms_.flush();
  }
  // Replace the rownrs in the buffer which is needed if in a later
  // step the MS gets updated.
  buffer.SetRowNumbers(rownrs);
}

void MSWriter::finish() {
  FinishMs();
  if (getNextStep()) getNextStep()->finish();
}

void MSWriter::StartNewMs() {
  // Create the MS.
  common::NSTimer::StartStop sstime(timer_);
  chunk_name_ = chunk_duration_ == 0.0
                    ? out_name_
                    : InsertNumberInFilename(out_name_, current_chunk_index_);
  CreateMs(chunk_name_, tile_size_, tile_n_chan_);
  // Write the parset info into the history.
  WriteHistory(ms_, parset_);
  ms_.flush(true, true);
  aocommon::Logger::Info << "Finished preparing output MS\n";
  GetWritableInfoOut().clearMetaChanged();

  use_write_thread_ = dynamic_cast<NullStep*>(getNextStep().get()) != nullptr;
  if (use_write_thread_) {
    is_write_thread_active_ = true;
    write_queue_thread_ = std::thread(&MSWriter::WriteQueueProcess, this);
  }
}

void MSWriter::FinishMs() {
  common::NSTimer::StartStop sstime(timer_);

  StopWriteThread();
  ms_.flush();

  // Create the VDS file.
  if (!cluster_desc_.empty()) {
    string vds_name = ms_.tableName() + ".vds";
    if (!vds_dir_.empty()) {
      if (vds_dir_[vds_dir_.size() - 1] != '/') {
        vds_dir_.append("/");
      }
      vds_name = vds_dir_ + string(casacore::Path(vds_name).baseName());
    }
    // Create VDS file without detailed time info.
    dp3::common::VdsMaker::create(ms_.tableName(), vds_name, cluster_desc_, "",
                                  false);
  }
  addToMS(chunk_name_);
}

void MSWriter::updateInfo(const DPInfo& info_in) {
  OutputStep::updateInfo(info_in);
  if (tile_n_chan_ <= 0 || tile_n_chan_ > info_in.nchan()) {
    tile_n_chan_ = info_in.nchan();
  }
  StartNewMs();
}

void MSWriter::show(std::ostream& os) const {
  os << "MSWriter " << name_ << '\n';
  os << "  output MS:      " << ms_.tableName() << '\n';
  os << "  nchan:          " << getInfoOut().nchan() << '\n';
  os << "  ncorrelations:  " << getInfoOut().ncorr() << '\n';
  os << "  nbaselines:     " << getInfoOut().nbaselines() << '\n';
  os << "  ntimes:         " << getInfoOut().ntime() << '\n';
  os << "  time interval:  " << getInfoOut().timeInterval() << '\n';
  os << "  DATA column:    " << data_col_name_ << '\n';
  os << "  FLAG column:    " << flag_col_name_ << '\n';
  os << "  WEIGHT column:  " << weight_col_name_ << '\n';
  os << GetCompressionString(st_man_keys_);
  os << "  scalar flags:   " << std::boolalpha << scalar_flags_ << '\n';
  os << "  use thread:     " << std::boolalpha << use_write_thread_ << '\n';
}

void MSWriter::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MSWriter " << name_ << '\n';

  duration = timer_.getElapsed();
  if (use_write_thread_) {
    os << "    ";
    FlagCounter::showPerc1(os, create_task_timer_.getElapsed(), duration);
    os << " Creating task\n";
  }
  os << "    ";
  FlagCounter::showPerc1(os, writer_timer_.getElapsed(), duration);
  os << (use_write_thread_ ? " Writing (threaded)\n" : " Writing\n");
}

void MSWriter::CreateMs(const std::string& out_name, unsigned int tile_size,
                        unsigned int tile_n_chan) {
  // Determine the data shape.
  IPosition data_shape(2, getInfoOut().ncorr(), getInfoOut().nchan());
  // Obtain the MS description.
  casacore::Table original_table(getInfoOut().msName());
  TableDesc tdesc(original_table.tableDesc());
  // Create the output table without the columns depending
  // on the nr of channels.
  // FLAG_CATEGORY is taken, but ignored when writing.
  std::vector<std::string> fixed_columns{
      "FLAG_CATEGORY", "WEIGHT",      "SIGMA",    "ARRAY_ID",
      "DATA_DESC_ID",  "EXPOSURE",    "FEED1",    "FEED2",
      "FIELD_ID",      "FLAG_ROW",    "INTERVAL", "OBSERVATION_ID",
      "PROCESSOR_ID",  "SCAN_NUMBER", "STATE_ID", "TIME",
      "TIME_CENTROID"};
  if (!uvw_compression_) fixed_columns.emplace_back("UVW");
  if (!antenna_compression_) {
    fixed_columns.emplace_back("ANTENNA1");
    fixed_columns.emplace_back("ANTENNA2");
  }
  Block<casacore::String> fixed_columns_block(fixed_columns.size());
  for (size_t i = 0; i != fixed_columns.size(); ++i)
    fixed_columns_block[i] = fixed_columns[i];
  Table temptable = original_table.project(fixed_columns_block);
  TableDesc newdesc = temptable.tableDesc();
  // Now quite some 'magic' is done to get the storage managers right.
  // Most of the columns do not change much and should be stored with
  // the IncrementalStMan. The ANTENNA columns change more often and
  // can best be stored with StandardStMan.
  // The new storage managers are only used for MSs stored with LofarStMan.
  // For 'normal' MSs they won't change.
  //
  // If needed, make WEIGHT and SIGMA a fixed shape, direct column.
  // In this way they are not written if the values are the same.
  {
    ColumnDesc& cdesc = newdesc.rwColumnDesc("WEIGHT");
    if (cdesc.shape().empty()) {
      cdesc.setShape(IPosition(1, getInfoOut().ncorr()), true);
    }
  }
  {
    ColumnDesc& cdesc = newdesc.rwColumnDesc("SIGMA");
    if (cdesc.shape().empty()) {
      cdesc.setShape(IPosition(1, getInfoOut().ncorr()), true);
    }
  }
  // Remove possible hypercolumn definitions.
// Test for casacore version 3.1.1 or smaller
#if CASACORE_MAJOR_VERSION < 3 ||    \
    (CASACORE_MAJOR_VERSION == 3 &&  \
     (CASACORE_MINOR_VERSION == 0 || \
      (CASACORE_MINOR_VERSION == 1 && CASACORE_PATCH_VERSION < 2)))
  newdesc.adjustHypercolumns(
      casacore::SimpleOrderedMap<casacore::String, casacore::String>(
          casacore::String()));
#else
  newdesc.adjustHypercolumns(std::map<casacore::String, casacore::String>());
#endif
  // Set data manager info.
  Record dminfo = temptable.dataManagerInfo();
  // Determine the DATA tile shape. Use all corrs and the given #channels.
  // The given tile size (in kbytes) determines the nr of rows in a tile .
  IPosition tile_shape(3, getInfoOut().ncorr(), tile_n_chan, 1);
  tile_shape[2] = tile_size * 1024 / (8 * tile_shape[0] * tile_shape[1]);
  if (tile_shape[2] < 1) {
    tile_shape[2] = 1;
  }
  // Replace all non-writable storage managers (i.e. LofarStMan) by ISM.
  dminfo = DataManInfo::adjustStMan(dminfo, "IncrementalStMan");

  if (!antenna_compression_) {
    // Remove ANTENNA1 and ANTENNA2 from the dminfo.
    // Don't remove them if already stored with StandardStMan.
    casacore::Vector<casacore::String> remove_cols(2);
    remove_cols[0] = "ANTENNA1";
    remove_cols[1] = "ANTENNA2";
    DataManInfo::removeDminfoColumns(dminfo, remove_cols, "StandardStMan");
  }
  // Configure UVW column.
  if (!uvw_compression_) {
    // Use as many rows as used for the DATA columns, but minimal 1024.
    int tsmnrow = tile_shape[2];
    if (tsmnrow < 1024) {
      tsmnrow = 1024;
    }
    // Make sure that the UVW column is stored as a fixed shape column.
    {
      ColumnDesc& cdesc = newdesc.rwColumnDesc("UVW");
      if (cdesc.shape().empty()) {
        cdesc.setShape(IPosition(1, 3), true);
      }
    }
    DataManInfo::setTiledStMan(
        dminfo, casacore::Vector<casacore::String>(1, "UVW"),
        "TiledColumnStMan", "TiledUVW", IPosition(2, 3, tsmnrow));
  }
  // Test if SSMVar already exists.
  bool has_ssm_var = false;
  for (unsigned int i = 0; i < dminfo.nfields(); ++i) {
    if (dminfo.subRecord(i).asString("NAME") == "SSMVar") {
      has_ssm_var = true;
      break;
    }
  }
  // Setup table creation. std::runtime_error is thrown if it exists already.
  Table::TableOption opt = overwrite_ ? Table::New : Table::NewNoReplace;
  SetupNewTable newtab(out_name, newdesc, opt);

  // First bind all columns to SSM.
  // For all columns defined in dminfo the bindings will be overwritten.
  // In this way variable columns like ANTENNA1/2 are bound to SSM.
  // Only do it if SSMVar does not exist (otherwise duplicate StMan name).
  if (!has_ssm_var) {
    casacore::StandardStMan ssm("SSMVar", 32768);
    newtab.bindAll(ssm);
  }

  // Bind all columns according to dminfo.
  newtab.bindCreate(dminfo);
  ms_ = Table(newtab);

  if (antenna_compression_) {
    std::unique_ptr<DataManager> ant_st_man =
        MakeStMan("AntennaPairStMan", "AntennaPairStMan");
    ms_.addColumn(tdesc["ANTENNA1"], *ant_st_man);
    ms_.addColumn(tdesc["ANTENNA2"], "AntennaPairStMan", false);
  }

  const std::string stdman_class_name(
      st_man_keys_.GetStorageManagerClassName());
  const bool use_custom_stman = !stdman_class_name.empty();
  const casacore::DataManagerCtor data_constructor =
      use_custom_stman ? DataManager::getCtor(stdman_class_name) : nullptr;
  const bool skip_dysco_for_data =
      st_man_keys_.storage_manager_name == "dysco" &&
      st_man_keys_.dysco_data_bit_rate == 0;
  if (use_custom_stman && !skip_dysco_for_data) {
    const Record data_man_specification = st_man_keys_.GetSpecification();
    std::unique_ptr<DataManager> data_manager(
        data_constructor("DataStMan", data_man_specification));
    MakeArrayColumn(tdesc["DATA"], data_shape, data_manager.get(), ms_, true);
  } else {
    TiledColumnStMan tsm("TiledData", tile_shape);
    MakeArrayColumn(tdesc["DATA"], data_shape, &tsm, ms_);
  }

  if (uvw_compression_) {
    std::unique_ptr<DataManager> uvw_st_man =
        MakeStMan("UvwStMan", "CompressedUvw");
    MakeArrayColumn(tdesc["UVW"], casacore::IPosition{3}, uvw_st_man.get(), ms_,
                    true);
  }

  // Add FLAG column
  if (scalar_flags_) {
    std::unique_ptr<DataManager> stokes_i_st_man =
        MakeStMan("StokesIStMan", "StokesIFlag");
    MakeArrayColumn(tdesc["FLAG"], data_shape, stokes_i_st_man.get(), ms_,
                    true);
  } else {
    // Use larger tile shape because flags are stored as bits.
    IPosition tile_shape_f(tile_shape);
    tile_shape_f[2] *= 8;
    TiledColumnStMan tsmf("TiledFlag", tile_shape_f);
    MakeArrayColumn(tdesc["FLAG"], data_shape, &tsmf, ms_);
  }

  if (st_man_keys_.storage_manager_name == "dysco" &&
      st_man_keys_.dysco_weight_bit_rate != 0) {
    // Add WEIGHT_SPECTRUM column using Dysco stman.
    std::unique_ptr<DataManager> dysco_st_man(
        data_constructor("DyscoWeightSpectrum", st_man_keys_.GetDyscoSpec()));
    ArrayColumnDesc<float> wsdesc("WEIGHT_SPECTRUM", "weight per corr/chan",
                                  data_shape,
                                  ColumnDesc::FixedShape | ColumnDesc::Direct);
    MakeArrayColumn(wsdesc, data_shape, dysco_st_man.get(), ms_, true);
  } else {
    // Add WEIGHT_SPECTRUM column using tsm.
    TiledColumnStMan tsmw("TiledWeightSpectrum", tile_shape);
    ArrayColumnDesc<float> wsdesc("WEIGHT_SPECTRUM", "weight per corr/chan",
                                  data_shape, ColumnDesc::FixedShape);
    MakeArrayColumn(wsdesc, data_shape, &tsmw, ms_);
  }
  // If present handle the CORRECTED_DATA and MODEL_DATA column.
  if (!tdesc.isColumn("CORRECTED_DATA")) {
    copy_corr_data_ = false;
  }
  if (!tdesc.isColumn("MODEL_DATA")) {
    copy_model_data_ = false;
  }
  if (copy_corr_data_) {
    TiledColumnStMan tsmc("CorrectedData", tile_shape);
    MakeArrayColumn(tdesc["CORRECTED_DATA"], data_shape, &tsmc, ms_);

    IPosition iw_shape(1, data_shape[1]);
    IPosition iw_shape_tile(2, tile_shape[1], tile_shape[2]);
    TiledColumnStMan tsmw("TiledImagingWeight", iw_shape_tile);
    ColumnDesc iwdesc(ArrayColumnDesc<float>("IMAGING_WEIGHT"));
    MakeArrayColumn(iwdesc, iw_shape, &tsmw, ms_);
  }
  if (copy_model_data_) {
    ColumnDesc mdesc = tdesc.columnDesc("MODEL_DATA");
    TableRecord& keyset = mdesc.rwKeywordSet();
    // Redefine possible keywords used by the CASA VisSet classes.
    if (keyset.isDefined("CHANNEL_SELECTION")) {
      keyset.removeField("CHANNEL_SELECTION");
    }
    Matrix<int> selection(2, 1);
    selection(0, 0) = 0;
    selection(1, 0) = getInfoOut().nchan();
    keyset.define("CHANNEL_SELECTION", selection);
    TiledColumnStMan tsmm("ModelData", tile_shape);
    MakeArrayColumn(mdesc, data_shape, &tsmm, ms_);
  }
  aocommon::Logger::Info << " copying info and subtables ...\n";
  // Copy the info and subtables.
  TableCopy::copyInfo(ms_, temptable);
  TableRecord& keyset = ms_.rwKeywordSet();
  if (keyset.isDefined(base::DP3MS::kBDAFactorsTable)) {
    keyset.removeField(base::DP3MS::kBDAFactorsTable);
  }
  CopySubTables(original_table);

  // Adjust the SPECTRAL_WINDOW and DATA_DESCRIPTION table as needed.
  UpdateSpw(out_name, getInfoOut());
  // Adjust the OBSERVATION table as needed.
  UpdateObs(out_name, getInfoOut());
  // Adjust the FIELD table as needed.
  if (getInfoOut().originalPhaseCenter().getValue() !=
      getInfoOut().phaseCenter().getValue()) {
    UpdatePhaseCentre(out_name, getInfoOut().phaseCenter());
  }
  UpdateBeam(ms_, "DATA", getInfoOut());
}

void MSWriter::CopySubTables(casacore::Table& original_table) {
  // AST-1186: This function is a re-implementation of
  // casacore::TableCopy::copySubTables.
  // When a subtable contains symbolic links, casacore::TableCopy::copySubTable
  // copies the symbolic links. This function copies the subtables by value.
  const bool kValueCopy = true;

  const casacore::TableRecord& keys_in = original_table.keywordSet();
  casacore::TableRecord& keys_out = ms_.rwKeywordSet();

  for (casacore::uInt i = 0; i < keys_in.nfields(); ++i) {
    if (keys_in.type(i) != casacore::TpTable) continue;

    const std::string subtable_name = keys_in.name(i);

    // Skip a subtable that has to be omitted.
    if (subtable_name == base::DP3MS::kBDATimeAxisTable ||
        subtable_name == base::DP3MS::kBDAFactorsTable)
      continue;

    Table subtable_in = keys_in.asTable(i);

    // Lock the subtable / keep the lock if already locked.
    casacore::TableLocker locker(subtable_in, casacore::FileLocker::Read);

    const casacore::String subtable_path_out =
        ms_.tableName() + "/" + subtable_name.c_str();
    subtable_in.deepCopy(subtable_path_out, Table::New, kValueCopy);

    keys_out.defineTable(subtable_name, Table(subtable_path_out));
  }
}

void MSWriter::UpdateSpw(const std::string& out_name,
                         const base::DPInfo& info) {
  // Fix the SPECTRAL_WINDOW values by updating the values in the subtable.
  IPosition shape(1, info.nchan());
  Table original_table(info.msName());
  Table in_spw = original_table.keywordSet().asTable("SPECTRAL_WINDOW");
  Table out_spw = Table(out_name + "/SPECTRAL_WINDOW", Table::Update);
  Table out_dd = Table(out_name + "/DATA_DESCRIPTION", Table::Update);
  if (out_spw.nrow() != out_dd.nrow())
    throw std::runtime_error(
        "nrow in SPECTRAL_WINDOW table is not the same as nrow in "
        "DATA_DESCRIPTION table");
  // Remove all rows before and after the selected band.
  // Do it from the end, otherwise row numbers change.
  for (int i = int(out_spw.nrow()) - 1; i >= 0; --i) {
    if (i != info.spectralWindow()) {
      out_spw.removeRow(i);
      out_dd.removeRow(i);
    }
  }
  // Set nr of channels.
  ScalarColumn<int> channum(out_spw, "NUM_CHAN");
  channum.fillColumn(info.nchan());
  // Change the column shapes.
  TableDesc tdesc = in_spw.tableDesc();
  MakeArrayColumn(tdesc["CHAN_FREQ"], shape, nullptr, out_spw);
  MakeArrayColumn(tdesc["CHAN_WIDTH"], shape, nullptr, out_spw);
  MakeArrayColumn(tdesc["EFFECTIVE_BW"], shape, nullptr, out_spw);
  MakeArrayColumn(tdesc["RESOLUTION"], shape, nullptr, out_spw);
  // Create the required column objects.
  ArrayColumn<double> out_freq(out_spw, "CHAN_FREQ");
  ArrayColumn<double> out_width(out_spw, "CHAN_WIDTH");
  ArrayColumn<double> out_bw(out_spw, "EFFECTIVE_BW");
  ArrayColumn<double> out_resolution(out_spw, "RESOLUTION");
  ScalarColumn<double> out_totalbw(out_spw, "TOTAL_BANDWIDTH");
  ScalarColumn<double> out_reffreq(out_spw, "REF_FREQUENCY");
  out_freq.put(0, casacore::Vector<double>(info.chanFreqs()));
  out_width.put(0, casacore::Vector<double>(info.chanWidths()));
  out_bw.put(0, casacore::Vector<double>(info.effectiveBW()));
  out_resolution.put(0, casacore::Vector<double>(info.resolutions()));
  out_totalbw.put(0, info.totalBW());
  out_reffreq.put(0, info.refFreq());
  // Adjust the spwid in the DATA_DESCRIPTION.
  ScalarColumn<int> spw_col(out_dd, "SPECTRAL_WINDOW_ID");
  spw_col.put(0, 0);
}

void MSWriter::UpdateObs(const std::string& out_name,
                         const base::DPInfo& info) {
  Table out_obs = Table(out_name + "/OBSERVATION", Table::Update);
  // Set nr of channels.
  ArrayColumn<double> time_range(out_obs, "TIME_RANGE");
  casacore::Vector<double> times(2);
  times[0] = info.firstTime() - 0.5 * info.timeInterval();
  times[1] = info.lastTime() + 0.5 * info.timeInterval();
  // There should be one row, but loop in case of.
  for (unsigned int i = 0; i < out_obs.nrow(); ++i) {
    time_range.put(i, times);
  }
}

void MSWriter::UpdatePhaseCentre(const std::string& out_name,
                                 const casacore::MDirection& new_phase_dir) {
  Table out_field = Table(out_name + "/FIELD", Table::Update);
  // Write new phase center.
  ArrayMeasColumn<MDirection> phase_col(out_field, "PHASE_DIR");
  // If a moving reference type like AZELGEO was used in the original MS, and
  // the phase centre is changed (with a phaseshift), the ref frame of the
  // column must be reset:
  phase_col.setDescRefCode(new_phase_dir.getRefPtr()->getType(), false);
  const casacore::Vector<MDirection> dir(1, new_phase_dir);
  phase_col.put(0, dir);
}

void MSWriter::UpdateBeam(Table& main_table, const std::string& out_col_name,
                          const DPInfo& info) {
  const char *beam_mode_field_name = "LOFAR_APPLIED_BEAM_MODE",
             *beam_dir_field_name = "LOFAR_APPLIED_BEAM_DIR";

  ArrayColumn<casacore::Complex> data_column(main_table, out_col_name);
  bool fields_exist = data_column.keywordSet().isDefined(beam_mode_field_name);
  const std::string mode_str = everybeam::ToString(
      static_cast<everybeam::CorrectionMode>(info.beamCorrectionMode()));
  // If no beam correction has been applied and the LOFAR beam fields don't
  // exist, we have to do nothing (no fields implies no beam correction).
  // If they do exist, we have to make sure they are set to indicate
  // no beam correction.
  if (fields_exist || info.beamCorrectionMode() !=
                          static_cast<int>(everybeam::CorrectionMode::kNone)) {
    data_column.rwKeywordSet().define(beam_mode_field_name, mode_str);
    Record record;
    MeasureHolder(info.beamCorrectionDir()).toRecord(record);
    data_column.rwKeywordSet().defineRecord(beam_dir_field_name, record);
  }
}

void MSWriter::WriteHistory(Table& ms, const common::ParameterSet& parset) {
  Table histtab(ms.keywordSet().asTable("HISTORY"));
  histtab.reopenRW();
  ScalarColumn<double> time(histtab, "TIME");
  ScalarColumn<int> obs_id(histtab, "OBSERVATION_ID");
  ScalarColumn<casacore::String> message(histtab, "MESSAGE");
  ScalarColumn<casacore::String> application(histtab, "APPLICATION");
  ScalarColumn<casacore::String> priority(histtab, "PRIORITY");
  ScalarColumn<casacore::String> origin(histtab, "ORIGIN");
  ArrayColumn<casacore::String> parms(histtab, "APP_PARAMS");
  ArrayColumn<casacore::String> cli(histtab, "CLI_COMMAND");
  // Put all parset entries in a Vector<casacore::String>.
  // Some WSRT MSs have a FixedShape APP_PARAMS and CLI_COMMAND column.
  // For them, put all params in a single vector element (with newlines).
  bool fixed_shaped =
      (parms.columnDesc().options() & ColumnDesc::FixedShape) != 0;
  casacore::Vector<casacore::String> appvec;
  casacore::Vector<casacore::String> clivec;
  if (fixed_shaped) {
    appvec.resize(1);
    clivec.resize(1);
    std::ostringstream ostr;
    parset.writeStream(ostr);
    appvec[0] = ostr.str();
  } else {
    appvec.resize(parset.size());
    Array<casacore::String>::contiter viter = appvec.cbegin();
    for (common::ParameterSet::const_iterator iter = parset.begin();
         iter != parset.end(); ++iter, ++viter) {
      *viter = iter->first + '=' + iter->second.get();
    }
  }
  unsigned int rownr = histtab.nrow();
  histtab.addRow();
  time.put(rownr, casacore::Time().modifiedJulianDay() * 24. * 3600.);
  obs_id.put(rownr, 0);
  message.put(rownr, "parameters");
  application.put(rownr, "DP3");
  priority.put(rownr, "NORMAL");
  origin.put(rownr, DP3Version::AsString());
  parms.put(rownr, appvec);
  cli.put(rownr, clivec);
}

void MSWriter::WriteData(Table& out, DPBuffer& buf) {
  if (buf.GetData().size() == 0) {
    return;
  }

  // If compressing, flagged values need to be set to NaN, and flagged
  // weights to zero, to decrease the dynamic range
  if (st_man_keys_.storage_manager_name == "dysco") {
    auto data_iterator = buf.GetData().begin();
    auto weights_iterator = buf.GetWeights().begin();
    for (bool flag : buf.GetFlags()) {
      if (flag) {
        *data_iterator =
            casacore::Complex(std::numeric_limits<float>::quiet_NaN(),
                              std::numeric_limits<float>::quiet_NaN());
        *weights_iterator = 0.0;
      }
      ++data_iterator;
      ++weights_iterator;
    }
  }

  // Write DATA, WEIGHT_SPECTRUM and FLAG
  ArrayColumn<casacore::Complex> data_col(out, data_col_name_);
  ArrayColumn<bool> flag_col(out, "FLAG");
  ArrayColumn<float> weightCol(out, "WEIGHT_SPECTRUM");
  const casacore::IPosition shape(3, getInfoOut().ncorr(), getInfoOut().nchan(),
                                  getInfoOut().nbaselines());
  const Cube<casacore::Complex> data(shape, buf.GetData().data(),
                                     casacore::SHARE);
  const Cube<float> weights(shape, buf.GetWeights().data(), casacore::SHARE);
  const Cube<bool> flags(shape, buf.GetFlags().data(), casacore::SHARE);
  data_col.putColumn(data);
  weightCol.putColumn(weights);
  flag_col.putColumn(flags);

  // A row is flagged if no flags in the row are False.
  auto c = partialNFalse(flags, IPosition(2, 0, 1));
  casacore::Vector<bool> row_flags(c == decltype(c)::value_type(0));
  ScalarColumn<bool> flag_row_col(out, "FLAG_ROW");
  flag_row_col.putColumn(row_flags);

  // Write UVW
  ArrayColumn<double> uvw_col(out, "UVW");
  const casacore::IPosition shape_uvw(2, 3, getInfoOut().nbaselines());
  const Matrix<double> uvws(shape_uvw, buf.GetUvw().data(), casacore::SHARE);
  uvw_col.putColumn(uvws);
}

void MSWriter::WriteMeta(Table& out, const DPBuffer& buf) {
  // Fill ANTENNA1/2.
  ScalarColumn<int> ant1col(out, "ANTENNA1");
  ScalarColumn<int> ant2col(out, "ANTENNA2");
  if (antenna_compression_) {
    // When using the compressing mgr, antenna 1 and 2 have to be written
    // directly after each other.
    const size_t n = getInfoOut().getAnt1().size();
    for (size_t i = 0; i != n; ++i) {
      ant1col.put(i, getInfoOut().getAnt1()[i]);
      ant2col.put(i, getInfoOut().getAnt2()[i]);
    }
  } else {
    ant1col.putColumn(casacore::Vector<int>(getInfoOut().getAnt1()));
    ant2col.putColumn(casacore::Vector<int>(getInfoOut().getAnt2()));
  }
  // Fill all rows that do not change.
  FillSca<double>(buf.GetTime(), out, "TIME");
  FillSca<double>(buf.GetTime(), out, "TIME_CENTROID");
  FillSca<double>(buf.GetExposure(), out, "EXPOSURE");
  FillSca<double>(getInfoOut().timeInterval(), out, "INTERVAL");
  FillSca<int>(0, out, "FEED1");
  FillSca<int>(0, out, "FEED2");
  FillSca<int>(0, out, "DATA_DESC_ID");
  FillSca<int>(0, out, "PROCESSOR_ID");
  FillSca<int>(0, out, "FIELD_ID");
  FillSca<int>(0, out, "SCAN_NUMBER");
  FillSca<int>(0, out, "ARRAY_ID");
  FillSca<int>(0, out, "OBSERVATION_ID");
  FillSca<int>(0, out, "STATE_ID");
  Array<float> arr(IPosition(1, getInfoOut().ncorr()));
  arr = 1;
  FillArr<float>(arr, out, "SIGMA");
  FillArr<float>(arr, out, "WEIGHT");
}

void MSWriter::CopyMeta(const Table& in, Table& out, bool copy_time_info) {
  // Copy all rows that do not change.
  CopySca<int>(in, out, "ANTENNA1");
  CopySca<int>(in, out, "ANTENNA2");
  CopySca<int>(in, out, "FEED1");
  CopySca<int>(in, out, "FEED2");
  CopySca<int>(in, out, "PROCESSOR_ID");
  CopySca<int>(in, out, "FIELD_ID");
  CopySca<int>(in, out, "SCAN_NUMBER");
  CopySca<int>(in, out, "ARRAY_ID");
  CopySca<int>(in, out, "OBSERVATION_ID");
  CopySca<int>(in, out, "STATE_ID");
  CopyArr<float>(in, out, "SIGMA");
  CopyArr<float>(in, out, "WEIGHT");
  if (copy_time_info) {
    CopySca<double>(in, out, "TIME");
    CopySca<double>(in, out, "TIME_CENTROID");
    CopySca<double>(in, out, "INTERVAL");
    CopySca<double>(in, out, "EXPOSURE");
    CopyArr<double>(in, out, "UVW");
  }
}

void MSWriter::StopWriteThread() {
  if (!is_write_thread_active_) {
    return;
  }
  write_queue_.write_end();
  write_queue_thread_.join();
  write_queue_.clear();
  is_write_thread_active_ = false;
}

void MSWriter::WriteQueueProcess() {
  std::unique_ptr<base::DPBuffer> buffer;
  while (write_queue_.read(buffer)) {
    ProcessBuffer(*buffer);
  }
}

void MSWriter::CreateTask(std::unique_ptr<base::DPBuffer> buffer) {
  const common::NSTimer::StartStop timer(create_task_timer_);

  write_queue_.write(std::move(buffer));
}

std::unique_ptr<DataManager> MakeStMan(const std::string& type_name,
                                       const std::string& instance_name,
                                       const Record& record) {
  casacore::DataManagerCtor constructor = DataManager::getCtor(type_name);
  std::unique_ptr<DataManager> st_man(constructor(instance_name, record));
  if (!st_man)
    throw std::runtime_error("Storage manager " + type_name +
                             " requested, but it is not available in "
                             "casacore.");
  return st_man;
}

}  // namespace steps
}  // namespace dp3
