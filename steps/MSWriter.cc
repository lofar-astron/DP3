// MSWriter.cc: DPPP step writing to an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MSWriter.h"

#include "InputStep.h"
#include "NullStep.h"

#include <Version.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../base/DPLogger.h"
#include "../base/MS.h"

#include "../common/VdsMaker.h"
#include "../common/ParameterSet.h"

#include <casacore/tables/Tables/TableCopy.h>
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
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/casa/version.h>

#include <iostream>
#include <limits>

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

MSWriter::MSWriter(InputStep& reader, const std::string& out_name,
                   const common::ParameterSet& parset,
                   const std::string& prefix)
    : reader_(reader),
      name_(prefix),
      out_name_(out_name),
      parset_(parset),
      nr_done_(0) {
  // Get tile size (default 1024 KBytes).
  tile_size_ = parset.getUint(prefix + "tilesize", 1024);
  tile_n_chan_ = parset.getUint(prefix + "tilenchan", 0);
  chunk_duration_ = parset.getDouble(prefix + "chunkduration", 0.0);
  overwrite_ = parset.getBool(prefix + "overwrite", false);
  nr_times_flush_ = parset.getUint(prefix + "flush", 60);
  copy_corr_data_ = parset.getBool(prefix + "copycorrecteddata", false);
  copy_model_data_ = parset.getBool(prefix + "copymodeldata", false);
  write_full_res_flags_ = parset.getBool(prefix + "writefullresflag", true);
  data_col_name_ = parset.getString(prefix + "datacolumn", "DATA");
  flag_col_name_ = parset.getString(prefix + "flagcolumn", "FLAG");
  weight_col_name_ =
      parset.getString(prefix + "weightcolumn", "WEIGHT_SPECTRUM");
  vds_dir_ = parset.getString(prefix + "vdsdir", string());
  cluster_desc_ = parset.getString(prefix + "clusterdesc", string());
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

  st_man_keys_.Set(parset, prefix);
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

bool MSWriter::process(const DPBuffer& buf) {
  if (chunk_start_time_ == 0.0) chunk_start_time_ = buf.getTime();

  if (chunk_duration_ != 0.0 &&
      buf.getTime() - chunk_start_time_ >= chunk_duration_) {
    FinishMs();
    ++current_chunk_index_;
    chunk_start_time_ = buf.getTime();
    StartNewMs();
  }

  common::NSTimer::StartStop sstime(timer_);

  UpdateInternalBuffer(buf);
  if (use_write_thread_) {
    CreateTask();
  } else {
    ProcessBuffer(internal_buffer_);
    getNextStep()->process(internal_buffer_);
  }

  return true;
}

void MSWriter::ProcessBuffer(DPBuffer& buffer) {
  const common::NSTimer::StartStop timer(writer_timer_);

  // Form the vector of the output table containing new rows.
  casacore::Vector<common::rownr_t> rownrs(nr_bl_);
  indgen(rownrs, ms_.nrow());
  // Add the necessary rows to the table.
  ms_.addRow(nr_bl_);
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
  buffer.setRowNrs(rownrs);
}

void MSWriter::finish() {
  FinishMs();
  getNextStep()->finish();
}

void MSWriter::addToMS(const string&) { getPrevStep()->addToMS(out_name_); }

void MSWriter::StartNewMs() {
  // Create the MS.
  common::NSTimer::StartStop sstime(timer_);
  const std::string chunk_name =
      chunk_duration_ == 0.0
          ? out_name_
          : InsertNumberInFilename(out_name_, current_chunk_index_);
  CreateMs(chunk_name, info(), tile_size_, tile_n_chan_);
  // Write the parset info into the history.
  WriteHistory(ms_, parset_);
  ms_.flush(true, true);
  DPLOG_INFO("Finished preparing output MS", false);
  info().clearMetaChanged();

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
}

void MSWriter::updateInfo(const DPInfo& info_in) {
  info() = info_in;
  interval_ = info().timeInterval();
  nr_corr_ = info().ncorr();
  nr_chan_ = info().nchan();
  nr_bl_ = info().nbaselines();
  nr_times_ = info().ntime();
  // Input can already be averaged, so take that into account.
  n_chan_avg_ = reader_.nchanAvgFullRes() * info().nchanAvg();
  n_time_avg_ = reader_.ntimeAvgFullRes() * info().ntimeAvg();
  if (tile_n_chan_ <= 0 || tile_n_chan_ > getInfo().nchan()) {
    tile_n_chan_ = getInfo().nchan();
  }
  StartNewMs();
}

void MSWriter::show(std::ostream& os) const {
  os << "MSWriter " << name_ << '\n';
  os << "  output MS:      " << ms_.tableName() << '\n';
  os << "  nchan:          " << nr_chan_ << '\n';
  os << "  ncorrelations:  " << nr_corr_ << '\n';
  os << "  nbaselines:     " << nr_bl_ << '\n';
  os << "  ntimes:         " << nr_times_ << '\n';
  os << "  time interval:  " << interval_ << '\n';
  os << "  DATA column:    " << data_col_name_ << '\n';
  os << "  FLAG column:    " << flag_col_name_ << '\n';
  os << "  WEIGHT column:  " << weight_col_name_ << '\n';
  if (st_man_keys_.stManName == "dysco") {
    os << "  Compressed:     yes\n"
       << "  Data bitrate:   " << st_man_keys_.dyscoDataBitRate << '\n'
       << "  Weight bitrate: " << st_man_keys_.dyscoWeightBitRate << '\n'
       << "  Dysco mode:     " << st_man_keys_.dyscoNormalization << ' '
       << st_man_keys_.dyscoDistribution << '('
       << st_man_keys_.dyscoDistTruncation << ")\n";
  } else {
    os << "  Compressed:     no\n";
  }
  os << "  use thread:     " << std::boolalpha << use_write_thread_ << '\n';
}

void MSWriter::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " MSWriter " << name_ << '\n';

  duration = timer_.getElapsed();
  os << "    ";
  FlagCounter::showPerc1(os, update_buffer_timer_.getElapsed(), duration);
  os << " Updating buffer\n";
  if (use_write_thread_) {
    os << "    ";
    FlagCounter::showPerc1(os, create_task_timer_.getElapsed(), duration);
    os << " Creating task\n";
  }
  os << "    ";
  FlagCounter::showPerc1(os, writer_timer_.getElapsed(), duration);
  os << (use_write_thread_ ? " Writing (threaded)\n" : " Writing\n");
}

void MSWriter::MakeArrayColumn(ColumnDesc desc, const IPosition& ipos,
                               DataManager* dm, Table& table,
                               bool make_direct_column) {
  desc.setOptions(0);
  desc.setShape(ipos);
  if (make_direct_column) {
    desc.setOptions(ColumnDesc::Direct | ColumnDesc::FixedShape);
  } else {
    desc.setOptions(ColumnDesc::FixedShape);
  }
  if (table.tableDesc().isColumn(desc.name())) {
    table.removeColumn(desc.name());
  }
  // Use storage manager if given.
  if (dm == 0) {
    table.addColumn(desc);
  } else {
    table.addColumn(desc, *dm);
  }
}

void MSWriter::CreateMs(const string& out_name, const DPInfo& info,
                        unsigned int tile_size, unsigned int tile_n_chan) {
  // Determine the data shape.
  IPosition data_shape(2, nr_corr_, nr_chan_);
  // Obtain the MS description.
  TableDesc tdesc(reader_.table().tableDesc());
  // Create the output table without the columns depending
  // on the nr of channels.
  // FLAG_CATEGORY is taken, but ignored when writing.
  Block<casacore::String> fixed_columns(20);
  fixed_columns[0] = "UVW";
  fixed_columns[1] = "FLAG_CATEGORY";
  fixed_columns[2] = "WEIGHT";
  fixed_columns[3] = "SIGMA";
  fixed_columns[4] = "ANTENNA1";
  fixed_columns[5] = "ANTENNA2";
  fixed_columns[6] = "ARRAY_ID";
  fixed_columns[7] = "DATA_DESC_ID";
  fixed_columns[8] = "EXPOSURE";
  fixed_columns[9] = "FEED1";
  fixed_columns[10] = "FEED2";
  fixed_columns[11] = "FIELD_ID";
  fixed_columns[12] = "FLAG_ROW";
  fixed_columns[13] = "INTERVAL";
  fixed_columns[14] = "OBSERVATION_ID";
  fixed_columns[15] = "PROCESSOR_ID";
  fixed_columns[16] = "SCAN_NUMBER";
  fixed_columns[17] = "STATE_ID";
  fixed_columns[18] = "TIME";
  fixed_columns[19] = "TIME_CENTROID";
  Table temptable = reader_.table().project(fixed_columns);
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
      cdesc.setShape(IPosition(1, nr_corr_), true);
    }
  }
  {
    ColumnDesc& cdesc = newdesc.rwColumnDesc("SIGMA");
    if (cdesc.shape().empty()) {
      cdesc.setShape(IPosition(1, nr_corr_), true);
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
  IPosition tile_shape(3, nr_corr_, tile_n_chan, 1);
  tile_shape[2] = tile_size * 1024 / (8 * tile_shape[0] * tile_shape[1]);
  if (tile_shape[2] < 1) {
    tile_shape[2] = 1;
  }
  // Replace all non-writable storage managers (i.e. LofarStMan) by ISM.
  dminfo = DataManInfo::adjustStMan(dminfo, "IncrementalStMan");
  // Remove ANTENNA1 and ANTENNA2 from the dminfo.
  // Don't remove them if already stored with StandardStMan.
  casacore::Vector<casacore::String> remove_cols(2);
  remove_cols[0] = "ANTENNA1";
  remove_cols[1] = "ANTENNA2";
  DataManInfo::removeDminfoColumns(dminfo, remove_cols, "StandardStMan");
  // Use TiledStMan for UVW.
  // Use as many rows as used for the DATA columns, but minimal 1024.
  int tsmnrow = tile_shape[2];
  if (tsmnrow < 1024) {
    tsmnrow = 1024;
  }
  DataManInfo::setTiledStMan(
      dminfo, casacore::Vector<casacore::String>(1, "UVW"), "TiledColumnStMan",
      "TiledUVW", IPosition(2, 3, tsmnrow));
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
  casacore::DataManagerCtor dysco_constructor = 0;
  Record dysco_spec;
  if (st_man_keys_.stManName == "dysco") {
    dysco_spec = st_man_keys_.GetDyscoSpec();
    dysco_constructor = DataManager::getCtor("DyscoStMan");
  }
  ms_ = Table(newtab);

  if (st_man_keys_.stManName == "dysco" && st_man_keys_.dyscoDataBitRate != 0) {
    // Add DATA column using Dysco stman.
    std::unique_ptr<DataManager> dysco_st_man(
        dysco_constructor("DyscoData", dysco_spec));
    MakeArrayColumn(tdesc["DATA"], data_shape, dysco_st_man.get(), ms_, true);
  } else {
    // Add DATA column using tsm.
    TiledColumnStMan tsm("TiledData", tile_shape);
    MakeArrayColumn(tdesc["DATA"], data_shape, &tsm, ms_);
  }

  // Add FLAG column using tsm.
  // Use larger tile shape because flags are stored as bits.
  IPosition tile_shape_f(tile_shape);
  tile_shape_f[2] *= 8;
  TiledColumnStMan tsmf("TiledFlag", tile_shape_f);
  MakeArrayColumn(tdesc["FLAG"], data_shape, &tsmf, ms_);

  if (write_full_res_flags_) {
    // Add LOFAR_FULL_RES_FLAG column using tsm.
    // The input can already be averaged and averaging can be done in
    // this run, so the full resolution is the combination of both.
    unsigned int orignchan = nr_chan_ * n_chan_avg_;
    IPosition data_shape_f(2, (orignchan + 7) / 8, n_time_avg_);
    IPosition tile_shape_f(3, (orignchan + 7) / 8, 1024, tile_shape[2]);
    TiledColumnStMan tsmf("TiledFullResFlag", tile_shape_f);
    ArrayColumnDesc<unsigned char> padesc("LOFAR_FULL_RES_FLAG",
                                          "flags in original full resolution",
                                          data_shape_f, ColumnDesc::FixedShape);
    MakeArrayColumn(padesc, data_shape_f, &tsmf, ms_);
  }
  if (st_man_keys_.stManName == "dysco" &&
      st_man_keys_.dyscoWeightBitRate != 0) {
    // Add WEIGHT_SPECTRUM column using Dysco stman.
    std::unique_ptr<DataManager> dysco_st_man(
        dysco_constructor("DyscoWeightSpectrum", dysco_spec));
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
    selection(1, 0) = nr_chan_;
    keyset.define("CHANNEL_SELECTION", selection);
    TiledColumnStMan tsmm("ModelData", tile_shape);
    MakeArrayColumn(mdesc, data_shape, &tsmm, ms_);
  }
  DPLOG_INFO(" copying info and subtables ...", false);
  // Copy the info and subtables.
  TableCopy::copyInfo(ms_, temptable);
  TableRecord& keyset = ms_.rwKeywordSet();
  if (keyset.isDefined(base::DP3MS::kBDAFactorsTable)) {
    keyset.removeField(base::DP3MS::kBDAFactorsTable);
  }
  casacore::Block<casacore::String> omitted_subtables(2);
  omitted_subtables[0] = base::DP3MS::kBDATimeAxisTable;
  omitted_subtables[1] = base::DP3MS::kBDAFactorsTable;
  TableCopy::copySubTables(ms_, temptable, false, omitted_subtables);
  // Adjust the SPECTRAL_WINDOW and DATA_DESCRIPTION table as needed.
  UpdateSpw(out_name, info);
  // Adjust the OBSERVATION table as needed.
  UpdateObs(out_name);
  // Adjust the FIELD table as needed.
  if (!info.phaseCenterIsOriginal()) {
    UpdatePhaseCentre(out_name, info);
  }
  UpdateBeam(out_name, "DATA", info);
}

void MSWriter::UpdateSpw(const string& out_name, const DPInfo& info) {
  // Fix the SPECTRAL_WINDOW values by updating the values in the subtable.
  IPosition shape(1, nr_chan_);
  Table in_spw = reader_.table().keywordSet().asTable("SPECTRAL_WINDOW");
  Table out_spw = Table(out_name + "/SPECTRAL_WINDOW", Table::Update);
  Table out_dd = Table(out_name + "/DATA_DESCRIPTION", Table::Update);
  if (out_spw.nrow() != out_dd.nrow())
    throw std::runtime_error(
        "nrow in SPECTRAL_WINDOW table is not the same as nrow in "
        "DATA_DESCRIPTION table");
  unsigned int spw = reader_.spectralWindow();
  // Remove all rows before and after the selected band.
  // Do it from the end, otherwise row numbers change.
  for (unsigned int i = out_spw.nrow(); i > 0;) {
    if (--i != spw) {
      out_spw.removeRow(i);
      out_dd.removeRow(i);
    }
  }
  // Set nr of channels.
  ScalarColumn<int> channum(out_spw, "NUM_CHAN");
  channum.fillColumn(nr_chan_);
  // Change the column shapes.
  TableDesc tdesc = in_spw.tableDesc();
  MakeArrayColumn(tdesc["CHAN_FREQ"], shape, 0, out_spw);
  MakeArrayColumn(tdesc["CHAN_WIDTH"], shape, 0, out_spw);
  MakeArrayColumn(tdesc["EFFECTIVE_BW"], shape, 0, out_spw);
  MakeArrayColumn(tdesc["RESOLUTION"], shape, 0, out_spw);
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

void MSWriter::UpdateObs(const string& out_name) {
  Table out_obs = Table(out_name + "/OBSERVATION", Table::Update);
  // Set nr of channels.
  ArrayColumn<double> time_range(out_obs, "TIME_RANGE");
  casacore::Vector<double> times(2);
  times[0] = reader_.firstTime() - 0.5 * reader_.getInfo().timeInterval();
  times[1] = reader_.lastTime() + 0.5 * reader_.getInfo().timeInterval();
  // There should be one row, but loop in case of.
  for (unsigned int i = 0; i < out_obs.nrow(); ++i) {
    time_range.put(i, times);
  }
}

void MSWriter::UpdatePhaseCentre(const string& out_name, const DPInfo& info) {
  Table out_field = Table(out_name + "/FIELD", Table::Update);
  // Write new phase center.
  ArrayMeasColumn<MDirection> phase_col(out_field, "PHASE_DIR");
  // If a moving reference type like AZELGEO was used in the original MS, and
  // the phase centre is changed (with a phaseshift), the ref frame of the
  // column must be reset:
  phase_col.setDescRefCode(info.phaseCenter().getRefPtr()->getType(), false);
  casacore::Vector<MDirection> dir(1, info.phaseCenter());
  phase_col.put(0, dir);
}

void MSWriter::UpdateBeam(const std::string& out_name,
                          const std::string& out_col_name, const DPInfo& info) {
  const char *beam_mode_field_name = "LOFAR_APPLIED_BEAM_MODE",
             *beam_dir_field_name = "LOFAR_APPLIED_BEAM_DIR";

  Table main_table(out_name, Table::Update);
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

void MSWriter::UpdateInternalBuffer(const base::DPBuffer& buffer) {
  const common::NSTimer::StartStop timer(update_buffer_timer_);

  internal_buffer_.referenceFilled(buffer);
  reader_.fetchWeights(buffer, internal_buffer_, timer_);
  reader_.fetchUVW(buffer, internal_buffer_, timer_);
  if (write_full_res_flags_) {
    reader_.fetchFullResFlags(buffer, internal_buffer_, timer_);
  }
}

void MSWriter::WriteData(Table& out, const DPBuffer& buf) {
  ArrayColumn<casacore::Complex> data_col(out, data_col_name_);
  ArrayColumn<bool> flag_col(out, "FLAG");
  ScalarColumn<bool> flag_row_col(out, "FLAG_ROW");

  if (buf.getData().empty()) {
    return;
  }

  // Write WEIGHT_SPECTRUM and DATA
  ArrayColumn<float> weightCol(out, "WEIGHT_SPECTRUM");
  const Array<float>& weights = buf.getWeights();

  // If compressing, flagged values need to be set to NaN, and flagged
  // weights to zero, to decrease the dynamic range
  if (st_man_keys_.stManName == "dysco") {
    Cube<casacore::Complex> data_copy = buf.getData().copy();
    Cube<casacore::Complex>::iterator data_iter = data_copy.begin();
    Cube<float> weights_copy = weights.copy();
    Cube<float>::iterator weights_iter = weights_copy.begin();
    for (Cube<bool>::const_iterator flag_iter = buf.getFlags().begin();
         flag_iter != buf.getFlags().end(); ++flag_iter) {
      if (*flag_iter) {
        *data_iter = casacore::Complex(std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN());
        *weights_iter = 0.;
      }
      ++data_iter;
      ++weights_iter;
    }
    data_col.putColumn(data_copy);
    weightCol.putColumn(weights_copy);
  } else {
    data_col.putColumn(buf.getData());
    weightCol.putColumn(weights);
  }

  flag_col.putColumn(buf.getFlags());
  // A row is flagged if no flags in the row are False.
  auto c = partialNFalse(buf.getFlags(), IPosition(2, 0, 1));
  casacore::Vector<bool> row_flags(c == decltype(c)::value_type(0));
  flag_row_col.putColumn(row_flags);
  if (write_full_res_flags_) {
    WriteFullResFlags(out, buf);
  }

  // Write UVW
  ArrayColumn<double> uvw_col(out, "UVW");
  const Array<double>& uvws = reader_.fetchUVW(buf, internal_buffer_, timer_);
  uvw_col.putColumn(uvws);
}

void MSWriter::WriteFullResFlags(Table& out, const DPBuffer& buf) {
  const Cube<bool>& flags = buf.getFullResFlags();
  const IPosition& of_shape = flags.shape();
  if ((unsigned int)(of_shape[0]) != n_chan_avg_ * nr_chan_)
    throw std::runtime_error(
        "Full Res Flags size " + std::to_string(of_shape[0]) +
        " does not equal " + std::to_string(n_chan_avg_) + "*" +
        std::to_string(nr_chan_) +
        ".\nTry setting \"msout.writefullresflag=false\" in input parset");
  if ((unsigned int)(of_shape[1]) != n_time_avg_)
    throw std::runtime_error(std::to_string(of_shape[1]) +
                             std::to_string(n_time_avg_));
  // Convert the bools to unsigned char bits.
  IPosition ch_shape(of_shape);
  ch_shape[0] = (of_shape[0] + 7) / 8;
  Cube<unsigned char> chars(ch_shape);
  // If their sizes match, do it all in one go.
  // Otherwise we have to iterate.
  if (of_shape[0] == ch_shape[0] * 8) {
    casacore::Conversion::boolToBit(chars.data(), flags.data(), flags.size());
  } else {
    if (of_shape[0] > ch_shape[0] * 8)
      throw std::runtime_error("Incorrect shape of full res flags");
    const bool* flags_ptr = flags.data();
    unsigned char* chars_ptr = chars.data();
    for (int i = 0; i < of_shape[1] * of_shape[2]; ++i) {
      casacore::Conversion::boolToBit(chars_ptr, flags_ptr, of_shape[0]);
      flags_ptr += of_shape[0];
      chars_ptr += ch_shape[0];
    }
  }
  ArrayColumn<unsigned char> full_res_col(out, "LOFAR_FULL_RES_FLAG");
  if (!full_res_col.keywordSet().isDefined("NCHAN_AVG")) {
    full_res_col.rwKeywordSet().define("NCHAN_AVG", int(n_chan_avg_));
    full_res_col.rwKeywordSet().define("NTIME_AVG", int(n_time_avg_));
  }
  full_res_col.putColumn(chars);
}

void MSWriter::WriteMeta(Table& out, const DPBuffer& buf) {
  // Fill ANTENNA1/2.
  ScalarColumn<int> ant1col(out, "ANTENNA1");
  ScalarColumn<int> ant2col(out, "ANTENNA2");
  ant1col.putColumn(casacore::Vector<int>(getInfo().getAnt1()));
  ant2col.putColumn(casacore::Vector<int>(getInfo().getAnt2()));
  // Fill all rows that do not change.
  FillSca<double>(buf.getTime(), out, "TIME");
  FillSca<double>(buf.getTime(), out, "TIME_CENTROID");
  FillSca<double>(buf.getExposure(), out, "EXPOSURE");
  FillSca<double>(interval_, out, "INTERVAL");
  FillSca<int>(0, out, "FEED1");
  FillSca<int>(0, out, "FEED2");
  FillSca<int>(0, out, "DATA_DESC_ID");
  FillSca<int>(0, out, "PROCESSOR_ID");
  FillSca<int>(0, out, "FIELD_ID");
  FillSca<int>(0, out, "SCAN_NUMBER");
  FillSca<int>(0, out, "ARRAY_ID");
  FillSca<int>(0, out, "OBSERVATION_ID");
  FillSca<int>(0, out, "STATE_ID");
  Array<float> arr(IPosition(1, nr_corr_));
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
  base::DPBuffer buffer;
  while (write_queue_.read(buffer)) {
    ProcessBuffer(buffer);
  }
}

void MSWriter::CreateTask() {
  const common::NSTimer::StartStop timer(create_task_timer_);

  base::DPBuffer buffer;
  buffer.copy(internal_buffer_);
  write_queue_.write(std::move(buffer));
}

}  // namespace steps
}  // namespace dp3
