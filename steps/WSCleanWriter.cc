// WSCleanWriter.cc: DP3 step writing a reordered MS
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MSReorder.h"
#include "WSCleanWriter.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <limits>
#include <fstream>

#include "../common/ParameterSet.h"

#include <aocommon/logger.h>
#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <boost/filesystem/path.hpp>
#include <memory>
#include <stdexcept>

using aocommon::Logger;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {

namespace steps {

WSCleanWriter::WSCleanWriter(const common::ParameterSet& parset,
                             const std::string& prefix)
    : name_(prefix),
      parset_(parset),
      temporary_directory_(
          parset_.getString(prefix + "temporaryDirectory", "")),
      selected_row_count_(0) {
  out_name_ = parset_.getString(prefix + "name", "");
  if (out_name_.empty() && parset.isDefined("msout.name"))
    out_name_ = parset_.getString("msout.name", "");
  if (out_name_.empty() && parset.isDefined("msout"))
    out_name_ = parset_.getString("msout", "");
  if (out_name_.empty() || out_name_ == ".")
    out_name_ = parset_.getString("msin", "");
  if (out_name_.empty())
    throw std::runtime_error("Could not infer name of the reordered files");

  const std::string pols_string =
      parset_.getString(prefix + "polarization", "I");
  if (pols_string == "instr") {
    pols_out_ = std::set{aocommon::Polarization::Instrumental};
  } else if (pols_string == "diag_instr") {
    pols_out_ = std::set{aocommon::Polarization::DiagonalInstrumental};
  } else {
    pols_out_ = aocommon::Polarization::ParseList(pols_string);
  }
  nr_polarizations_ = pols_out_.size();
}

WSCleanWriter::~WSCleanWriter() = default;

void WSCleanWriter::StartReorder() {
  files_.resize(nr_polarizations_);

  start_time_ = dp_info_.startTime() / 86400;
  channel_start_ = dp_info_.startchan();
  channel_count_ = dp_info_.nchan();
  data_desc_id_ = dp_info_.spectralWindow();

  size_t file_index = 0;
  size_t part = 0;
  for (aocommon::PolarizationEnum pol : pols_out_) {
    ReorderFile& file = files_[file_index];
    const std::string part_prefix = reorder::GetPartPrefix(
        out_name_, part, pol, data_desc_id_, temporary_directory_);
    file.data = std::make_unique<std::ofstream>(part_prefix + ".tmp");
    file.weight = std::make_unique<std::ofstream>(part_prefix + "-w.tmp");
    file.data->seekp(reorder::PartHeader::BINARY_SIZE, std::ios::beg);
    ++file_index;
  }

  const std::string meta_filename =
      reorder::GetMetaFilename(out_name_, temporary_directory_, data_desc_id_);
  meta_file_ptr_ = std::make_unique<std::ofstream>(meta_filename);
  reorder::MetaHeader meta_header;
  meta_header.selected_row_count = 0;  // not yet known
  meta_header.filename_length = out_name_.size();
  meta_header.start_time = 0;
  meta_header.Write(*meta_file_ptr_);
  meta_file_ptr_->write(out_name_.c_str(), out_name_.size());
  if (!meta_file_ptr_->good())
    throw std::runtime_error("Error writing to temporary file " +
                             meta_filename);
}

bool WSCleanWriter::process(std::unique_ptr<dp3::base::DPBuffer> buffer) {
  const common::NSTimer::StartStop sstime(timer_);

  ReorderBuffer(*buffer);
  getNextStep()->process(std::move(buffer));
  return true;
}

void WSCleanWriter::ReorderBuffer(dp3::base::DPBuffer& buffer) {
  const common::NSTimer::StartStop timer(writer_timer_);

  // Get DP3 time frame details.

  const dp3::base::DPBuffer::FlagsType& buff_flags = buffer.GetFlags();
  const dp3::base::DPBuffer::UvwType& buff_uvw = buffer.GetUvw();
  const dp3::base::DPBuffer::WeightsType& buff_weights = buffer.GetWeights();
  const dp3::base::DPBuffer::DataType& buff_data = buffer.GetData();

  const size_t n_baselines = buff_data.shape(0);
  const size_t n_channels = buff_data.shape(1);

  const std::vector<int>& antenna1_list = dp_info_.getAnt1();
  const std::vector<int>& antenna2_list = dp_info_.getAnt2();
  const std::set<aocommon::PolarizationEnum> pols_in = dp_info_.polarizations();

  const size_t polarizations_per_file =
      aocommon::Polarization::GetVisibilityCount(*pols_out_.begin());

  std::vector<std::complex<float>> data_buffer_(
      n_channels * polarizations_per_file, 0.0);
  std::vector<float> weight_buffer_(n_channels * polarizations_per_file);

  for (size_t bl = 0; bl < n_baselines; bl++) {
    // Skip self-correlations, in WSClean this is done using an MSSelection
    if (antenna1_list[bl] == antenna2_list[bl]) continue;

    reorder::MetaRecord meta;
    const bool* flag_ptr = &buff_flags(bl, 0, 0);
    const float* weight_ptr = &buff_weights(bl, 0, 0);
    const std::complex<float>* data_ptr = &buff_data(bl, 0, 0);

    meta.u = buff_uvw(bl, 0);
    meta.v = buff_uvw(bl, 1);
    meta.w = buff_uvw(bl, 2);
    meta.antenna1 = antenna1_list[bl];
    meta.antenna2 = antenna2_list[bl];
    meta.field_id = 0;
    meta.time = buffer.GetTime();

    ++selected_row_count_;

    meta.Write(*meta_file_ptr_);
    if (!meta_file_ptr_->good())
      throw std::runtime_error("Error writing to temporary file");

    size_t file_index = 0;
    for (aocommon::PolarizationEnum pol : pols_out_) {
      ReorderFile& file = files_[file_index];

      reorder::ExtractData(data_buffer_.data(), 0, n_channels, pols_in,
                           data_ptr, pol);
      file.data->write(
          reinterpret_cast<char*>(data_buffer_.data()),
          n_channels * polarizations_per_file * sizeof(std::complex<float>));
      if (!file.data->good())
        throw std::runtime_error("Error writing to temporary data file");

      reorder::ExtractWeights(weight_buffer_.data(), 0, n_channels, pols_in,
                              data_ptr, weight_ptr, flag_ptr, pol);
      file.weight->write(reinterpret_cast<char*>(weight_buffer_.data()),
                         n_channels * polarizations_per_file * sizeof(float));
      if (!file.weight->good())
        throw std::runtime_error("Error writing to temporary weights file");

      ++file_index;
    }
  }
}

void WSCleanWriter::FinishReorder() {
  reorder::MetaHeader meta_header;
  meta_header.selected_row_count = selected_row_count_;
  meta_header.filename_length = out_name_.size();
  meta_header.start_time = start_time_;
  meta_file_ptr_->seekp(0);
  meta_header.Write(*meta_file_ptr_);
  meta_file_ptr_->write(out_name_.c_str(), out_name_.size());

  reorder::PartHeader header;
  header.has_model = false;
  header.channel_start = channel_start_;
  header.channel_count = channel_count_;
  header.data_desc_id = data_desc_id_;
  for (size_t file_index = 0; file_index < pols_out_.size(); file_index++) {
    ReorderFile& file = files_[file_index];
    file.data->seekp(0, std::ios::beg);
    header.Write(*file.data);
    if (!file.data->good())
      throw std::runtime_error("Error writing to temporary data file");

    file.data.reset();
    file.weight.reset();
  }
}

void WSCleanWriter::finish() {
  FinishReorder();
  if (getNextStep()) getNextStep()->finish();
}

void WSCleanWriter::show(std::ostream& os) const {
  os << "WSCleanWriter \n";
  os << "Reordering and writing to " << out_name_ << "\n";
  os << "Reordering to polarizations ";
  for (aocommon::PolarizationEnum pol : pols_out_) {
    os << aocommon::Polarization::TypeToFullString(pol) << ' ';
  }
  os << '\n';
}

void WSCleanWriter::updateInfo(const base::DPInfo& info_in) {
  OutputStep::updateInfo(info_in);
  dp_info_ = info_in;
  if (getInfo().metaChanged()) {
    Logger::Warn
        << "Meta data changes detected (by averaging, adding/removing"
        << " stations, phase shifting or upsampling) when reordering."
        << " Ensure that there is an output step to write out a Measurement"
        << " set to prevent data inconsistency.\n";
  }
  StartReorder();
}

void WSCleanWriter::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " WSCleanWriter " << name_ << '\n';

  duration = timer_.getElapsed();
  os << "    ";
  FlagCounter::showPerc1(os, writer_timer_.getElapsed(), duration);
  os << " Writing\n";
}

}  // namespace steps
}  // namespace dp3
