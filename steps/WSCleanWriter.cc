// WSCleanWriter.cc: DP3 step writing a reordered MS
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "WSCleanWriter.h"

#include <map>
#include <stdexcept>
#include <algorithm>

#include <schaapcommon/reordering/msselection.h>
#include <schaapcommon/reordering/handledata.h>
#include <schaapcommon/reordering/filewriter.h>

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>
#include <aocommon/multibanddata.h>

#include "../common/ParameterSet.h"

using aocommon::Logger;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using schaapcommon::reordering::FileWriter;

namespace dp3 {

namespace steps {

WSCleanWriter::WSCleanWriter(const common::ParameterSet& parset,
                             const std::string& prefix)
    : name_(prefix),
      parset_(parset),
      temporary_directory_(
          parset_.getString(prefix + "temporary_directory", "")),
      channels_per_file_(parset_.getUint32(prefix + "chanperfile", 0)) {
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
}

WSCleanWriter::~WSCleanWriter() = default;

aocommon::BandData WSCleanWriter::GetBand(size_t start_channel,
                                          size_t end_channel) const {
  std::vector<aocommon::ChannelInfo> channels;
  for (size_t ch = start_channel; ch != end_channel; ++ch)
    channels.emplace_back(getInfoOut().chanFreqs()[ch], 0.0);
  return aocommon::BandData(channels, channels[0].Frequency());
}

void WSCleanWriter::StartReorder() {
  data_desc_id_ = getInfoOut().spectralWindow();

  std::vector<schaapcommon::reordering::ChannelRange> channel_ranges;

  std::vector<aocommon::MultiBandData> bands_per_part;
  const size_t start_channel = getInfoOut().startchan();
  const size_t end_channel = getInfoOut().startchan() + getInfoOut().nchan();
  if (channels_per_file_ == 0) {
    channel_ranges.push_back({data_desc_id_, start_channel, end_channel});

    aocommon::MultiBandData& bands = bands_per_part.emplace_back();
    bands.SetBand(data_desc_id_, GetBand(start_channel, end_channel));
  } else {
    // TODO when the nr of channels are not integer divisable, the output
    // channels can have widely different nr of channels, which is not good for
    // image quality. The best would be to use WSClean's code for dividing the
    // channels.
    size_t part_start_channel = start_channel;
    for (; part_start_channel + channels_per_file_ < end_channel;
         part_start_channel += channels_per_file_) {
      const size_t part_end_channel = part_start_channel + channels_per_file_;
      channel_ranges.push_back(
          {data_desc_id_, part_start_channel, part_end_channel});
      aocommon::MultiBandData& bands = bands_per_part.emplace_back();
      bands.SetBand(data_desc_id_,
                    GetBand(part_start_channel, part_end_channel));
    }

    if (part_start_channel < end_channel) {
      channel_ranges.push_back(
          {data_desc_id_, part_start_channel, end_channel});
      aocommon::MultiBandData& bands = bands_per_part.emplace_back();
      bands.SetBand(data_desc_id_, GetBand(part_start_channel, end_channel));
    }
  }

  Logger::Info << "Dividing channels into " << channel_ranges.size()
               << " files\n";

  schaapcommon::reordering::MSSelection selection;

  handle_data_ = std::make_unique<schaapcommon::reordering::HandleData>(
      out_name_, getInfoOut().dataColumnName(), "MODEL_DATA",
      schaapcommon::reordering::StorageManagerType::Default,
      temporary_directory_,
      schaapcommon::reordering::MakeRegularChannelMap(channel_ranges), false,
      false, pols_out_, selection, bands_per_part,
      getInfoOut().antennaNames().size(), true,
      [](schaapcommon::reordering::HandleData) {});

  std::vector<aocommon::OptionalNumber<size_t>> data_desc_ids;
  std::tie(handle_data_->metadata_indices_, data_desc_ids) =
      schaapcommon::reordering::MakeMetaFilesMap(handle_data_->channels_);
  std::map<size_t, std::set<aocommon::PolarizationEnum>> pol_per_data_desc_id{
      {data_desc_id_, getInfoOut().polarizations()}};

  writer_ = std::make_unique<FileWriter>(*handle_data_, pol_per_data_desc_id,
                                         data_desc_ids,
                                         getInfoOut().startTime() / 86400);
}

bool WSCleanWriter::process(std::unique_ptr<dp3::base::DPBuffer> buffer) {
  const common::NSTimer::StartStop sstime(timer_);

  ReorderBuffer(*buffer);
  getNextStep()->process(std::move(buffer));
  return true;
}

void WSCleanWriter::ReorderBuffer(dp3::base::DPBuffer& buffer) {
  const common::NSTimer::StartStop timer(writer_timer_);

  const dp3::base::DPBuffer::FlagsType& buff_flags = buffer.GetFlags();
  const dp3::base::DPBuffer::UvwType& buff_uvw = buffer.GetUvw();
  const dp3::base::DPBuffer::WeightsType& buff_weights = buffer.GetWeights();
  const dp3::base::DPBuffer::DataType& buff_data = buffer.GetData();

  const size_t n_baselines = buff_data.shape(0);

  const std::vector<int>& antenna1_list = getInfoOut().getAnt1();
  const std::vector<int>& antenna2_list = getInfoOut().getAnt2();
  const std::set<aocommon::PolarizationEnum> pols_in =
      getInfoOut().polarizations();

  for (size_t bl = 0; bl < n_baselines; bl++) {
    // Skip self-correlations, in WSClean this is done using an MSSelection
    if (antenna1_list[bl] == antenna2_list[bl]) continue;

    const bool* flag_ptr = &buff_flags(bl, 0, 0);
    const float* weight_ptr = &buff_weights(bl, 0, 0);
    const std::complex<float>* data_ptr = &buff_data(bl, 0, 0);

    writer_->WriteMetaRow(buff_uvw(bl, 0), buff_uvw(bl, 1), buff_uvw(bl, 2),
                          buffer.GetTime(), data_desc_id_, antenna1_list[bl],
                          antenna2_list[bl], 0);

    writer_->WriteDataRow(data_ptr, nullptr, weight_ptr, flag_ptr,
                          data_desc_id_);
  }
}

void WSCleanWriter::FinishReorder() {
  writer_->UpdateMetaHeaders();
  writer_->UpdatePartHeaders(false);
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
  if (getInfoOut().metaChanged()) {
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
