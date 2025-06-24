// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SVPInput.h"
#include <base/SubtableWriter.h>

#include <boost/filesystem.hpp>

#include <sys/socket.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#include <sys/un.h>

#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>

#include <aocommon/logger.h>

using aocommon::Logger;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {
SVPInput::SVPInput([[maybe_unused]] const common::ParameterSet& parset,
                   const std::string& prefix)
    : name_(prefix),
      receiver_socket_(parset.getString(prefix + "socket", kSocketName)) {
  ConnectServer(receiver_socket_.c_str());
  InitializeInfo();
  CreateInitialSubtables();
}

void SVPInput::finish() {
  DisconnectServer();
  getNextStep()->finish();
}

bool SVPInput::process(std::unique_ptr<dp3::base::DPBuffer> buffer) {
  const size_t n_bl = getInfoOut().nbaselines();
  const size_t n_ch = getInfoOut().nchan();
  const size_t n_cr = getInfoOut().ncorr();

  buffer->GetData().resize({n_bl, n_ch, n_cr});

  int recv;
  double time_centroid = 0.0, time_duration = 1.0;
  recv = ReceiveItem<double>(time_centroid);
  recv = ReceiveItem<double>(time_duration);
  // Wait till data arrives
  size_t data_size;
  recv = ReceiveItem<size_t>(data_size);
  uint data_bits;
  recv = ReceiveItem<uint>(data_bits);

  if (recv == -1) {
    aocommon::Logger::Debug
        << "DP3 SVPInput Error receiving data size, skipping";
    return false;
  }
  assert(data_bits == 16 or data_bits == 32);

  float inv_scale = (float)1.0 / metadata_.scale_factor_;

  if (data_bits == 16) {
    std::vector<int16_t> recv_buffer(2 * n_bl * n_ch * n_cr);
    recv = ReceiveBytes(reinterpret_cast<char*>(recv_buffer.data()),
                        2 * n_bl * n_ch * n_cr * sizeof(int16_t));
    if (recv == -1) {
      aocommon::Logger::Debug << "DP3 SVPInput Error receiving data buffer";
    }

    for (size_t bl = 0; bl < n_bl; bl++) {
      for (size_t ch = 0; ch < n_ch; ch++) {
        std::complex<float>* data_pointer = &buffer->GetData()(bl, ch, 0);

        size_t recv_index = bl * n_ch * n_cr + ch * n_cr + 0;
        int16_t recv_r = recv_buffer[2 * recv_index];
        int16_t recv_i = recv_buffer[2 * recv_index + 1];

        data_pointer[0] = std::complex<float>((float)recv_r * inv_scale,
                                              (float)recv_i * inv_scale);

        recv_index = bl * n_ch * n_cr + ch * n_cr + 1;
        recv_r = recv_buffer[2 * recv_index];
        recv_i = recv_buffer[2 * recv_index + 1];

        data_pointer[1] = std::complex<float>((float)recv_r * inv_scale,
                                              (float)recv_i * inv_scale);

        if (n_cr == 4) {
          recv_index = bl * n_ch * n_cr + ch * n_cr + 2;
          recv_r = recv_buffer[2 * recv_index];
          recv_i = recv_buffer[2 * recv_index + 1];

          data_pointer[2] = std::complex<float>((float)recv_r * inv_scale,
                                                (float)recv_i * inv_scale);

          recv_index = bl * n_ch * n_cr + ch * n_cr + 3;
          recv_r = recv_buffer[2 * recv_index];
          recv_i = recv_buffer[2 * recv_index + 1];

          data_pointer[3] = std::complex<float>((float)recv_r * inv_scale,
                                                (float)recv_i * inv_scale);
        }
      }
    }
  } else if (data_bits == 32) {
    std::vector<int32_t> recv_buffer(2 * n_bl * n_ch * n_cr);
    recv = ReceiveBytes(reinterpret_cast<char*>(recv_buffer.data()),
                        2 * n_bl * n_ch * n_cr * sizeof(int32_t));
    if (recv == -1) {
      aocommon::Logger::Debug << "DP3 SVPInput Error receiving data buffer";
    }

    for (size_t bl = 0; bl < n_bl; bl++) {
      for (size_t ch = 0; ch < n_ch; ch++) {
        std::complex<float>* data_pointer = &buffer->GetData()(bl, ch, 0);
        size_t recv_index = bl * n_ch * n_cr + ch * n_cr + 0;
        int32_t recv_r = recv_buffer[2 * recv_index];
        int32_t recv_i = recv_buffer[2 * recv_index + 1];

        data_pointer[0] = std::complex<float>((float)recv_r * inv_scale,
                                              (float)recv_i * inv_scale);

        recv_index = bl * n_ch * n_cr + ch * n_cr + 1;
        recv_r = recv_buffer[2 * recv_index];
        recv_i = recv_buffer[2 * recv_index + 1];

        data_pointer[1] = std::complex<float>((float)recv_r * inv_scale,
                                              (float)recv_i * inv_scale);

        if (n_cr == 4) {
          recv_index = bl * n_ch * n_cr + ch * n_cr + 2;
          recv_r = recv_buffer[2 * recv_index];
          recv_i = recv_buffer[2 * recv_index + 1];

          data_pointer[2] = std::complex<float>((float)recv_r * inv_scale,
                                                (float)recv_i * inv_scale);

          recv_index = bl * n_ch * n_cr + ch * n_cr + 3;
          recv_r = recv_buffer[2 * recv_index];
          recv_i = recv_buffer[2 * recv_index + 1];

          data_pointer[3] = std::complex<float>((float)recv_r * inv_scale,
                                                (float)recv_i * inv_scale);
        }
      }
    }
  }

  buffer->SetTime(time_centroid);
  buffer->SetExposure(time_duration);
  // Add aux data
  buffer->GetWeights().resize({n_bl, n_ch, n_cr});
  buffer->GetFlags().resize({n_bl, n_ch, n_cr});
  buffer->GetWeights().fill(1.0);
  buffer->GetFlags().fill(false);

  // calculate uvw
  buffer->GetUvw().resize({n_bl, 3});
  const std::vector<int>& antenna1 = getInfoOut().getAnt1();
  const std::vector<int>& antenna2 = getInfoOut().getAnt2();
  for (size_t bl = 0; bl < n_bl; bl++) {
    xt::view(buffer->GetUvw(), bl, xt::all()) = xt::adapt(
        uvw_calculator_->getUVW(antenna2[bl], antenna1[bl], time_centroid));
  }

  getNextStep()->process(std::move(buffer));
  return true;
}

void SVPInput::show(std::ostream& os) const {
  os << "SVPInput " << std::endl;
  os << "Temporary MS :" << temp_out_ms_ << std::endl;
  os << "Socket :" << receiver_socket_ << std::endl;
}

std::string SVPInput::msName() const { return temp_out_ms_; }

void SVPInput::InitializeInfo() {
  // Get (static) metadata
  // channels
  int recv = ReceiveItem<size_t>(metadata_.nr_channels_);
  // antennas
  recv = ReceiveItem<size_t>(metadata_.nr_antennas_);
  // polarizations
  recv = ReceiveItem<int>(metadata_.nr_polarizations_);
  // cross corr scale factor
  recv = ReceiveItem<double>(metadata_.scale_factor_);
  // start time
  recv = ReceiveItem<double>(metadata_.start_time_);
  // time interval (integration time)
  recv = ReceiveItem<double>(metadata_.integration_time_);
  // number of integrations
  recv = ReceiveItem<size_t>(metadata_.nr_times_);
  // source name (length and name)
  int name_len = 0;
  recv = ReceiveItem<int>(name_len);
  if (name_len != 0) {
    metadata_.source_name_.resize(name_len);
    recv = ReceiveBytes(reinterpret_cast<char*>(metadata_.source_name_.data()),
                        name_len * sizeof(char));
  } else {
    metadata_.source_name_ = "TEST";
    aocommon::Logger::Debug
        << "DP3 SVPInput did not receive source name, using default "
        << metadata_.source_name_;
  }
  // source coords (phase center, also add proper motion)
  recv = ReceiveItem<double>(metadata_.ra_);
  recv = ReceiveItem<double>(metadata_.dec_);
  // frame name (length and name)
  name_len = 0;
  recv = ReceiveItem<int>(name_len);
  if (name_len != 0) {
    metadata_.frame_name_.resize(name_len);
    recv = ReceiveBytes(reinterpret_cast<char*>(metadata_.frame_name_.data()),
                        name_len * sizeof(char));
  } else {
    metadata_.frame_name_ = "J2000";
  }
  // telecope name (length and name)
  name_len = 0;
  recv = ReceiveItem<int>(name_len);
  if (name_len != 0) {
    metadata_.telescope_name_.resize(name_len);
    recv =
        ReceiveBytes(reinterpret_cast<char*>(metadata_.telescope_name_.data()),
                     name_len * sizeof(char));
  } else {
    metadata_.telescope_name_ = "ALMA";
  }

  // antenna coords
  std::vector<double> ant_pos(metadata_.nr_antennas_ * 3);
  recv = ReceiveBytes(reinterpret_cast<char*>(ant_pos.data()),
                      3 * metadata_.nr_antennas_ * sizeof(double));

  // antenna diameters
  metadata_.antenna_diameters_.resize(metadata_.nr_antennas_);
  recv =
      ReceiveBytes(reinterpret_cast<char*>(metadata_.antenna_diameters_.data()),
                   metadata_.nr_antennas_ * sizeof(double));

  // channel frequencies
  metadata_.chan_freqs_.resize(metadata_.nr_channels_);
  metadata_.chan_widths_.resize(metadata_.nr_channels_);
  recv = ReceiveBytes(reinterpret_cast<char*>(metadata_.chan_freqs_.data()),
                      metadata_.nr_channels_ * sizeof(double));
  recv = ReceiveBytes(reinterpret_cast<char*>(metadata_.chan_widths_.data()),
                      metadata_.nr_channels_ * sizeof(double));

  if (recv == -1) {
    throw std::runtime_error("DP3 SVPInput did not receive metadata.");
  }

  metadata_.antenna_positions_.resize(metadata_.nr_antennas_);
  for (size_t ci = 0; ci < metadata_.nr_antennas_; ci++) {
    casacore::Vector<double> vals(3);
    vals[0] = ant_pos[3 * ci];
    vals[1] = ant_pos[3 * ci + 1];
    vals[2] = ant_pos[3 * ci + 2];
    metadata_.antenna_positions_[ci] = casacore::MPosition(
        casacore::Quantum<casacore::Vector<double>>(vals, "m"),
        casacore::MPosition::ITRF);
  }
  // baselines, exclude auto-corr
  size_t nbl = metadata_.nr_antennas_ * (metadata_.nr_antennas_ - 1) / 2;
  // +metadata_.nr_antennas_; for autocorr
  std::vector<int> ant1(nbl);
  std::vector<int> ant2(nbl);
  std::size_t bl = 0;
  for (std::size_t st2 = 0; st2 < metadata_.nr_antennas_; st2++) {
    for (std::size_t st1 = 0; st1 < st2; st1++) {
      ant1[bl] = st1;
      ant2[bl] = st2;
      bl++;
    }
  }

  // Set antenna names
  metadata_.antenna_names_.resize(metadata_.nr_antennas_);
  for (size_t ci = 0; ci < metadata_.nr_antennas_; ci++) {
    metadata_.antenna_names_[ci] =
        metadata_.telescope_name_ + std::to_string(ci);
  }

  // Set integration time and number of time steps
  metadata_.end_time_ = metadata_.start_time_ + metadata_.integration_time_ *
                                                    (double)metadata_.nr_times_;

  casacore::MDirection::Types frame;
  if (metadata_.frame_name_ == "ICRS") {
    frame = casacore::MDirection::ICRS;
  } else {
    frame = casacore::MDirection::J2000;
  }

  phase_direction_ =
      casacore::MDirection(casacore::Quantity(metadata_.ra_, "rad"),
                           casacore::Quantity(metadata_.dec_, "rad"), frame);

  DPInfo dataset_info(metadata_.nr_polarizations_, metadata_.nr_channels_);

  dataset_info.setAntennas(metadata_.antenna_names_,
                           metadata_.antenna_diameters_,
                           metadata_.antenna_positions_, ant1, ant2);
  dataset_info.setChannels(std::vector(metadata_.chan_freqs_),
                           std::vector(metadata_.chan_widths_));

  dataset_info.setArrayInformation(metadata_.antenna_positions_[0],
                                   phase_direction_, phase_direction_,
                                   phase_direction_);

  uvw_calculator_ = std::make_unique<base::UVWCalculator>(
      phase_direction_, metadata_.antenna_positions_[0],
      metadata_.antenna_positions_);

  dataset_info.setTimeIntervalAndSteps(metadata_.integration_time_,
                                       metadata_.nr_times_);

  GetWritableInfoOut() = dataset_info;
}

void SVPInput::CreateInitialSubtables() {
  temp_out_ms_ = (boost::filesystem::temp_directory_path() /
                  boost::filesystem::unique_path("TEMP%%%%%%%.MS"))
                     .string();

  GetWritableInfoOut().setMsNames(temp_out_ms_, kDataColumnName,
                                  kFlagColumnName, kWeightColumnName);
  std::unique_ptr<base::SubtableWriter> writer =
      std::make_unique<base::SubtableWriter>(temp_out_ms_,
                                             metadata_.nr_channels_);

  CreateAntennaTable(metadata_.start_time_, *writer);

  CreateSpectralWindowTable(metadata_.nr_channels_, metadata_.chan_freqs_,
                            metadata_.chan_widths_, *writer);
  CreateSourceTable(metadata_.source_name_, metadata_.start_time_,
                    metadata_.end_time_, *writer);
  CreateFieldTable(metadata_.source_name_, metadata_.start_time_, *writer);
  CreateObservationTable(metadata_.start_time_, metadata_.end_time_, *writer);
  CreatePolarizationTable(*writer, metadata_.nr_polarizations_);
}

void SVPInput::CreateAntennaTable(double time, base::SubtableWriter& writer) {
  std::vector<base::SubtableWriter::AntennaInfo> antennas(
      metadata_.nr_antennas_);
  for (size_t ant = 0; ant != metadata_.nr_antennas_; ++ant) {
    base::SubtableWriter::AntennaInfo& antenna_info = antennas[ant];
    antenna_info.name = metadata_.antenna_names_[ant];
    antenna_info.station = metadata_.telescope_name_;
    antenna_info.type = "GROUND-BASED";
    antenna_info.mount = "ALT-AZ";
    antenna_info.x = metadata_.antenna_positions_[ant].getValue()(0);
    antenna_info.y = metadata_.antenna_positions_[ant].getValue()(1);
    antenna_info.z = metadata_.antenna_positions_[ant].getValue()(2);
    antenna_info.diameter = metadata_.antenna_diameters_[ant];
    antenna_info.flag = false;
  }

  writer.WriteAntennas(antennas, metadata_.axes_, time);
}

void SVPInput::CreateSpectralWindowTable(int nchannels,
                                         std::vector<double> const& chan_freqs,
                                         std::vector<double> const& chan_widths,
                                         base::SubtableWriter& writer) {
  std::vector<base::SubtableWriter::ChannelInfo> channels(nchannels);
  for (size_t ch = 0; ch != channels.size(); ++ch) {
    base::SubtableWriter::ChannelInfo& channel = channels[ch];
    channel.channel_frequency = chan_freqs[ch];
    channel.channel_width = chan_widths[ch];
    channel.effective_bandwidth = chan_widths[ch];
    channel.resolution = chan_widths[ch];
  }
  std::string band_name = metadata_.telescope_name_ + std::string("BAND");
  writer.WriteBandInfo(band_name, channels, chan_freqs[0], chan_widths[0],
                       false);
}

void SVPInput::CreateSourceTable(std::string& name, double startTime,
                                 double endTime, base::SubtableWriter& writer) {
  base::SubtableWriter::SourceInfo source;

  source.source_id = 0;
  source.time = startTime;
  source.interval = endTime;
  source.spectral_window_id = 0;
  source.num_lines = 0;
  source.name = name;
  source.calibration_group = 0;
  source.code = "";
  source.ra = metadata_.ra_;  // (in radians)
  source.dec = metadata_.dec_;
  source.proper_motion[0] = 0.0;
  source.proper_motion[1] = 0.0;
  writer.WriteSource(source);
}

void SVPInput::CreateFieldTable(std::string& name, double time,
                                base::SubtableWriter& writer) {
  base::SubtableWriter::FieldInfo field;
  field.name = name;
  field.code = std::string();
  field.time = time;
  field.num_poly = 0;
  field.delay_direction_ra = metadata_.ra_;  // (in radians)
  field.delay_direction_dec = metadata_.dec_;
  field.phase_direction_ra = field.delay_direction_ra;
  field.phase_direction_dec = field.delay_direction_dec;
  field.reference_direction_ra = field.delay_direction_ra;
  field.reference_direction_dec = field.delay_direction_dec;
  field.source_id = -1;
  field.flag_row = false;
  writer.WriteField(field);
}

void SVPInput::CreateObservationTable(double startTime, double endTime,
                                      base::SubtableWriter& writer) {
  base::SubtableWriter::ObservationInfo observation;

  observation.telescope_name = metadata_.telescope_name_;
  observation.start_time = startTime;
  observation.end_time = endTime;
  observation.observer = "Unknown";
  observation.schedule_type = metadata_.telescope_name_;
  ;
  observation.project = "Unknown";
  observation.release_date = 0;
  observation.flag_row = false;

  observation.antenna_type = metadata_.telescope_name_;
  observation.rcu_mode = 0;
  observation.flag_window_size = 0;
  writer.WriteObservation(observation);
}

void SVPInput::CreatePolarizationTable(base::SubtableWriter& writer,
                                       const int nr_polarizations) {
  writer.WriteLinearPolarizations(false, nr_polarizations);
}

void SVPInput::ConnectServer(const char* path) {
  if ((sock_fd_ = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
    throw std::runtime_error("DP3 SVPInput Cannot create unix socket.");
  }

  struct sockaddr_un saddr;
  saddr.sun_family = AF_UNIX;
  strcpy(saddr.sun_path, path);
  int ssize = offsetof(struct sockaddr_un, sun_path) + strlen(saddr.sun_path);

  int conn = connect(sock_fd_, (struct sockaddr*)&saddr, ssize);
  if (conn == -1) {
    close(sock_fd_);
    throw std::runtime_error("DP3 SVPInput Cannot connect unix socket.");
  }
}

template <typename T>
int SVPInput::ReceiveItem(T& message) {
  assert(sock_fd_ != -1);
  T ret;
  char* data = (char*)&ret;

  int left = sizeof(ret);
  int rc;
  /* loop while all data is not received */
  do {
    rc = read(sock_fd_, data, left);
    if (rc <= 0) {
      return -1;
    } else {
      data += rc;
      left -= rc;
    }
  } while (left > 0);

  message = ret;
  return rc;
}

int SVPInput::ReceiveBytes(char* message, const int message_size) {
  assert(sock_fd_ != -1);
  char* data = (char*)message;

  int left = sizeof(char) * (size_t)message_size;
  int rc;
  /* loop while all data is not received */
  do {
    rc = read(sock_fd_, data, left);
    if (rc <= 0) {
      return -1;
    } else {
      data += rc;
      left -= rc;
    }
  } while (left > 0);

  return rc;
}

}  // namespace steps
}  // namespace dp3
