// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_AARTFAAC_FILE_H_
#define DP3_AARTFAAC_FILE_H_

#include <complex>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "base/RcuMode.h"

namespace dp3::aartfaac {

struct Timestep {
  double startTime, endTime;
};

/**
 * CASA times are not in MJD, but in s.
 */
constexpr double TimeToCasa(double timestamp) {
  return timestamp + ((2440587.5 - 2400000.5) * 86400.0);
}

// An AARTFAAC file has a header followed by the data, which is written as:
// std::complex<float>
// visibilities[nr_baselines][nr_channels][nr_pols][nr_pols];
class File {
 public:
  File(const char *filename, base::RcuMode mode)
      : filename_(filename), mode_(mode) {
    file_.open(filename);

    file_.seekg(0, std::ios::end);
    filesize_ = file_.tellg();

    // Read first header
    file_.seekg(0, std::ios::beg);
    file_.read(reinterpret_cast<char *>(&header_), sizeof(PacketHeader));

    header_.Check();

    blockSize_ = sizeof(std::complex<float>) * header_.VisPerTimestep();
    bandwidth_ = mode.Bandwidth();
    frequency_ = header_.FirstChannelFrequency();

    SeekToTimestep(0);
  }

  explicit File(const char *filename)
      : filename_(filename), mode_(base::RcuMode::Unused) {
    file_.seekg(0, std::ios::end);
    filesize_ = file_.tellg();

    // Read first header
    file_.seekg(0, std::ios::beg);
    file_.read(reinterpret_cast<char *>(&header_), sizeof(PacketHeader));

    header_.Check();

    blockSize_ = sizeof(std::complex<float>) * header_.VisPerTimestep();

    SeekToTimestep(0);
  }

  void SkipTimesteps(int count) {
    file_.seekg(count * (sizeof(PacketHeader) + blockSize_), std::ios::cur);
    blockPosition_ += count;
  }

  void SeekToTimestep(size_t timestep) {
    file_.seekg(timestep * (sizeof(PacketHeader) + blockSize_), std::ios::beg);
    blockPosition_ = timestep;
  }

  Timestep ReadTimestep(std::complex<float> *buffer) {
    file_.read(reinterpret_cast<char *>(&header_), sizeof(PacketHeader));
    file_.read(reinterpret_cast<char *>(buffer), blockSize_);
    if (!file_) throw std::runtime_error("Error reading file");
    ++blockPosition_;

    return Timestep{TimeToCasa(header_.StartTime()),
                    TimeToCasa(header_.EndTime())};
  }

  Timestep ReadMetadata() {
    PacketHeader h;
    file_.read(reinterpret_cast<char *>(&h), sizeof(PacketHeader));
    if (!file_) throw std::runtime_error("Error reading file");
    SeekToTimestep(blockPosition_);

    return Timestep{TimeToCasa(h.StartTime()), TimeToCasa(h.EndTime())};
  }

  bool HasMore() const { return blockPosition_ < (filesize_ / blockSize_); }

  size_t NTimesteps() const { return filesize_ / blockSize_; }

  size_t VisPerTimestep() const { return header_.VisPerTimestep(); }

  size_t NChannels() const { return header_.NrChannels(); }

  size_t NAntennas() const { return header_.NrReceivers(); }

  size_t NPolarizations() const { return header_.NrPolarizations(); }

  uint8_t CorrelationMode() const { return header_.CorrelationMode(); }

  double Bandwidth() const { return bandwidth_; }
  double StartTime() const { return TimeToCasa(header_.StartTime()); }
  double Frequency() const { return frequency_; }
  double IntegrationTime() const {
    return header_.EndTime() - header_.StartTime();
  }

  const PacketHeader &Header() const { return header_; }

  size_t BlockPos() const { return blockPosition_; }

 private:
  std::ifstream file_;
  size_t blockSize_;
  size_t filesize_;
  size_t blockPosition_ = 0;
  std::string filename_;
  PacketHeader header_;
  base::RcuMode mode_;
  double frequency_;
  double bandwidth_;
};
}  // namespace dp3::aartfaac
#endif  // AARTFAAC_FILE_H_
