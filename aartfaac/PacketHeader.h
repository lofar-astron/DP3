// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_AARTFAAC_PACKETHEADER_H_
#define DP3_AARTFAAC_PACKETHEADER_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>

namespace dp3::aartfaac {
namespace details {

template <typename T>
void PrintArray(std::ostringstream& str, T* data, size_t n) {
  T previous = data[0];
  int count = 1;

  for (size_t i = 1; i < n; ++i) {
    if (data[i] == previous) {
      ++count;
    } else {
      str << previous;
      if (count > 1) {
        str << "x" << count;
      }
      str << ", ";

      previous = data[i];
      count = 1;
    }
  }

  str << previous;
  if (count > 1) {
    str << "x" << count;
  }
}
}  // namespace details

struct PacketHeader {
  uint32_t magic;
  uint16_t nr_receivers;
  uint8_t nr_polarizations;
  uint8_t correlation_mode;
  double start_time;
  double end_time;
  uint32_t weights[300];  // Fixed-sized field, independent of #stations!
  uint32_t nr_samples_per_integration;
  uint16_t nr_channels;
  char pad0[2];
  double first_channel_frequency;
  double channel_bandwidth;
  char pad1[288];

  PacketHeader()
      : magic(kHeaderMagicTag),
        nr_receivers(0),
        nr_polarizations(0),
        correlation_mode(0),
        start_time(0),
        end_time(0),
        nr_samples_per_integration(0),
        nr_channels(0),
        first_channel_frequency(0),
        channel_bandwidth(0) {
    std::fill_n(weights, std::size(weights), 0);
    std::fill_n(pad0, std::size(pad0), 0);
    std::fill_n(pad1, std::size(pad1), 0);
  }

  uint32_t Magic() const { return magic; }

  double NrReceivers() const { return nr_receivers; }

  uint8_t NrPolarizations() const { return nr_polarizations; }

  uint8_t CorrelationMode() const { return correlation_mode; }

  double StartTime() const { return start_time; }

  double EndTime() const { return end_time; }

  size_t NrWeights() const { return 300; }

  uint32_t* Weights() { return weights; }

  uint16_t NrSamplesPerIntegration() const {
    return nr_samples_per_integration;
  }

  uint16_t NrChannels() const { return nr_channels; }

  double FirstChannelFrequency() const { return first_channel_frequency; }

  double ChannelBandwidth() const { return channel_bandwidth; }

  template <typename T>
  void CopyWeights(T* destination) {
    for (unsigned i = 0; i < NrWeights(); i++) {
      destination[i] = weights[i];
    }
  }

  size_t NrBaselines() const { return NrReceivers() * (NrReceivers() + 1) / 2; }

  size_t VisPerTimestep() const {
    return NrBaselines() * NrChannels() * NrPolarizations();
  }

  void Check() {
    if (magic != Magic()) {
      ThrowMagicError();
    }
    if (correlation_mode != 15) {
      ThrowCorrelationModeError(correlation_mode);
    }
  }

  void ThrowMagicError() {
    throw std::runtime_error(
        "This file does not start with the standard header prefix. It is not "
        "a supported Aartfaac correlation file or is damaged.");
  }

  void ThrowCorrelationModeError(uint8_t correlationMode) {
    std::ostringstream str;
    str << "This Aartfaac file specifes a correlation mode of '"
        << int(correlationMode)
        << "'. This tool can only handle sets with 4 polarizations (mode "
           "15).";
    throw std::runtime_error(str.str());
  }

  std::string ToString() {
    std::ostringstream str;
    str << "magic = " << Magic() << '\n';
    str << "nr receivers = " << NrReceivers() << '\n'
        << "nr polarizations = " << static_cast<unsigned>(NrPolarizations())
        << '\n'
        << "correlation mode = " << static_cast<unsigned>(CorrelationMode())
        << '\n'
        << std::fixed << "start time = " << StartTime() << '\n'
        << "end time = " << EndTime()
        << " (total: " << (EndTime() - StartTime()) << " s)\n"
        << "weights = [";  //
    details::PrintArray(str, Weights(), NrWeights());
    str << "] (" << NrWeights() << "x)\n";
    str << "nr samples per integration = " << NrSamplesPerIntegration() << '\n'
        << "nr channels = " << NrChannels() << '\n';
    str << std::fixed << std::setprecision(3)
        << "first channel frequency = " << first_channel_frequency / 1e6
        << " MHz\n"
        << "channel bandwidth = " << channel_bandwidth / 1e6 << " MHz\n";
    return str.str();
  }

 private:
  /// Magic number in raw corr visibilities.
  static inline constexpr uint32_t kHeaderMagicTag = 0x3B98F003;
};

static_assert(sizeof(PacketHeader) == 1536,
              "Header should be of size 1536 bytes");
}  // namespace dp3::aartfaac
#endif  // AARTFAAC_PACKETHEADER_H_
