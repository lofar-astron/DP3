// BdaExpander.cc: DP3 step class to expand a BdaBuffer to a DPBuffer (BDA data
// to regular data)
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Chiara Salvoni

#include "BdaExpander.h"

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

#include <casacore/casa/BasicMath/Math.h>

#include <dp3/common/Types.h>

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

BdaExpander::BdaExpander(const std::string& prefix)
    : timer_("BDA Expander"), step_name_(prefix) {}

void BdaExpander::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);

  if (!info_in.isBDAIntervalFactorInteger()) {
    throw std::invalid_argument(
        "Invalid info in input data INTEGER_INTERVAL_FACTORS must be true");
  }

  // Calculate single channel non-averaged
  assert(!info_in.BdaChanWidths().empty());
  const std::vector<double>& chan_widths = info_in.chanWidths(0);
  const double total_bw =
      std::accumulate(chan_widths.begin(), chan_widths.end(), 0.0);
  double single_channel_bw = total_bw / info_in.nchan();

  // Define center frequencies and widths
  std::vector<double> freqs;
  freqs.reserve(info_in.nchan());
  std::vector<double> widths(info_in.nchan(), single_channel_bw);
  freqs.push_back(info_in.BdaChanFreqs().front().front() -
                  0.5 * info_in.BdaChanWidths().front().front() +
                  0.5 * single_channel_bw);

  for (unsigned int k = 1; k < info_in.nchan(); k++) {
    freqs.push_back(freqs[k - 1] + single_channel_bw);
  }

  next_time_slot_to_process_ = 0;
  channels_mapping_.resize(info_in.nbaselines());
  for (unsigned int k = 0; k < info_in.nbaselines(); k++) {
    int next_chan = 0;
    // Store the mapping between averaged and non averaged channels in the
    // variable channels_mapping_
    for (size_t i = 0; i < info_in.chanWidths(k).size(); i++) {
      const int n_averaged_channels =
          std::round(info_in.chanWidths(k)[i] / single_channel_bw);
      if (n_averaged_channels == 1) {
        // non-bda-averaged channel
        channels_mapping_[k].push_back(next_chan);

      } else {
        // bda-averaged channels
        for (int p = 0; p < n_averaged_channels; p++) {
          // The same bda-averaged channel will be copied over multiple channels
          // in the regular output
          channels_mapping_[k].push_back(next_chan);
        }
      }
      next_chan++;
    }
  }

  infoOut().setChannels(std::move(freqs), std::move(widths));
  infoOut().setMetaChanged();
}

void BdaExpander::show(std::ostream& os) const {
  os << "BdaExpander " << step_name_ << '\n';
}

void BdaExpander::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " BdaExpander " << step_name_ << '\n';
}

bool BdaExpander::process(std::unique_ptr<base::BdaBuffer> bda_buffer) {
  timer_.start();

  std::vector<BdaBuffer::Row> rows = bda_buffer->GetRows();

  for (const BdaBuffer::Row& row : bda_buffer->GetRows()) {
    const double current_timeslot_start = row.time - row.interval / 2;
    unsigned int current_timeslot_index = static_cast<unsigned int>(
        round((current_timeslot_start - getInfoOut().startTime()) /
              getInfoOut().timeInterval()));

    // RB_elements has one map per each time interval (as we expect it DPBuffer)
    // Check if current interval is already in the map RB_elements
    // Check should not happen for averaged data -> these are treated
    // differently
    if (casacore::near(row.interval, getInfoOut().timeInterval())) {
      auto it = RB_elements.find(current_timeslot_index);
      if (it == RB_elements.end()) {
        // create new element if a RegularBufferElement for the
        // current_time_centroid does not exist yet
        std::tie(it, std::ignore) = RB_elements.emplace(
            current_timeslot_index,
            RegularBufferElement(getInfoOut().nbaselines(),
                                 getInfoOut().ncorr(), getInfoOut().nchan(),
                                 row.time, row.exposure));
      }

      // if time interval is equal to the smallest time interval, the data is
      // not averaged -> add data "as is" to the right DPBuffer
      RegularBufferElement& rb_element = it->second;
      rb_element.baseline[row.baseline_nr] = true;
      CopyData(*bda_buffer, row, rb_element.regular_buffer);
    } else {
      // If time interval is different than original, the data is averaged ->
      // copy data (deep copy) to multiple DPBuffer if not, change time and
      // interval and spread across the number of slots allowed
      double slots_to_fill = row.interval / getInfoOut().timeInterval();
      for (double i = 0; i < slots_to_fill; i++) {
        const double timeslot_start =
            current_timeslot_start + i * getInfoOut().timeInterval();
        const double timeslot_center =
            timeslot_start + getInfoOut().timeInterval() / 2;
        const unsigned int timeslot_index = static_cast<unsigned int>(
            round((timeslot_start - getInfoOut().startTime()) /
                  getInfoOut().timeInterval()));

        auto it = RB_elements.find(timeslot_index);
        if (it == RB_elements.end()) {
          // create new element if a RegularBufferElement for the
          // current_time_centroid does not exist yet
          std::tie(it, std::ignore) = RB_elements.emplace(
              timeslot_index,
              RegularBufferElement(getInfoOut().nbaselines(),
                                   getInfoOut().ncorr(), getInfoOut().nchan(),
                                   timeslot_center,
                                   getInfoOut().timeInterval()));
        }
        RegularBufferElement& rb_element = it->second;
        rb_element.baseline[row.baseline_nr] = true;
        CopyData(*bda_buffer, row, rb_element.regular_buffer, slots_to_fill);
      }
    }
  }

  // checks if the DPBuffer for the next time slot has data for each baseline
  // If true, sends the DPBuffer to the next processing step
  auto RB = RB_elements.find(next_time_slot_to_process_);
  while (RB != RB_elements.end()) {
    auto it =
        find(RB->second.baseline.begin(), RB->second.baseline.end(), false);
    if (it == RB->second.baseline.end()) {
      getNextStep()->process(std::move(RB->second.regular_buffer));

      RB_elements.erase(next_time_slot_to_process_);
      next_time_slot_to_process_++;
      RB = RB_elements.find(next_time_slot_to_process_);
    } else {
      break;
    }
  }
  timer_.stop();

  return false;
}

void BdaExpander::finish() {
  // Check that internal buffer is empty
  if (!RB_elements.size() == 0) {
    RB_elements.clear();
  }

  // Let the next steps finish.
  getNextStep()->finish();
}

void BdaExpander::CopyData(const BdaBuffer& bda_buffer,
                           const BdaBuffer::Row& bda_row,
                           std::unique_ptr<DPBuffer>& buf_out,
                           float time_averaging_factor) {
  DPBuffer::DataType& data = buf_out->GetData();
  DPBuffer::WeightsType& weights = buf_out->GetWeights();
  DPBuffer::FlagsType& flags = buf_out->GetFlags();
  DPBuffer::UvwType& uvw = buf_out->GetUvw();
  const std::size_t current_bl = bda_row.baseline_nr;

  for (unsigned int chan = 0; chan < getInfoOut().nchan(); ++chan) {
    // Set the pointers to the right value: when channel averaging happens, the
    // pointers will have the same values for multiple loops.
    const std::size_t offset =
        bda_row.offset +
        channels_mapping_[current_bl][chan] * getInfoOut().ncorr();

    const std::complex<float>* bda_data = bda_buffer.GetData();
    const float* bda_weights = bda_buffer.GetWeights();
    const bool* bda_flags = bda_buffer.GetFlags();

    if (bda_data) {
      for (unsigned int corr = 0; corr < getInfoOut().ncorr(); ++corr) {
        data(current_bl, chan, corr) = bda_data[offset + corr];
      }
    }
    if (bda_weights) {
      float channel_averaging_factor =
          std::count(channels_mapping_[current_bl].begin(),
                     channels_mapping_[current_bl].end(),
                     channels_mapping_[current_bl][chan]);

      for (unsigned int corr = 0; corr < getInfoOut().ncorr(); ++corr) {
        weights(current_bl, chan, corr) = bda_weights[offset + corr] /
                                          time_averaging_factor /
                                          channel_averaging_factor;
      }
    }
    if (bda_flags) {
      for (unsigned int corr = 0; corr < getInfoOut().ncorr(); ++corr) {
        flags(current_bl, chan, corr) = bda_flags[offset + corr];
      }
    }

    uvw(current_bl, 0) = bda_row.uvw[0];
    uvw(current_bl, 1) = bda_row.uvw[1];
    uvw(current_bl, 2) = bda_row.uvw[2];
  }
}

BdaExpander::RegularBufferElement::RegularBufferElement(size_t n_baselines,
                                                        unsigned int n_corr,
                                                        unsigned int n_chan,
                                                        double current_time,
                                                        double current_exposure)
    : baseline(std::vector<bool>(n_baselines, false)),
      regular_buffer(
          std::make_unique<DPBuffer>(current_time, current_exposure)) {
  const std::array<std::size_t, 3> shape{n_baselines, n_chan, n_corr};
  regular_buffer->GetData().resize(shape);
  regular_buffer->GetWeights().resize(shape);
  regular_buffer->GetFlags().resize(shape);
  regular_buffer->GetUvw().resize({n_baselines, 3});

  regular_buffer->GetData().fill(0.0);
  regular_buffer->GetWeights().fill(0.0);
  regular_buffer->GetFlags().fill(false);
  regular_buffer->GetUvw().fill(0.0);
}
}  // namespace steps
}  // namespace dp3
