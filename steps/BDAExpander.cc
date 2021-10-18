// BDAExpander.cc: DPPP step class to expand a BDABuffer to a DPBuffer (BDA data
// to regular data) Copyright (C) 2021 ASTRON (Netherlands Institute for Radio
// Astronomy) SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Chiara Salvoni

#include "BDAExpander.h"

#include <iostream>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/Types.h"

#include <boost/make_unique.hpp>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Cube.h>

#include <sstream>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

using dp3::base::BDABuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

BDAExpander::BDAExpander(const string &prefix)
    : timer_("BDA Expander"), step_name_(prefix) {}

BDAExpander::~BDAExpander() {}

void BDAExpander::updateInfo(const DPInfo &_info) {
  info() = _info;
  info().setNeedVisData();
  info().setWriteData();

  if (!info().isBDAIntervalFactorInteger()) {
    throw std::invalid_argument(
        "Invalid info in input data INTEGER_INTERVAL_FACTORS must be true");
  }

  if (info().metaChanged()) {
    throw std::invalid_argument("Update step " + step_name_ +
                                " is not possible because meta data changes");
  }

  // Calculate single channel non-averaged
  assert(!info().BdaChanWidths().empty());
  const std::vector<double> &chan_widths = info().chanWidths(0);
  const double total_bw =
      std::accumulate(chan_widths.begin(), chan_widths.end(), 0.0);
  double single_channel_bw = total_bw / info().nchan();

  // Define center frequencies and widths
  std::vector<double> freqs;
  freqs.reserve(info().nchan());
  std::vector<double> widths(info().nchan(), single_channel_bw);
  freqs.push_back(info().BdaChanFreqs().front().front() -
                  0.5 * info().BdaChanWidths().front().front() +
                  0.5 * single_channel_bw);

  for (unsigned int k = 1; k < info().nchan(); k++) {
    freqs.push_back(freqs[k - 1] + single_channel_bw);
  }

  next_time_slot_to_process_ = 0;
  channels_mapping_.resize(info().nbaselines());
  for (unsigned int k = 0; k < info().nbaselines(); k++) {
    int next_chan = 0;
    // Store the mapping between averaged and non averaged channels in the
    // variable channels_mapping_
    for (size_t i = 0; i < info().chanWidths(k).size(); i++) {
      const int n_averaged_channels =
          std::round(info().chanWidths(k)[i] / single_channel_bw);
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

  info().set(std::move(freqs), std::move(widths));
}

void BDAExpander::show(std::ostream &os) const {
  os << "BDAExpander " << step_name_ << '\n';
}

void BDAExpander::showTimings(std::ostream &os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, timer_.getElapsed(), duration);
  os << " BDAExpander " << step_name_ << '\n';
}

bool BDAExpander::process(std::unique_ptr<base::BDABuffer> bda_buffer) {
  timer_.start();

  std::vector<BDABuffer::Row> rows = bda_buffer->GetRows();

  for (std::size_t row_nr = 0; row_nr < rows.size(); ++row_nr) {
    double current_interval = rows[row_nr].interval;
    double current_timeslot_centroid = rows[row_nr].time;
    double current_timeslot_start =
        current_timeslot_centroid - current_interval / 2;
    unsigned int current_timeslot_index = static_cast<unsigned int>(round(
        (current_timeslot_start - info().startTime()) / info().timeInterval()));
    double current_exposure = rows[row_nr].exposure;
    unsigned int current_bl = rows[row_nr].baseline_nr;

    // RB_elements has one map per each time interval (as we expect it DPBuffer)
    // Check if current interval is already in the map RB_elements
    // Check should not happen for averaged data -> these are treated
    // differently
    if (casacore::near(current_interval, info().timeInterval())) {
      auto it = RB_elements.find(current_timeslot_index);
      if (it == RB_elements.end()) {
        // create new element if a RegularBufferElement for the
        // current_time_centroid does not exist yet
        RegularBufferElement RB = RegularBufferElement(
            info().nbaselines(), info().ncorr(), info().nchan(),
            current_timeslot_centroid, current_exposure);
        RB_elements.insert(std::pair<unsigned int, RegularBufferElement>(
            current_timeslot_index, RB));
      }
    }

    if (casacore::near(current_interval, info().timeInterval())) {
      // if time interval is equal to the smallest time interval, the data is
      // not averaged -> add data "as is" to the right DPBuffer
      RB_elements[current_timeslot_index].baseline_[rows[row_nr].baseline_nr] =
          true;
      CopyData(rows[row_nr], RB_elements[current_timeslot_index].regular_buffer,
               current_bl);
    } else {
      // If time interval is different than original, the data is averaged ->
      // copy data (deep copy) to multiple DPBuffer if not, change time and
      // interval and spread across the number of slots allowed
      double slots_to_fill = current_interval / info().timeInterval();
      for (double i = 0; i < slots_to_fill; i++) {
        double timeslot_center =
            (current_timeslot_centroid - current_interval / 2) +
            i * info().timeInterval() + info().timeInterval() / 2;

        double timeslot_start = timeslot_center - info().timeInterval() / 2;
        unsigned int timeslot_index = static_cast<unsigned int>(round(
            (timeslot_start - info().startTime()) / info().timeInterval()));

        auto it = RB_elements.find(timeslot_index);
        if (it == RB_elements.end()) {
          // create new element if a RegularBufferElement for the
          // current_time_centroid does not exist yet
          RegularBufferElement RB = RegularBufferElement(
              info().nbaselines(), info().ncorr(), info().nchan(),
              timeslot_center, info().timeInterval());
          RB_elements.insert(std::pair<unsigned int, RegularBufferElement>(
              timeslot_index, RB));
        }
        RB_elements.at(timeslot_index).baseline_[rows[row_nr].baseline_nr] =
            true;
        CopyData(rows[row_nr], RB_elements[timeslot_index].regular_buffer,
                 current_bl, slots_to_fill);
      }
    }
  }

  // checks if the DPBuffer for the next time slot has data for each baseline
  // If true, sends the DPBuffer to the next processing step
  auto RB = RB_elements.find(next_time_slot_to_process_);
  while (RB != RB_elements.end()) {
    auto it =
        find(RB->second.baseline_.begin(), RB->second.baseline_.end(), false);
    if (it == RB->second.baseline_.end()) {
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

void BDAExpander::finish() {
  // Check that internal buffer is empty
  if (!RB_elements.size() == 0) {
    RB_elements.clear();
  }

  // Let the next steps finish.
  getNextStep()->finish();
}

void BDAExpander::CopyData(const BDABuffer::Row &bda_row, DPBuffer &buf_out,
                           unsigned int current_bl,
                           float time_averaging_factor) {
  casacore::Cube<casacore::Complex> &data = buf_out.getData();
  casacore::Cube<float> &weights = buf_out.getWeights();
  casacore::Cube<bool> &flags = buf_out.getFlags();
  casacore::Matrix<double> &uvw = buf_out.getUVW();

  for (unsigned int chan = 0; chan < info().nchan(); ++chan) {
    // Set the pointers to the right value: when channel averaging happens, the
    // pointers will have the same values for multiple loops.
    const std::complex<float> *pointer_to_data =
        bda_row.data + channels_mapping_[current_bl][chan] * info().ncorr();
    const float *pointer_to_weights =
        bda_row.weights + channels_mapping_[current_bl][chan] * info().ncorr();
    const bool *pointer_to_flags =
        bda_row.flags + channels_mapping_[current_bl][chan] * info().ncorr();

    if (bda_row.data) {
      for (unsigned int corr = 0; corr < info().ncorr(); ++corr) {
        data(corr, chan, current_bl) = *pointer_to_data;
        ++pointer_to_data;
      }
    }
    if (bda_row.weights) {
      float channel_averaging_factor =
          std::count(channels_mapping_[current_bl].begin(),
                     channels_mapping_[current_bl].end(),
                     channels_mapping_[current_bl][chan]);

      for (unsigned int corr = 0; corr < info().ncorr(); ++corr) {
        weights(corr, chan, current_bl) = *pointer_to_weights /
                                          time_averaging_factor /
                                          channel_averaging_factor;
        ++pointer_to_weights;
      }
    }
    if (bda_row.flags) {
      for (unsigned int corr = 0; corr < info().ncorr(); ++corr) {
        flags(corr, chan, current_bl) = *pointer_to_flags;
        ++pointer_to_flags;
      }
    }

    uvw(0, current_bl) = bda_row.uvw[0];
    uvw(1, current_bl) = bda_row.uvw[1];
    uvw(2, current_bl) = bda_row.uvw[2];
  }
}

BDAExpander::RegularBufferElement::RegularBufferElement(
    int n_baseline, unsigned int n_corr, unsigned int n_chan,
    double current_time, double current_exposure) {
  std::vector<bool> baseline(n_baseline, false);
  this->baseline_ = baseline;

  const std::vector<std::size_t> kChannelCounts(n_chan, 1);
  casacore::Cube<casacore::Complex> data(n_corr, n_chan, n_baseline, 0.0);
  casacore::Cube<bool> flags(data.shape(), false);
  casacore::Cube<float> weights(data.shape(), 0.0);
  casacore::Cube<bool> full_res_flags(kChannelCounts.size(), 1, n_baseline,
                                      false);
  casacore::Matrix<double> uvw(3, n_baseline, 0.0);
  regular_buffer.setData(data);
  regular_buffer.setWeights(weights);
  regular_buffer.setFlags(flags);
  regular_buffer.setFullResFlags(full_res_flags);
  regular_buffer.setUVW(uvw);
  regular_buffer.setTime(current_time);
  regular_buffer.setExposure(current_exposure);
}
}  // namespace steps
}  // namespace dp3
