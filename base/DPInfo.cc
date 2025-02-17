// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include <dp3/base/DPInfo.h>

#include <algorithm>
#include <cmath>
#include <numeric>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/STLIO.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <EveryBeam/correctionmode.h>

#include "../common/Epsilon.h"

using casacore::MDirection;
using casacore::MPosition;
using casacore::MS;

namespace dp3 {
namespace base {

DPInfo::DPInfo(unsigned int n_correlations, unsigned int original_n_channels,
               std::string antenna_set)
    : meta_changed_(false),
      ms_name_(),
      data_column_name_(MS::columnName(MS::DATA)),
      flag_column_name_(MS::columnName(MS::FLAG)),
      weight_column_name_(MS::columnName(MS::WEIGHT_SPECTRUM)),
      antenna_set_(std::move(antenna_set)),
      n_correlations_(n_correlations),
      start_channel_(0),
      original_n_channels_(original_n_channels),
      n_channels_(original_n_channels),
      channel_averaging_factor_(1),
      time_averaging_factors_({1}),
      first_time_(0.0),
      last_time_(0.0),
      time_interval_(1.0),
      n_times_(1),
      beam_correction_mode_(static_cast<int>(everybeam::CorrectionMode::kNone)),
      channel_frequencies_(1),  // Ensure that chanFreqs(), chanWidths(),
      channel_widths_(1),       // resolutions() and effectiveBW()
      resolutions_(1),          // can retrieve a first list.
      effective_bandwidth_(1),
      total_bandwidth_(0.0),
      spectral_window_(0),
      polarizations_() {}

void DPInfo::setTimes(double first_time, double last_time,
                      double time_interval) {
  if (last_time < first_time) {
    throw std::invalid_argument("Last time is before the first time");
  }
  if (time_interval <= 0.0) {
    throw std::invalid_argument("Time interval is zero or negative");
  }

  first_time_ = first_time;
  last_time_ = last_time;
  time_interval_ = time_interval;

  // Calculate the amount of timeslots. The 1.5 comes from:
  // - rounding (0.5)
  // - (last_time_ - first_time_) gives the distance between the centroids of
  //   the first and last time slot. We need to add 0.5 interval at the
  //   beginning and 0.5 interval at the end to obtain the entire time span.
  n_times_ = (last_time - first_time) / time_interval + 1.5;
}

void DPInfo::setMsNames(const std::string& ms_name,
                        const std::string& data_column_name,
                        const std::string& flag_column_name,
                        const std::string& weight_column_name) {
  ms_name_ = ms_name;
  data_column_name_ = data_column_name;
  flag_column_name_ = flag_column_name;
  weight_column_name_ = weight_column_name;
}

void DPInfo::setChannels(std::vector<double>&& chan_freqs,
                         std::vector<double>&& chan_widths,
                         std::vector<double>&& resolutions,
                         std::vector<double>&& effective_bw, double ref_freq,
                         int spectral_window) {
  if (chan_freqs.size() != chan_widths.size()) {
    throw std::invalid_argument(
        "Channel width count does not match frequency count");
  }
  if (resolutions.empty()) {
    resolutions = chan_widths;
  } else if (resolutions.size() != chan_freqs.size()) {
    throw std::invalid_argument(
        "Channel resolution count does not match frequency count");
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  } else if (effective_bw.size() != chan_freqs.size()) {
    throw std::invalid_argument(
        "Channel effective bandwidth count does not match frequency count");
  }

  n_channels_ = chan_freqs.size();
  if (ref_freq == 0) {
    // Takes mean of middle elements if n is even; takes middle if odd.
    ref_freq =
        0.5 * (chan_freqs[(n_channels_ - 1) / 2] + chan_freqs[n_channels_ / 2]);
  }
  reference_frequency_ = ref_freq;

  channel_frequencies_.clear();
  channel_widths_.clear();
  resolutions_.clear();
  effective_bandwidth_.clear();

  channel_frequencies_.push_back(std::move(chan_freqs));
  channel_widths_.push_back(std::move(chan_widths));
  resolutions_.push_back(std::move(resolutions));
  effective_bandwidth_.push_back(std::move(effective_bw));

  total_bandwidth_ = std::accumulate(effective_bandwidth_.front().begin(),
                                     effective_bandwidth_.front().end(), 0.0);
  spectral_window_ = spectral_window;
}

void DPInfo::setChannels(std::vector<std::vector<double>>&& chan_freqs,
                         std::vector<std::vector<double>>&& chan_widths,
                         std::vector<std::vector<double>>&& resolutions,
                         std::vector<std::vector<double>>&& effective_bw,
                         double ref_freq, int spectral_window) {
  if (resolutions.empty()) {
    resolutions = chan_widths;
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  }
  if (chan_freqs.size() != nbaselines() || chan_widths.size() != nbaselines() ||
      resolutions.size() != nbaselines() ||
      effective_bw.size() != nbaselines()) {
    throw std::invalid_argument(
        "Invalid baseline count while setting frequency info");
  }
  for (std::size_t i = 0; i < nbaselines(); ++i) {
    if (chan_freqs[i].size() != chan_widths[i].size()) {
      throw std::invalid_argument(
          "Channel width count does not match frequency count for baseline " +
          std::to_string(i));
    }
    if (chan_freqs[i].size() != resolutions[i].size()) {
      throw std::invalid_argument(
          "Channel resolution count does not match frequency count for "
          "baseline " +
          std::to_string(i));
    }
    if (chan_freqs[i].size() != effective_bw[i].size()) {
      throw std::invalid_argument(
          "Channel effective bandwidth count does not match frequency count "
          "for baseline " +
          std::to_string(i));
    }
  }

  const double total_bw = std::accumulate(effective_bw.front().begin(),
                                          effective_bw.front().end(), 0.0);
  for (std::vector<double>& eff_bw_bl : effective_bw) {
    if (std::accumulate(eff_bw_bl.begin(), eff_bw_bl.end(), 0.0) != total_bw) {
      throw std::runtime_error("Total BW is not equal for all baselines");
    }
  }

  // Find the baseline with the most channels.
  auto comp = [](const std::vector<double>& left,
                 const std::vector<double>& right) {
    return left.size() < right.size();
  };
  auto it = std::max_element(chan_freqs.begin(), chan_freqs.end(), comp);
  n_channels_ = it->size();

  if (ref_freq == 0) {
    // Takes mean of middle elements if n is even; takes middle if odd.
    ref_freq = 0.5 * ((*it)[(n_channels_ - 1) / 2] + (*it)[n_channels_ / 2]);
  }
  reference_frequency_ = ref_freq;

  channel_frequencies_ = std::move(chan_freqs);
  channel_widths_ = std::move(chan_widths);
  resolutions_ = std::move(resolutions);
  effective_bandwidth_ = std::move(effective_bw);
  total_bandwidth_ = total_bw;
  spectral_window_ = spectral_window;
}

bool DPInfo::channelsAreRegular() const {
  if (channel_frequencies_.empty()) {
    return true;
  }

  // Check that all baselines have equal channel layouts.
  const double kTolerance = 1.0;  // Hz
  for (std::size_t bl = 1; bl < channel_frequencies_.size(); ++bl) {
    if (!common::EpsilonEqual(channel_frequencies_.front(),
                              channel_frequencies_[bl], kTolerance) ||
        !common::EpsilonEqual(channel_widths_.front(), channel_widths_[bl],
                              kTolerance) ||
        !common::EpsilonEqual(resolutions_.front(), resolutions_[bl],
                              kTolerance) ||
        !common::EpsilonEqual(effective_bandwidth_.front(),
                              effective_bandwidth_[bl], kTolerance)) {
      return false;
    }
  }

  // Check that channels are evenly spaced.
  const std::vector<double>& freqs = channel_frequencies_.front();
  const std::vector<double>& widths = channel_widths_.front();
  if (freqs.size() > 1) {
    const double freqstep0 = freqs[1] - freqs[0];
    const double kTolerance = 1.e3;  // Compare up to 1kHz accuracy.
    for (std::size_t i = 1; i < freqs.size(); ++i) {
      if ((std::abs(freqs[i] - freqs[i - 1] - freqstep0) >= kTolerance) ||
          (std::abs(widths[i] - widths[0]) >= kTolerance)) {
        return false;
      }
    }
  }

  return true;
}

void DPInfo::setArrayInformation(const MPosition& arrayPos,
                                 const MDirection& phaseCenter,
                                 const MDirection& delayCenter,
                                 const MDirection& tileBeamDir) {
  array_position_ = arrayPos;
  original_phase_center_ = phaseCenter;
  phase_center_ = phaseCenter;
  delay_center_ = delayCenter;
  tile_beam_direction_ = tileBeamDir;
}

void DPInfo::setAntennas(const std::vector<std::string>& antNames,
                         const std::vector<double>& antDiam,
                         const std::vector<MPosition>& antPos,
                         const std::vector<int>& ant1,
                         const std::vector<int>& ant2) {
  if (antNames.size() != antDiam.size() || antNames.size() != antPos.size())
    throw std::invalid_argument(
        "The name, diameter and position arrays are not of the same size");
  if (ant1.size() != ant2.size())
    throw std::invalid_argument(
        "The ant1 and ant2 arrays are not of the same size");
  antenna_names_ = antNames;
  antenna_diameters_ = antDiam;
  antenna_positions_ = antPos;
  antenna1_ = ant1;
  antenna2_ = ant2;
  // Set which antennae are used.
  setAntUsed();
}

void DPInfo::setAntUsed() {
  antennas_used_.clear();
  antenna_map_.resize(antenna_names_.size());
  std::fill(antenna_map_.begin(), antenna_map_.end(), -1);
  for (unsigned int i = 0; i < antenna1_.size(); ++i) {
    if (antenna1_[i] >= static_cast<int>(antenna_map_.size()) ||
        antenna2_[i] >= static_cast<int>(antenna_map_.size()))
      throw std::runtime_error("Antenna map has an inconsistent size");
    antenna_map_[antenna1_[i]] = 0;
    antenna_map_[antenna2_[i]] = 0;
  }
  antennas_used_.reserve(antenna_names_.size());
  for (unsigned int i = 0; i < antenna_map_.size(); ++i) {
    if (antenna_map_[i] == 0) {
      antenna_map_[i] = antennas_used_.size();
      antennas_used_.push_back(i);
    }
  }
}

casacore::MeasureHolder DPInfo::copyMeasure(
    const casacore::MeasureHolder fromMeas) {
  casacore::Record rec;
  casacore::String msg;
  if (!fromMeas.toRecord(msg, rec))
    throw std::runtime_error("Could not copy MeasureHolder record to record");
  casacore::MeasureHolder mh2;
  if (!mh2.fromRecord(msg, rec))
    throw std::runtime_error("Could not copy record to MeasureHolder");
  return mh2;
}

unsigned int DPInfo::update(unsigned int chanAvg, unsigned int timeAvg) {
  if (channel_frequencies_.size() != 1) {
    throw std::runtime_error("Averaging does not support BDA");
  }

  if (chanAvg > n_channels_) {
    chanAvg = n_channels_;
  }
  if (timeAvg > n_times_) {
    timeAvg = n_times_;
  }
  if (n_channels_ % chanAvg != 0)
    throw std::runtime_error(
        "When averaging, nr of channels must divide integrally; "
        "nr of channels = " +
        std::to_string(n_channels_) +
        " averaging factor = " + std::to_string(chanAvg));
  channel_averaging_factor_ *= chanAvg;
  n_channels_ = (n_channels_ + chanAvg - 1) / chanAvg;
  time_averaging_factors_.front() *= timeAvg;
  n_times_ = (n_times_ + timeAvg - 1) / timeAvg;
  // Adjust first_time_ to be the centroid of the first averaged interval.
  // Subtract 0.5 * old interval; Add 0.5 * new interval. Same for last_time_.
  const double time_adjustment = 0.5 * (timeAvg - 1) * time_interval_;
  first_time_ += time_adjustment;
  last_time_ -= time_adjustment;
  time_interval_ *= timeAvg;

  std::vector<double> freqs(n_channels_);
  std::vector<double> widths(n_channels_, 0.0);
  std::vector<double> resols(n_channels_, 0.0);
  std::vector<double> effBWs(n_channels_, 0.0);
  double totBW = 0;
  for (unsigned int i = 0; i < n_channels_; ++i) {
    freqs[i] = 0.5 * (channel_frequencies_.front()[i * chanAvg] +
                      channel_frequencies_.front()[(i + 1) * chanAvg - 1]);
    for (unsigned int j = 0; j < chanAvg; ++j) {
      widths[i] += channel_widths_.front()[i * chanAvg + j];
      resols[i] += resolutions_.front()[i * chanAvg + j];
      effBWs[i] += effective_bandwidth_.front()[i * chanAvg + j];
    }
    totBW += effBWs[i];
  }
  channel_frequencies_.front() = std::move(freqs);
  channel_widths_.front() = std::move(widths);
  resolutions_.front() = std::move(resols);
  effective_bandwidth_.front() = std::move(effBWs);
  total_bandwidth_ = totBW;
  return chanAvg;
}

void DPInfo::update(std::vector<unsigned int>&& timeAvg) {
  time_averaging_factors_ = std::move(timeAvg);
}

void DPInfo::SelectChannels(unsigned int start_channel,
                            unsigned int n_channels) {
  if (start_channel == 0 && n_channels == n_channels_) return;

  if (channel_frequencies_.size() != 1) {
    throw std::runtime_error("Channel selection does not support BDA");
  }

  assert(channel_frequencies_.front().size() ==
             channel_widths_.front().size() &&
         "The number of elements of the channel frequencies and channel widths "
         "should be equal.");
  assert(channel_frequencies_.front().size() == resolutions_.front().size() &&
         "The number of elements of the channel frequencies and resolutions "
         "should be equal.");
  assert(channel_frequencies_.front().size() ==
             effective_bandwidth_.front().size() &&
         "The number of elements of the channel frequencies and effective "
         "bandwidths should be equal.");
  if (start_channel + n_channels > channel_frequencies_.front().size()) {
    throw std::invalid_argument("Channel range is out of bounds.");
  }

  auto freqs_begin = channel_frequencies_.front().begin() + start_channel;
  auto widths_begin = channel_widths_.front().begin() + start_channel;
  auto resol_begin = resolutions_.front().begin() + start_channel;
  auto effbw_begin = effective_bandwidth_.front().begin() + start_channel;
  channel_frequencies_.front() =
      std::vector<double>(freqs_begin, freqs_begin + n_channels);
  channel_widths_.front() =
      std::vector<double>(widths_begin, widths_begin + n_channels);
  resolutions_.front() =
      std::vector<double>(resol_begin, resol_begin + n_channels);
  effective_bandwidth_.front() =
      std::vector<double>(effbw_begin, effbw_begin + n_channels);

  // Add the new start channel to an existing start_channel, so MSUpdater can
  // still update the correct channel(s) in the original input MS.
  start_channel_ += start_channel;
  n_channels_ = n_channels;
}

void DPInfo::SelectBaselines(const std::vector<unsigned int>& baselines) {
  std::vector<int> ant1(baselines.size());
  std::vector<int> ant2(baselines.size());
  for (unsigned int i = 0; i < baselines.size(); ++i) {
    ant1[i] = antenna1_[baselines[i]];
    ant2[i] = antenna2_[baselines[i]];
  }
  antenna1_ = std::move(ant1);
  antenna2_ = std::move(ant2);
  // Clear; they'll be recalculated if needed.
  baseline_lengths_.clear();
  auto_correlation_indices_.clear();

  setAntUsed();
}

void DPInfo::RemoveUnusedAntennas() {
  if (antennas_used_.size() < antenna_map_.size()) {
    // First remove stations.
    std::vector<std::string> names(antennas_used_.size());
    std::vector<double> diameters(antennas_used_.size());
    std::vector<MPosition> positions;
    positions.reserve(antennas_used_.size());
    for (unsigned int i = 0; i < antennas_used_.size(); ++i) {
      names[i] = antenna_names_[antennas_used_[i]];
      diameters[i] = antenna_diameters_[antennas_used_[i]];
      positions.push_back(antenna_positions_[antennas_used_[i]]);
    }
    // Use the new vectors.
    antenna_names_ = std::move(names);
    antenna_diameters_ = std::move(diameters);
    antenna_positions_ = std::move(positions);
    // Renumber the baselines.
    for (unsigned int i = 0; i < antenna1_.size(); ++i) {
      antenna1_[i] = antenna_map_[antenna1_[i]];
      antenna2_[i] = antenna_map_[antenna2_[i]];
    }

    setMetaChanged();

    // Now fill the antennas_used_ and antenna_map_ vectors again.
    setAntUsed();
    // Clear; they'll be recalculated if needed.
    baseline_lengths_.clear();
    auto_correlation_indices_.clear();
  }
}

const std::vector<std::string> DPInfo::GetUsedAntennaNames() const {
  std::vector<std::string> used_antenna_names;
  used_antenna_names.reserve(antennas_used_.size());
  for (size_t used_antenna : antennas_used_) {
    used_antenna_names.emplace_back(antenna_names_[used_antenna]);
  }
  return used_antenna_names;
}

const std::vector<double>& DPInfo::getBaselineLengths() const {
  // Calculate the baseline lengths if not done yet.
  if (baseline_lengths_.empty()) {
    // First get the antenna positions.
    const std::vector<MPosition>& antPos = antennaPos();
    std::vector<casacore::Vector<double>> antVec;
    antVec.reserve(antPos.size());
    for (std::vector<MPosition>::const_iterator iter = antPos.begin();
         iter != antPos.end(); ++iter) {
      // Convert to ITRF and keep as x,y,z in m.
      antVec.push_back(
          MPosition::Convert(*iter, MPosition::ITRF)().getValue().getValue());
    }
    // Fill in the length of each baseline.
    std::vector<double> blength;
    baseline_lengths_.reserve(antenna1_.size());
    for (unsigned int i = 0; i < antenna1_.size(); ++i) {
      casacore::Array<double> diff(antVec[antenna2_[i]] - antVec[antenna1_[i]]);
      baseline_lengths_.push_back(sqrt(sum(diff * diff)));
    }
  }
  return baseline_lengths_;
}

const std::vector<int>& DPInfo::getAutoCorrIndex() const {
  if (auto_correlation_indices_.empty()) {
    int nant =
        1 + std::max(*std::max_element(antenna1_.begin(), antenna1_.end()),
                     *std::max_element(antenna2_.begin(), antenna2_.end()));
    auto_correlation_indices_.resize(nant);
    std::fill(auto_correlation_indices_.begin(),
              auto_correlation_indices_.end(), -1);
    // Keep the baseline table index for the autocorrelations.
    for (unsigned int i = 0; i < antenna1_.size(); ++i) {
      if (antenna1_[i] == antenna2_[i]) {
        auto_correlation_indices_[antenna1_[i]] = i;
      }
    }
  }
  return auto_correlation_indices_;
}

base::Direction DPInfo::phaseCenterDirection() const {
  const MDirection j2000_direction(
      MDirection::Convert(phase_center_, MDirection::J2000)());
  const casacore::Quantum<casacore::Vector<double>> angles =
      j2000_direction.getAngle();
  return {angles.getBaseValue()[0], angles.getBaseValue()[1]};
}

}  // namespace base
}  // namespace dp3
