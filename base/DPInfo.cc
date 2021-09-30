// DPInfo.cc: General info about DPPP data processing attributes like averaging
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DPInfo.h"
#include "Exceptions.h"

#include "../common/Epsilon.h"

#include "../steps/InputStep.h"

#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/STLIO.h>

#include <aocommon/threadpool.h>

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace base {

DPInfo::DPInfo()
    : need_data_(false),
      write_data_(false),
      write_flags_(false),
      write_weights_(false),
      meta_changed_(false),
      n_correlations_(0),
      start_channel_(0),
      n_channels_(0),
      channel_averaging_factor_(1),
      n_times_(0),
      time_averaging_factors_({1}),
      start_time_(0),
      time_interval_(0),
      phase_center_is_original_(true),
      beam_correction_mode_(everybeam::CorrectionMode::kNone),
      n_threads_(aocommon::ThreadPool::NCPUs()) {}

void DPInfo::init(unsigned int ncorr, unsigned int startChan,
                  unsigned int nchan, unsigned int ntime, double startTime,
                  double timeInterval, const string& msName,
                  const string& antennaSet) {
  n_correlations_ = ncorr;
  start_channel_ = startChan;
  n_channels_ = nchan;
  original_n_channels_ = nchan;
  n_times_ = ntime;
  start_time_ = startTime;
  time_interval_ = timeInterval;
  ms_name_ = msName;
  antenna_set_ = antennaSet;
}

void DPInfo::set(std::vector<double>&& chan_freqs,
                 std::vector<double>&& chan_widths,
                 std::vector<double>&& resolutions,
                 std::vector<double>&& effective_bw, double ref_freq) {
  if (resolutions.empty()) {
    resolutions = chan_widths;
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  }

  if (ref_freq == 0) {
    int n = chan_freqs.size();
    // Takes mean of middle elements if n is even; takes middle if odd.
    ref_freq = 0.5 * (chan_freqs[(n - 1) / 2] + chan_freqs[n / 2]);
  }

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
  reference_frequency_ = ref_freq;
}

void DPInfo::set(std::vector<std::vector<double>>&& chan_freqs,
                 std::vector<std::vector<double>>&& chan_widths,
                 std::vector<std::vector<double>>&& resolutions,
                 std::vector<std::vector<double>>&& effective_bw,
                 double ref_freq) {
  if (resolutions.empty()) {
    resolutions = chan_widths;
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  }
  if (chan_freqs.size() != nbaselines() || chan_widths.size() != nbaselines() ||
      resolutions.size() != nbaselines() ||
      effective_bw.size() != nbaselines()) {
    throw Exception("Invalid baseline count while setting frequency info");
  }

  const double total_bw = std::accumulate(effective_bw.front().begin(),
                                          effective_bw.front().end(), 0.0);
  for (std::vector<double>& eff_bw_bl : effective_bw) {
    if (std::accumulate(eff_bw_bl.begin(), eff_bw_bl.end(), 0.0) != total_bw) {
      throw Exception("Total BW is not equal for all baselines");
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

  channel_frequencies_ = std::move(chan_freqs);
  channel_widths_ = std::move(chan_widths);
  resolutions_ = std::move(resolutions);
  effective_bandwidth_ = std::move(effective_bw);
  total_bandwidth_ = total_bw;
  reference_frequency_ = ref_freq;
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

void DPInfo::set(const MPosition& arrayPos, const MDirection& phaseCenter,
                 const MDirection& delayCenter, const MDirection& tileBeamDir) {
  array_position_ = arrayPos;
  phase_center_ = phaseCenter;
  delay_center_ = delayCenter;
  tile_beam_direction_ = tileBeamDir;
}

void DPInfo::set(const Vector<casacore::String>& antNames,
                 const Vector<Double>& antDiam, const vector<MPosition>& antPos,
                 const Vector<Int>& ant1, const Vector<Int>& ant2) {
  if (antNames.size() != antDiam.size() || antNames.size() != antPos.size())
    throw std::invalid_argument(
        "The name, diameter and position arrays are not of the same size");
  if (ant1.size() != ant2.size())
    throw std::invalid_argument(
        "The ant1 and ant2 arrays are not of the same size");
  antenna_names_.reference(antNames);
  antenna_diameters_.reference(antDiam);
  antenna_positions_ = antPos;
  antenna1_ = std::vector<std::size_t>(ant1.begin(), ant1.end());
  antenna2_ = std::vector<std::size_t>(ant2.begin(), ant2.end());
  // Set which antennae are used.
  setAntUsed();
}

void DPInfo::setAntUsed() {
  antennas_used_.clear();
  antenna_map_.resize(antenna_names_.size());
  std::fill(antenna_map_.begin(), antenna_map_.end(), -1);
  for (unsigned int i = 0; i < antenna1_.size(); ++i) {
    if (antenna1_[i] >= antenna_map_.size() ||
        antenna2_[i] >= antenna_map_.size())
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

MeasureHolder DPInfo::copyMeasure(const MeasureHolder fromMeas) {
  Record rec;
  String msg;
  if (!fromMeas.toRecord(msg, rec))
    throw std::runtime_error("Could not copy MeasureHolder record to record");
  MeasureHolder mh2;
  if (!mh2.fromRecord(msg, rec))
    throw std::runtime_error("Could not copy record to MeasureHolder");
  return mh2;
}

unsigned int DPInfo::update(unsigned int chanAvg, unsigned int timeAvg) {
  if (channel_frequencies_.size() != 1) {
    throw Exception("Averaging does not support BDA");
  }

  if (chanAvg > n_channels_) {
    chanAvg = n_channels_;
  }
  if (timeAvg > n_times_) {
    timeAvg = n_times_;
  }
  if (n_channels_ % chanAvg != 0)
    throw Exception(
        "When averaging, nr of channels must divide integrally; "
        "nr of channels = " +
        std::to_string(n_channels_) +
        " averaging factor = " + std::to_string(chanAvg));
  channel_averaging_factor_ *= chanAvg;
  n_channels_ = (n_channels_ + chanAvg - 1) / chanAvg;
  time_averaging_factors_.front() *= timeAvg;
  n_times_ = (n_times_ + timeAvg - 1) / timeAvg;
  time_interval_ *= timeAvg;
  std::vector<double> freqs(n_channels_);
  std::vector<double> widths(n_channels_, 0.);
  std::vector<double> resols(n_channels_, 0.);
  std::vector<double> effBWs(n_channels_, 0.);
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

void DPInfo::update(unsigned int startChan, unsigned int nchan,
                    const vector<unsigned int>& baselines, bool removeAnt) {
  if (channel_frequencies_.size() != 1) {
    throw Exception("Channel selection does not support BDA");
  }
  start_channel_ = startChan;
  auto freqs_begin = channel_frequencies_.front().begin() + startChan;
  auto widths_begin = channel_widths_.front().begin() + startChan;
  auto resol_begin = resolutions_.front().begin() + startChan;
  auto effbw_begin = effective_bandwidth_.front().begin() + startChan;
  channel_frequencies_.front() =
      std::vector<double>(freqs_begin, freqs_begin + nchan);
  channel_widths_.front() =
      std::vector<double>(widths_begin, widths_begin + nchan);
  resolutions_.front() = std::vector<double>(resol_begin, resol_begin + nchan);
  effective_bandwidth_.front() =
      std::vector<double>(effbw_begin, effbw_begin + nchan);
  n_channels_ = nchan;
  // Keep only selected baselines.
  if (!baselines.empty()) {
    std::vector<std::size_t> ant1(baselines.size());
    std::vector<std::size_t> ant2(baselines.size());
    for (unsigned int i = 0; i < baselines.size(); ++i) {
      ant1[i] = antenna1_[baselines[i]];
      ant2[i] = antenna2_[baselines[i]];
    }
    antenna1_ = std::move(ant1);
    antenna2_ = std::move(ant2);
    // Clear; they'll be recalculated if needed.
    baseline_lengths_.clear();
    auto_correlation_indices_.clear();
  }
  setAntUsed();
  // If needed, remove the stations and renumber the baselines.
  if (removeAnt) {
    removeUnusedAnt();
  }
}

void DPInfo::removeUnusedAnt() {
  if (antennas_used_.size() < antenna_map_.size()) {
    // First remove stations.
    Vector<casacore::String> names(antennas_used_.size());
    Vector<Double> diameters(antennas_used_.size());
    vector<MPosition> positions;
    positions.reserve(antennas_used_.size());
    for (unsigned int i = 0; i < antennas_used_.size(); ++i) {
      names[i] = antenna_names_[antennas_used_[i]];
      diameters[i] = antenna_diameters_[antennas_used_[i]];
      positions.push_back(antenna_positions_[antennas_used_[i]]);
    }
    // Use the new vectors.
    antenna_names_.reference(names);
    antenna_diameters_.reference(diameters);
    antenna_positions_.swap(positions);
    // Renumber the baselines.
    for (unsigned int i = 0; i < antenna1_.size(); ++i) {
      antenna1_[i] = antenna_map_[antenna1_[i]];
      antenna2_[i] = antenna_map_[antenna2_[i]];
    }
    // Now fill the antennas_used_ and antenna_map_ vectors again.
    setAntUsed();
    // Clear; they'll be recalculated if needed.
    baseline_lengths_.clear();
    auto_correlation_indices_.clear();
  }
}

const vector<double>& DPInfo::getBaselineLengths() const {
  // Calculate the baseline lengths if not done yet.
  if (baseline_lengths_.empty()) {
    // First get the antenna positions.
    const vector<MPosition>& antPos = antennaPos();
    vector<Vector<double>> antVec;
    antVec.reserve(antPos.size());
    for (vector<MPosition>::const_iterator iter = antPos.begin();
         iter != antPos.end(); ++iter) {
      // Convert to ITRF and keep as x,y,z in m.
      antVec.push_back(
          MPosition::Convert(*iter, MPosition::ITRF)().getValue().getValue());
    }
    // Fill in the length of each baseline.
    vector<double> blength;
    baseline_lengths_.reserve(antenna1_.size());
    for (unsigned int i = 0; i < antenna1_.size(); ++i) {
      Array<double> diff(antVec[antenna2_[i]] - antVec[antenna1_[i]]);
      baseline_lengths_.push_back(sqrt(sum(diff * diff)));
    }
  }
  return baseline_lengths_;
}

const vector<int>& DPInfo::getAutoCorrIndex() const {
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

Record DPInfo::toRecord() const {
  Record rec;
  rec.define("NeedVisData", need_data_);
  rec.define("WriteData", write_data_);
  rec.define("WriteFlags", write_flags_);
  rec.define("WriteWeights", write_weights_);
  rec.define("MetaChanged", meta_changed_);
  rec.define("MSName", ms_name_);
  rec.define("AntennaSet", antenna_set_);
  rec.define("NCorr", n_correlations_);
  rec.define("StartChan", start_channel_);
  rec.define("OrigNChan", original_n_channels_);
  rec.define("NChan", n_channels_);
  rec.define("ChanAvg", channel_averaging_factor_);
  rec.define("NTime", n_times_);
  rec.define("TimeAvg", time_averaging_factors_.front());
  rec.define("StartTime", start_time_);
  rec.define("TimeInterval", time_interval_);
  rec.define("ChanFreqs",
             casacore::Vector<double>(channel_frequencies_.front()));
  rec.define("ChanWidths", casacore::Vector<double>(channel_widths_.front()));
  rec.define("Resolutions", casacore::Vector<double>(resolutions_.front()));
  rec.define("EffectiveBW",
             casacore::Vector<double>(effective_bandwidth_.front()));
  rec.define("TotalBW", total_bandwidth_);
  rec.define("RefFreq", reference_frequency_);
  rec.define("AntNames", antenna_names_);
  rec.define("AntDiam", antenna_diameters_);
  rec.define("AntUsed", Vector<int>(antennas_used_));
  rec.define("AntMap", Vector<int>(antenna_map_));
  rec.define("Ant1", Vector<int>(antenna1_));
  rec.define("Ant2", Vector<int>(antenna2_));
  rec.define("BLength", Vector<double>(baseline_lengths_));
  rec.define("AutoCorrIndex", Vector<int>(auto_correlation_indices_));
  return rec;
}

void DPInfo::fromRecord(const Record& rec) {
  if (rec.isDefined("NeedVisData")) {
    rec.get("NeedVisData", need_data_);
  }
  if (rec.isDefined("WriteData")) {
    rec.get("WriteData", write_data_);
  }
  if (rec.isDefined("WriteFlags")) {
    rec.get("WriteFlags", write_flags_);
  }
  if (rec.isDefined("WriteWeights")) {
    rec.get("WriteWeights", write_weights_);
  }
  if (rec.isDefined("MetaChanged")) {
    rec.get("MetaChanged", meta_changed_);
  }
  if (rec.isDefined("MSName")) {
    ms_name_ = rec.asString("MSName");
  }
  if (rec.isDefined("AntennaSet")) {
    antenna_set_ = rec.asString("AntennaSet");
  }
  if (rec.isDefined("NCorr")) {
    rec.get("NCorr", n_correlations_);
  }
  if (rec.isDefined("StartChan")) {
    rec.get("StartChan", start_channel_);
  }
  if (rec.isDefined("OrigNChan")) {
    rec.get("OrigNChan", original_n_channels_);
  }
  if (rec.isDefined("NChan")) {
    rec.get("NChan", n_channels_);
  }
  if (rec.isDefined("ChanAvg")) {
    rec.get("ChanAvg", channel_averaging_factor_);
  }
  if (rec.isDefined("NTime")) {
    rec.get("NTime", n_times_);
  }
  if (rec.isDefined("TimeAvg")) {
    rec.get("TimeAvg", time_averaging_factors_.front());
  }
  if (rec.isDefined("StartTime")) {
    rec.get("StartTime", start_time_);
  }
  if (rec.isDefined("TimeInterval")) {
    rec.get("TimeInterval", time_interval_);
  }
  if (rec.isDefined("ChanFreqs")) {
    channel_frequencies_.clear();
    channel_frequencies_.push_back(rec.toArrayDouble("ChanFreqs").tovector());
  }
  if (rec.isDefined("ChanWidths")) {
    channel_widths_.clear();
    channel_widths_.push_back(rec.toArrayDouble("ChanWidths").tovector());
  }
  if (rec.isDefined("Resolutions")) {
    resolutions_.clear();
    resolutions_.push_back(rec.toArrayDouble("Resolutions").tovector());
  }
  if (rec.isDefined("EffectiveBW")) {
    effective_bandwidth_.clear();
    effective_bandwidth_.push_back(rec.toArrayDouble("EffectiveBW").tovector());
  }
  if (rec.isDefined("TotalBW")) {
    rec.get("TotalBW", total_bandwidth_);
  }
  if (rec.isDefined("RefFreq")) {
    rec.get("RefFreq", reference_frequency_);
  }
  if (rec.isDefined("AntNames")) {
    rec.get("AntNames", antenna_names_);
  }
  if (rec.isDefined("AntDiam")) {
    rec.get("AntDiam", antenna_diameters_);
  }
  /// if (rec.isDefined ("AntUsed")) {
  /// antennas_used_ = rec.toArrayInt("AntUsed").tovector();
  ///}
  /// if (rec.isDefined ("AntMap")) {
  ///  antenna_map_ = rec.toArrayInt("AntMap").tovector();
  ///}
  if (rec.isDefined("Ant1")) {
    casacore::Vector<casacore::Int> ant1 = rec.toArrayInt("Ant1");
    antenna1_ = std::vector<std::size_t>(ant1.begin(), ant1.end());
  }
  if (rec.isDefined("Ant2")) {
    casacore::Vector<casacore::Int> ant2 = rec.toArrayInt("Ant2");
    antenna2_ = std::vector<std::size_t>(ant2.begin(), ant2.end());
  }
  /// if (rec.isDefined ("BLength")) {
  ///  baseline_lengths_ = rec.toArrayDouble("BLength").tovector();
  ///}
  /// if (rec.isDefined ("AutoCorrIndex")) {
  ///  auto_correlation_indices_ = rec.toArrayInt("AutoCorrIndex").tovector();
  ///}
}

}  // namespace base
}  // namespace dp3
