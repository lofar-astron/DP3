// AOFlaggerStep.cc: DPPP step class to flag data based on rficonsole
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Andre Offringa, Ger van Diepen

#include "AOFlaggerStep.h"

#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"

#include "../common/Memory.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/OS/File.h>

#include <aocommon/parallelfor.h>
#include <aoflagger.h>

#include <iostream>
#include <algorithm>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

AOFlaggerStep::AOFlaggerStep(InputStep* input,
                             const common::ParameterSet& parset,
                             const string& prefix)
    : name_(prefix),
      buffer_index_(0),
      n_times_(0),
      memory_needed_(0),
      flag_counter_(input->msName(), parset, prefix + "count."),
      move_time_(0),
      flag_time_(0),
      stats_time_(0),
      qstats_() {
  strategy_name_ = parset.getString(prefix + "strategy", string());
  if (strategy_name_.empty())
    strategy_name_ =
        aoflagger_.FindStrategyFile(aoflagger::TelescopeId::LOFAR_TELESCOPE);

  window_size_ = parset.getUint(prefix + "timewindow", 0);
  memory_ = parset.getUint(prefix + "memorymax", 0);
  memory_percentage_ = parset.getUint(prefix + "memoryperc", 0);
  overlap_ = parset.getUint(prefix + "overlapmax", 0);
  // Also look for keyword overlap for backward compatibility.
  if (overlap_ == 0) {
    overlap_ = parset.getUint(prefix + "overlap", 0);
  }
  overlap_percentage_ = parset.getDouble(prefix + "overlapperc", -1);
  flag_auto_correlations_ = parset.getBool(prefix + "autocorr", true);
  collect_statistics_ = parset.getBool(prefix + "keepstatistics", true);
}

AOFlaggerStep::~AOFlaggerStep() {}

void AOFlaggerStep::show(std::ostream& os) const {
  os << "AOFlaggerStep " << name_ << '\n';
  os << "  strategy:       " << strategy_name_ << '\n';
  os << "  timewindow:     " << window_size_ << '\n';
  os << "  overlap:        " << overlap_ << '\n';
  os << "  keepstatistics: " << collect_statistics_ << '\n';
  os << "  autocorr:       " << flag_auto_correlations_ << '\n';
  os << "  max memory used ";
  formatBytes(os, memory_needed_);
  os << '\n';
}

void AOFlaggerStep::formatBytes(std::ostream& os, double bytes) {
  int exp = 0;
  while (bytes >= 1024 && exp < 5) {
    bytes /= 1024;
    exp++;
  }

  unsigned int origPrec = os.precision();
  os.precision(1);

  if (exp == 0) {
    os << std::fixed << bytes << " "
       << "B";
  } else {
    os << std::fixed << bytes << " "
       << "KMGTPE"[exp - 1] << "B";
  }

  os.precision(origPrec);
}

void AOFlaggerStep::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteFlags();
  // Determine available memory.
  double availMemory = casacore::HostInfo::memoryTotal() * 1024.;
  // Determine how much memory can be used.
  double memory = common::AvailableMemory(memory_, memory_percentage_, false);

  // Determine how much buffer space is needed per time slot.
  // The flagger needs 3 extra work buffers (data+flags) per thread.
  double timeSize = (sizeof(casacore::Complex) + sizeof(bool)) *
                    (infoIn.nbaselines() + 3 * getInfo().nThreads()) *
                    infoIn.nchan() * infoIn.ncorr();
  // If no overlap percentage is given, set it to 1%.
  if (overlap_percentage_ < 0 && overlap_ == 0) {
    overlap_percentage_ = 1;
  }
  // If no time window given, determine it from the available memory.
  if (window_size_ == 0) {
    double nt = memory / timeSize;
    if (overlap_percentage_ > 0) {
      // Determine the overlap (add 0.5 for rounding).
      // If itsOverLap is also given, it is the maximum.
      double tw = nt / (1 + 2 * overlap_percentage_ / 100);
      unsigned int overlap =
          (unsigned int)(overlap_percentage_ * tw / 100 + 0.5);
      if (overlap_ == 0 || overlap < overlap_) {
        overlap_ = overlap;
      }
    }
    window_size_ = (unsigned int)(std::max(1., nt - 2 * overlap_));
    // Make the window size divide the nr of times nicely (if known).
    // In that way we cannot have a very small last window.
    if (infoIn.ntime() > 0) {
      unsigned int nwindow = 1 + (infoIn.ntime() - 1) / window_size_;
      window_size_ = 1 + (infoIn.ntime() - 1) / nwindow;
      if (overlap_percentage_ > 0) {
        unsigned int overlap =
            (unsigned int)(overlap_percentage_ * window_size_ / 100 + 0.5);
        if (overlap < overlap_) {
          overlap_ = overlap;
        }
      }
    }
  }
  if (overlap_ == 0) {
    overlap_ = (unsigned int)(overlap_percentage_ * window_size_ / 100);
  }
  // Check if it all fits in memory.
  memory_needed_ = (window_size_ + 2 * overlap_) * timeSize;
  if (memory_needed_ >= availMemory)
    throw std::runtime_error(
        "Timewindow " + std::to_string(window_size_) + " and/or overlap " +
        std::to_string(overlap_) + ' ' + std::to_string(memory) +
        " too large for available memory " + std::to_string(availMemory));
  // Size the buffer (need overlap on both sides).
  buffer_.resize(window_size_ + 2 * overlap_);
  // Initialize the flag counters.
  flag_counter_.init(getInfo());

  // Initialize metadata for aoflagger
  frequencies_ = infoIn.chanFreqs();
  aoflagger::Band band;
  band.id = 0;
  band.channels.resize(frequencies_.size());
  for (size_t ch = 0; ch != frequencies_.size(); ++ch) {
    band.channels[ch].frequency = frequencies_[ch];
    band.channels[ch].width = infoIn.chanWidths()[ch];
  }
  aoflagger_.SetBandList({band});

  const casacore::Vector<casacore::String>& ant_names = infoIn.antennaNames();
  const casacore::Vector<casacore::MPosition>& ant_pos = infoIn.antennaPos();
  std::vector<aoflagger::Antenna> antennas(ant_names.size());
  for (size_t i = 0; i != ant_names.size(); ++i) {
    antennas[i].id = i;
    antennas[i].name = ant_names[i];
    casacore::Vector<double> vec = ant_pos[i].getValue().get();
    antennas[i].x = vec[0];
    antennas[i].y = vec[1];
    antennas[i].z = vec[2];
  }
  aoflagger_.SetAntennaList(std::move(antennas));
}

void AOFlaggerStep::showCounts(std::ostream& os) const {
  os << "\nFlags set by AOFlaggerStep " << name_;
  os << "\n===========================\n";
  flag_counter_.showBaseline(os, n_times_);
  flag_counter_.showChannel(os, n_times_);
  flag_counter_.showCorrelation(os, n_times_);
}

void AOFlaggerStep::showTimings(std::ostream& os, double duration) const {
  double flagDur = timer_.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, flagDur, duration);
  os << " AOFlaggerStep " << name_ << '\n';
  os << "          ";
  // move time and flag time are sum of all threads.
  // Scale them to a single elapsed time.
  double factor =
      (compute_timer_.getElapsed() / (move_time_ + flag_time_ + stats_time_));
  FlagCounter::showPerc1(os, move_time_ * factor, flagDur);
  os << " of it spent in shuffling data" << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, flag_time_ * factor, flagDur);
  os << " of it spent in calculating flags" << '\n';
  if (collect_statistics_) {
    os << "          ";
    FlagCounter::showPerc1(
        os, stats_time_ * factor + quality_timer_.getElapsed(), flagDur);
    os << " of it spent in making quality statistics" << '\n';
  }
}

// Alternative strategy is to flag in windows
//  0 ..  n+2m
//  n .. 2n+2m
// 2n .. 3n+2m  etc.
// and also update the flags in the overlaps
bool AOFlaggerStep::process(const DPBuffer& buf) {
  timer_.start();
  // Accumulate in the time window until the window and overlap are full.
  n_times_++;
  buffer_[buffer_index_].copy(buf);
  ++buffer_index_;
  if (buffer_index_ == window_size_ + 2 * overlap_) {
    flag(2 * overlap_);
  }
  timer_.stop();
  return true;
}

void AOFlaggerStep::finish() {
  std::cerr << "  " << buffer_index_
            << " time slots to finish in AOFlaggerStep ...\n";
  timer_.start();
  // Set window size to all entries left.
  window_size_ = buffer_index_;
  if (window_size_ > 0) {
    // Flag the remaining time slots (without right overlap).
    flag(0);
  }
  buffer_.clear();
  timer_.stop();
  // Let the next step finish its processing.
  getNextStep()->finish();
}

void AOFlaggerStep::addToMS(const string& msName) {
  timer_.start();
  if (collect_statistics_) {
    quality_timer_.start();
    qstats_.WriteStatistics(msName);
    quality_timer_.stop();
  }
  timer_.stop();
  getPrevStep()->addToMS(msName);
}

void AOFlaggerStep::flag(unsigned int rightOverlap) {
  // Get the sizes of the axes.
  size_t nrbl = buffer_[0].getData().shape()[2];
  size_t ncorr = buffer_[0].getData().shape()[0];
  if (ncorr != 4)
    throw std::runtime_error(
        "AOFlaggerStep can only handle all 4 correlations");
  // Get antenna numbers in case applyautocorr is true.
  const casacore::Vector<int>& ant1 = getInfo().getAnt1();
  const casacore::Vector<int>& ant2 = getInfo().getAnt2();
  compute_timer_.start();

  aoflagger::Interval interval;
  interval.id = 0;
  interval.times.resize(buffer_.size());
  for (std::size_t i = 0; i != buffer_.size(); ++i) {
    interval.times[i] = buffer_[i].getTime();
  }
  aoflagger_.SetIntervalList({interval});

  // Now flag each baseline and collect statistics for this time window.
  // The baselines can be processed in parallel.

  struct ThreadData {
    FlagCounter counter;
    aoflagger::QualityStatistics qstats;
    aoflagger::Strategy strategy;
  };
  std::vector<ThreadData> threadData(getInfo().nThreads());
  // Create thread-private objects.
  for (size_t t = 0; t != getInfo().nThreads(); ++t) {
    threadData[t].counter.init(getInfo());
    // Create a statistics object for all polarizations.
    threadData[t].qstats = aoflagger_.MakeQualityStatistics(
        interval.times.data(), interval.times.size(), frequencies_.data(),
        frequencies_.size(), 4, false);
    threadData[t].strategy = aoflagger_.LoadStrategyFile(strategy_name_);
  }

  aocommon::ParallelFor<size_t> loop(getInfo().nThreads());
  loop.Run(0, nrbl, [&](size_t ib, size_t thread) {
    // Do autocorrelations only if told so.
    if (ant1[ib] == ant2[ib]) {
      if (flag_auto_correlations_) {
        flagBaseline(0, window_size_ + rightOverlap, 0, ib,
                     threadData[thread].counter, threadData[thread].strategy,
                     threadData[thread].qstats);
      }
    } else {
      flagBaseline(0, window_size_ + rightOverlap, 0, ib,
                   threadData[thread].counter, threadData[thread].strategy,
                   threadData[thread].qstats);
    }
  });  // end of parallel for

  // Add the counters to the overall object.
  for (size_t t = 0; t != getInfo().nThreads(); ++t) {
    flag_counter_.add(threadData[t].counter);
    if (collect_statistics_) {
      quality_timer_.stop();
      // Add the rfi statistics to the global object.
      if (t == 0) {
        qstats_ = std::move(threadData[t].qstats);
      } else {
        qstats_ += threadData[t].qstats;
      }
      quality_timer_.start();
    }
  }

  compute_timer_.stop();
  timer_.stop();
  // Let the next step process the buffers.
  // If possible, discard the buffer processed to minimize memory usage.
  for (unsigned int i = 0; i < window_size_; ++i) {
    getNextStep()->process(buffer_[i]);
    ///        itsBuf[i] = DPBuffer();
    /// cout << "cleared buffer " << i << '\n';
  }
  timer_.start();
  // Shift the buffers still needed to the beginning of the vector.
  // This is a bit easier than keeping a wrapped vector.
  // Note it is a cheap operation, because shallow copies are made.
  for (unsigned int i = 0; i < rightOverlap; ++i) {
    buffer_[i].copy(buffer_[i + window_size_]);
    /// cout << "moved buffer " <<i+itsWindowSize<<" to "<< i << '\n';
  }
  buffer_index_ = rightOverlap;
}

void AOFlaggerStep::flagBaseline(unsigned int leftOverlap,
                                 unsigned int windowSize,
                                 unsigned int rightOverlap, unsigned int bl,
                                 FlagCounter& counter,
                                 aoflagger::Strategy& strategy,
                                 aoflagger::QualityStatistics& rfiStats) {
  common::NSTimer moveTimer, flagTimer, qualTimer;
  moveTimer.start();
  // Get the sizes of the axes.
  unsigned int ntime = leftOverlap + windowSize + rightOverlap;
  unsigned int nchan = buffer_[0].getData().shape()[1];
  unsigned int blsize = nchan * buffer_[0].getData().shape()[0];
  // Fill the rficonsole buffers and flag.
  // Create the objects for the real and imaginary data of all corr.
  aoflagger::ImageSet imageSet = aoflagger_.MakeImageSet(ntime, nchan, 8);
  imageSet.SetAntennas(info().getAnt1()[bl], info().getAnt2()[bl]);
  imageSet.SetInterval(0);
  imageSet.SetBand(0);
  aoflagger::FlagMask origFlags = aoflagger_.MakeFlagMask(ntime, nchan);
  const unsigned int iStride = imageSet.HorizontalStride();
  const unsigned int fStride = origFlags.HorizontalStride();
  for (unsigned int i = 0; i < ntime; ++i) {
    const casacore::Complex* data = buffer_[i].getData().data() + bl * blsize;
    const bool* flags = buffer_[i].getFlags().data() + bl * blsize;
    for (unsigned int j = 0; j < nchan; ++j) {
      for (unsigned int p = 0; p != 4; ++p) {
        imageSet.ImageBuffer(p * 2)[i + j * iStride] = data->real();
        imageSet.ImageBuffer(p * 2 + 1)[i + j * iStride] = data->imag();
        data++;
      }
      origFlags.Buffer()[i + j * fStride] = *flags;
      flags += 4;
    }
  }
  // Execute the strategy to do the flagging.
  moveTimer.stop();
  flagTimer.start();
  aoflagger::FlagMask rfiMask = strategy.Run(imageSet);
  flagTimer.stop();
  // Put back the true flags and count newly set flags.
  moveTimer.start();
  for (unsigned int i = leftOverlap; i < windowSize + leftOverlap; ++i) {
    bool* flags = buffer_[i].getFlags().data() + bl * blsize;
    for (unsigned int j = 0; j < nchan; ++j) {
      // Only set if not already set.
      // If any corr is newly set, set all corr.
      if (!flags[0]) {
        bool setFlag = true;
        if (rfiMask.Buffer()[i + j * fStride]) {
          counter.incrCorrelation(0);
          counter.incrCorrelation(1);
          counter.incrCorrelation(2);
          counter.incrCorrelation(3);
        } else {
          setFlag = false;
        }
        if (setFlag) {
          counter.incrBaseline(bl);
          counter.incrChannel(j);
          for (int k = 0; k < 4; ++k) {
            flags[k] = true;
          }
        }
      }
      flags += 4;
    }
  }
  moveTimer.stop();
  // Update the RFI statistics if needed.
  if (collect_statistics_) {
    qualTimer.start();
    addStats(rfiStats, imageSet, rfiMask, origFlags, bl);
    qualTimer.stop();
  }
  std::lock_guard<std::mutex> lock(mutex_);
  // Add the timings.
  move_time_ += moveTimer.getElapsed();
  flag_time_ += flagTimer.getElapsed();
  stats_time_ += qualTimer.getElapsed();
}

void AOFlaggerStep::addStats(aoflagger::QualityStatistics& qStats,
                             const aoflagger::ImageSet& values,
                             const aoflagger::FlagMask& rfiMask,
                             const aoflagger::FlagMask& origMask, int bl) {
  qStats.CollectStatistics(values, rfiMask, origMask, getInfo().getAnt1()[bl],
                           getInfo().getAnt2()[bl]);
}

}  // namespace steps
}  // namespace dp3
