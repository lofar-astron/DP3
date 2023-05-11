// MadFlagger.cc: DP3 step class to flag data based on median filtering
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "MadFlagger.h"

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../base/DPLogger.h"

#include "../common/Median.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Containers/RecordField.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>

#include <aocommon/parallelfor.h>

#include <xtensor/xcomplex.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>

using casacore::Record;
using casacore::RecordFieldPtr;
using casacore::RecordGram;
using casacore::TableExprNode;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

using dp3::common::operator<<;

namespace dp3 {
namespace steps {

MadFlagger::MadFlagger(const common::ParameterSet& parset,
                       const std::string& prefix)
    : itsName(prefix),
      itsThresholdStr(parset.getString(prefix + "threshold", "1")),
      itsFreqWindowStr(parset.getString(prefix + "freqwindow", "1")),
      itsTimeWindowStr(parset.getString(prefix + "timewindow", "1")),
      itsNTimes(0),
      itsNTimesDone(0),
      itsFlagCounter(parset, prefix + "count."),
      itsMoveTime(0),
      itsMedianTime(0) {
  itsFlagCorr = parset.getUintVector(prefix + "correlations",
                                     std::vector<unsigned int>());
  itsApplyAutoCorr = parset.getBool(prefix + "applyautocorr", false);
  itsMinBLength = parset.getDouble(prefix + "blmin", -1);
  itsMaxBLength = parset.getDouble(prefix + "blmax", 1e30);
}

MadFlagger::~MadFlagger() {}

void MadFlagger::show(std::ostream& os) const {
  os << "MADFlagger " << itsName << '\n';
  os << "  freqwindow:     " << itsFreqWindowStr
     << "   (max = " << itsFreqWindow << ')' << '\n';
  os << "  timewindow:     " << itsTimeWindowStr
     << "   (max = " << itsTimeWindow << ')' << '\n';
  os << "  threshold:      " << itsThresholdStr << "   (max = " << itsThreshold
     << ')' << '\n';
  os << "  correlations:   " << itsFlagCorr << '\n';
  os << "  applyautocorr:  " << itsApplyAutoCorr
     << "   (nautocorr = " << itsNrAutoCorr << ')' << '\n';
  os << "  blmin:          " << itsMinBLength << " m" << '\n';
  os << "  blmax:          " << itsMaxBLength << " m" << '\n';
}

void MadFlagger::showCounts(std::ostream& os) const {
  os << "\nFlags set by MADFlagger " << itsName;
  os << "\n=======================\n";
  itsFlagCounter.showBaseline(os, itsNTimes);
  itsFlagCounter.showChannel(os, itsNTimes);
  itsFlagCounter.showCorrelation(os, itsNTimes);
}

void MadFlagger::showTimings(std::ostream& os, double duration) const {
  double flagDur = itsTimer.getElapsed();
  os << "  ";
  base::FlagCounter::showPerc1(os, flagDur, duration);
  os << " MADFlagger " << itsName << '\n';
  os << "          ";
  // move time and median time are sum of all threads.
  // Scale them to a single elapsed time.
  double factor =
      (itsComputeTimer.getElapsed() / (itsMoveTime + itsMedianTime));
  base::FlagCounter::showPerc1(os, itsMoveTime * factor, flagDur);
  os << " of it spent in shuffling data\n";
  os << "          ";
  base::FlagCounter::showPerc1(os, itsMedianTime * factor, flagDur);
  os << " of it spent in calculating medians\n";
}

void MadFlagger::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  // Get baseline indices of autocorrelations.
  itsAutoCorrIndex = info().getAutoCorrIndex();
  itsNrAutoCorr = 0;
  for (unsigned int i = 0; i < itsAutoCorrIndex.size(); ++i) {
    if (itsAutoCorrIndex[i] >= 0) {
      itsNrAutoCorr++;
    }
  }
  if (itsApplyAutoCorr && itsNrAutoCorr <= 0)
    throw std::runtime_error(
        "applyautocorr=True cannot be used if "
        "the data does not contain autocorrelations");
  // Calculate the baseline lengths.
  itsBLengths = info().getBaselineLengths();
  // Evaluate the window size expressions.
  getExprValues(infoIn.nchan(), infoIn.ntime());
  itsBuffers.resize(itsTimeWindow);
  itsAmplitudes.resize(itsTimeWindow);
  for (size_t i = 0; i < itsAmplitudes.size(); ++i) {
    itsAmplitudes[i].resize(
        {infoIn.nbaselines(), infoIn.nchan(), infoIn.ncorr()});
  }
  // Set or check the correlations to flag on.
  std::vector<unsigned int> flagCorr;
  unsigned int ncorr = infoIn.ncorr();
  if (itsFlagCorr.empty()) {
    // No correlations given means use them all.
    for (unsigned int i = 0; i < ncorr; ++i) {
      flagCorr.push_back(i);
    }
  } else {
    for (const unsigned int ip : itsFlagCorr) {
      // Only take valid corrrelations.
      if (itsFlagCorr[ip] < ncorr) {
        flagCorr.push_back(itsFlagCorr[ip]);
      }
    }
    // If no valid left, use first one.
    if (flagCorr.empty()) {
      DPLOG_INFO_STR("No valid correlations given in MadFlagger " + itsName +
                     "; first one will be used");
      flagCorr.push_back(0);
    }
  }
  itsFlagCorr = flagCorr;
  // Initialize the flag counters.
  itsFlagCounter.init(getInfo());
}

bool MadFlagger::process(std::unique_ptr<base::DPBuffer> buffer) {
  itsTimer.start();
  // Accumulate buffers in the time window.
  // The buffers are wrapped, thus oldest entries are overwritten.
  unsigned int index = itsNTimes % itsTimeWindow;
  itsBuffers[index] = std::move(buffer);
  base::DPBuffer& dBuffer = *itsBuffers[index];
  // Calculate amplitudes if needed.
  itsAmplitudes[index] = xt::abs(dBuffer.GetData());
  // Fill flags if needed.
  if (dBuffer.GetFlags().size() == 0) {
    dBuffer.ResizeFlags(dBuffer.GetData().shape());
    dBuffer.GetFlags().fill(false);
  }
  itsNTimes++;

  // Flag if there are enough time entries accumulated in itsBuffers.
  if (itsNTimes > itsTimeWindow / 2) {
    // Fill the vector telling which time entries to use for the medians.
    // If window is not entirely full, use copies as needed.
    // This is done as follows:
    // Suppose timewindow=9 and we only have entries 0,1,2,3,4.
    // The entries are mirrored, thus we get 4,3,2,1,0,1,2,3,4
    // to obtain sufficient time entries.
    // Furthermore, indices are arranged such that any width can be used; thus
    // first center, then one left and right, second left and right, etc.
    // So,        [4,3,2,1,0,1,2,3,4] (with center 0) yields [0,1,1,2,2,3,3,4,4]
    // Similarly, [0,1,2,3,4,5,6,7,8] (with center 4) yields [4,3,5,2,6,1,7,0,8]
    // and        [2,3,4,5,6,7,8,0,1] (with center 6) yields [6,5,7,4,8,3,0,2,1]
    std::vector<unsigned int> timeEntries;
    timeEntries.reserve(itsTimeWindow);
    timeEntries.push_back(itsNTimesDone % itsTimeWindow);  // center
    for (unsigned int i = 1; i <= itsTimeWindow / 2; ++i) {
      timeEntries.push_back(std::abs(int(itsNTimesDone) - int(i)) %
                            itsTimeWindow);
      timeEntries.push_back((itsNTimesDone + i) % itsTimeWindow);
    }
    flag(itsNTimesDone % itsTimeWindow, timeEntries);
    itsNTimesDone++;
  }
  itsTimer.stop();
  return true;
}

void MadFlagger::finish() {
  itsTimer.start();
  // Adjust window size if there are fewer time entries.
  if (itsNTimes < itsTimeWindow) {
    itsTimeWindow = 1 + ((itsNTimes - 1) / 2) * 2;  // make sure it is odd
  }
  unsigned int halfWindow = itsTimeWindow / 2;
  std::vector<unsigned int> timeEntries(itsTimeWindow);
  // Process possible leading entries.
  // This can happen if the window was larger than number of times.
  while (itsNTimesDone < itsNTimes - halfWindow) {
    // Process in the same way as in process.
    unsigned int inx = 0;
    timeEntries[inx++] = itsNTimesDone % itsTimeWindow;  // center
    for (unsigned int i = 1; i <= halfWindow; ++i) {
      timeEntries[inx++] =
          std::abs(int(itsNTimesDone) - int(i)) % itsTimeWindow;
      timeEntries[inx++] = (itsNTimesDone + i) % itsTimeWindow;
    }
    flag(itsNTimesDone, timeEntries);
    itsNTimesDone++;
  }
  assert(itsNTimes - itsNTimesDone == halfWindow);
  // Process the remaining time entries.
  while (itsNTimesDone < itsNTimes) {
    unsigned int inx = 0;
    timeEntries[inx++] = itsNTimesDone % itsTimeWindow;  // center
    for (unsigned int i = 1; i <= halfWindow; ++i) {
      timeEntries[inx++] =
          std::abs(int(itsNTimesDone) - int(i)) % itsTimeWindow;
      // Time entries might need to be mirrored at the end.
      unsigned int ri = itsNTimesDone + i;
      if (ri >= itsNTimes) {
        ri = 2 * (itsNTimes - 1) - ri;
      }
      timeEntries[inx++] = ri % itsTimeWindow;
    }
    flag(itsNTimesDone % itsTimeWindow, timeEntries);
    itsNTimesDone++;
  }
  itsTimer.stop();
  // Let the next step finish its processing.
  getNextStep()->finish();
}

void MadFlagger::flagBaseline(
    const std::vector<int>& ant1, const std::vector<int>& ant2,
    const std::vector<unsigned int>& timeEntries, unsigned int ib,
    unsigned int ncorr, unsigned int nchan, const float* bufferDataPtr,
    bool* bufferFlagPtr, float& Z1, float& Z2, std::vector<float>& tempBuffer,
    base::FlagCounter& counter, common::NSTimer& moveTimer,
    common::NSTimer& medianTimer) {
  unsigned int blsize = ncorr * nchan;
  float MAD = 1.4826;  ///< constant determined by Pandey

  double threshold = itsThresholdArr[ib];
  // Do only autocorrelations if told so.
  // Otherwise do baseline only if length within min-max.
  if ((!itsApplyAutoCorr && itsBLengths[ib] >= itsMinBLength &&
       itsBLengths[ib] <= itsMaxBLength) ||
      (itsApplyAutoCorr && ant1[ib] == ant2[ib])) {
    for (unsigned int ic = 0; ic < nchan; ++ic) {
      size_t offset = ib * blsize + ic * ncorr;
      bool* flagPtr = &bufferFlagPtr[offset];
      const float* dataPtr = &bufferDataPtr[offset];
      bool corrIsFlagged = false;

      // Iterate over given correlations.
      for (unsigned int i = 0; i < itsFlagCorr.size(); i++) {
        unsigned int ip = itsFlagCorr[i];
        // If one correlation is flagged, all of them will be flagged.
        // So no need to check others.
        if (flagPtr[ip]) {
          corrIsFlagged = true;
          break;
        }
        // Calculate values from the median.
        computeFactors(timeEntries, ib, ic, ip, nchan, ncorr, Z1, Z2,
                       tempBuffer, moveTimer,
                       medianTimer);  // main outputs: Z1 & Z2
        if (dataPtr[ip] > Z1 + threshold * Z2 * MAD) {
          corrIsFlagged = true;
          counter.incrBaseline(ib);
          counter.incrChannel(ic);
          counter.incrCorrelation(ip);
          break;
        }
      }
      if (corrIsFlagged) {
        for (unsigned int ip = 0; ip < ncorr; ++ip) {
          flagPtr[ip] = true;
        }
      }
    }
  }
}

void MadFlagger::flag(unsigned int index,
                      const std::vector<unsigned int>& timeEntries) {
  // Get antenna numbers in case applyautocorr is true.
  const std::vector<int>& ant1 = getInfo().getAnt1();
  const std::vector<int>& ant2 = getInfo().getAnt2();
  // MadFlagger requires non-flagged DPBuffers of multiple time steps, of which
  // it only needs the amplitudes of the data (visibilities) and the flags. The
  // amplitudes are stored separately in itsAmplitudes.
  std::unique_ptr<DPBuffer> bufferForNextStep = std::move(itsBuffers[index]);
  itsBuffers[index] =
      std::make_unique<DPBuffer>(*bufferForNextStep, kFlagsField);

  const std::array<size_t, 3> shape = bufferForNextStep->GetData().shape();
  unsigned int ncorr = shape[2];
  unsigned int nchan = shape[1];
  unsigned int blsize = ncorr * nchan;
  int nrbl = shape[0];
  // Get pointers to amplitudes (of visibility data) and flags.
  const float* bufferDataPtr = itsAmplitudes[index].data();
  bool* bufferFlagPtr = bufferForNextStep->GetFlags().data();
  itsComputeTimer.start();
  // Now flag each baseline, channel and correlation for this time window.
  // This can be done in parallel.
  struct ThreadData {
    std::vector<float> tempBuffer;
    base::FlagCounter counter;
    common::NSTimer moveTimer;
    common::NSTimer medianTimer;
    float Z1;  // median
    float Z2;  // median of absolute difference
  };
  std::vector<ThreadData> threadData(getInfo().nThreads());

  // Create thread-private counter.
  for (ThreadData& data : threadData) {
    data.tempBuffer.resize(itsFreqWindow * timeEntries.size());
    data.counter.init(getInfo());
  }

  // The for loop can be parallelized. This must be done dynamically,
  // because the execution time of each iteration can vary a lot.
  aocommon::ParallelFor<size_t> loop(getInfo().nThreads());
  loop.Run(0, nrbl, [&](size_t ib, size_t thread) {
    ThreadData& data = threadData[thread];
    flagBaseline(ant1, ant2, timeEntries, ib, ncorr, nchan, bufferDataPtr,
                 bufferFlagPtr, data.Z1, data.Z2, data.tempBuffer, data.counter,
                 data.moveTimer, data.medianTimer);
  });  // end of parallel loop

  // Add the counters to the overall object.
  for (ThreadData& data : threadData) {
    itsFlagCounter.add(data.counter);
    // Add the timings.
    itsMoveTime += data.moveTimer.getElapsed();
    itsMedianTime += data.medianTimer.getElapsed();
  }

  itsComputeTimer.stop();
  // Apply autocorrelations flags if needed.
  // Only to baselines with length within min-max.
  if (itsApplyAutoCorr) {
    for (int ib = 0; ib < nrbl; ++ib) {
      // Flag crosscorr if at least one autocorr is present.
      // Only if baseline length within min-max.
      if (ant1[ib] != ant2[ib] && itsBLengths[ib] >= itsMinBLength &&
          itsBLengths[ib] <= itsMaxBLength) {
        int inx1 = itsAutoCorrIndex[ant1[ib]];
        int inx2 = itsAutoCorrIndex[ant2[ib]];
        if (inx1 >= 0 || inx2 >= 0) {
          // Find flags of the autocorrelations of both antennae.
          // Use other autocorr if one autocorr does not exist.
          // In this way inx does not need to be tested in the inner loop.
          if (inx1 < 0) {
            inx1 = inx2;
          } else if (inx2 < 0) {
            inx2 = inx1;
          }
          // Flag if not flagged yet and if one of autocorr is flagged.
          for (unsigned int ic = 0; ic < nchan; ++ic) {
            bool* flagPtr = &bufferFlagPtr[ib * blsize + ic * ncorr];
            bool* flagAnt1 = &bufferFlagPtr[inx1 * nchan * ncorr + ic * ncorr];
            bool* flagAnt2 = &bufferFlagPtr[inx2 * nchan * ncorr + ic * ncorr];
            if (!*flagPtr && (*flagAnt1 || *flagAnt2)) {
              for (unsigned int ip = 0; ip < ncorr; ++ip) {
                flagPtr[ip] = true;
              }
              itsFlagCounter.incrBaseline(ib);
              itsFlagCounter.incrChannel(ic);
            }
          }
        }
      }
    }
  }
  // Process the result in the next step.
  itsTimer.stop();
  getNextStep()->process(std::move(bufferForNextStep));
  itsTimer.start();
}

void MadFlagger::computeFactors(const std::vector<unsigned int>& timeEntries,
                                unsigned int bl, int chan, int corr, int nchan,
                                int ncorr, float& Z1, float& Z2,
                                std::vector<float>& tempBuffer,
                                common::NSTimer& moveTimer,
                                common::NSTimer& medianTimer) {
  moveTimer.start();
  // Collect all non-flagged data points for given baseline, channel,
  // and correlation in the window around the channel.
  unsigned int np = 0;
  // At the beginning or end of the window the values are wrapped.
  // So we might need to move in two parts.
  // This little piece of code is tested in tMirror.cc.
  int hw = itsFreqWindowArr[bl] / 2;
  int s1 = chan - hw;
  int e1 = chan + hw + 1;
  int s2 = 1;
  int e2 = 1;
  if (s1 < 0) {
    e2 = -s1 + 1;
    s1 = 0;
  } else if (e1 > nchan) {
    s2 = nchan + nchan - e1 - 1;  // e1-nchan+1 too far, go back that amount
    e2 = nchan - 1;
    e1 = nchan;
  }
  // Iterate over all time entries.
  const unsigned int* iter = &(timeEntries[0]);
  const unsigned int* endIter = iter + itsTimeWindowArr[bl];
  for (; iter != endIter; ++iter) {
    const base::DPBuffer& inbuf = *itsBuffers[*iter];
    const xt::xtensor<float, 3>& ampl = itsAmplitudes[*iter];
    // Get pointers to given baseline and correlation.
    unsigned int offset = bl * nchan * ncorr + corr;
    const float* dataPtr = ampl.data() + offset;
    const bool* flagPtr = inbuf.GetFlags().data() + offset;
    // Now move data from the two channel parts.
    for (int i = s1 * ncorr; i < e1 * ncorr; i += ncorr) {
      if (!flagPtr[i]) {
        tempBuffer[np++] = dataPtr[i];
      }
    }
    for (int i = s2 * ncorr; i < e2 * ncorr; i += ncorr) {
      if (!flagPtr[i]) {
        tempBuffer[np++] = dataPtr[i];
      }
    }
  }
  moveTimer.stop();
  // If only flagged data, don't do anything.
  if (np == 0) {
    Z1 = -1.0;
    Z2 = 0.0;
  } else {
    medianTimer.start();
    // Get median of data and get median of absolute difference.
    Z1 = dp3::common::Median(tempBuffer);
    for (unsigned int i = 0; i < np; ++i) {
      tempBuffer[i] = std::abs(tempBuffer[i] - Z1);
    }
    Z2 = dp3::common::Median(tempBuffer);
    medianTimer.stop();
  }
}

void MadFlagger::getExprValues(int maxNChan, int maxNTime) {
  // Parse the expressions.
  // Baseline length can be used as 'bl' in the expressions.
  Record rec;
  rec.define("bl", double(0));
  TableExprNode node1(RecordGram::parse(rec, itsFreqWindowStr));
  TableExprNode node2(RecordGram::parse(rec, itsTimeWindowStr));
  TableExprNode node3(RecordGram::parse(rec, itsThresholdStr));
  // Size the arrays.
  size_t nrbl = itsBLengths.size();
  itsThresholdArr.reserve(nrbl);
  itsTimeWindowArr.reserve(nrbl);
  itsFreqWindowArr.reserve(nrbl);
  itsFreqWindow = 0;
  itsTimeWindow = 0;
  itsThreshold = -1e30;
  // Evaluate the expression for each baseline.
  double result;
  RecordFieldPtr<double> blref(rec, "bl");
  for (unsigned int i = 0; i < nrbl; ++i) {
    // Put the length of each baseline in the record used to evaluate.
    *blref = itsBLengths[i];
    // Evaluate freqwindow size and make it odd if needed.
    node1.get(rec, result);
    int freqWindow = std::min(std::max(1, int(result + 0.5)), maxNChan);
    if (freqWindow % 2 == 0) {
      freqWindow--;
    }
    itsFreqWindowArr.push_back(freqWindow);
    itsFreqWindow = std::max(itsFreqWindow, (unsigned int)(freqWindow));
    // Evaluate timewindow size and make it odd if needed.
    node2.get(rec, result);
    int timeWindow = std::max(1, int(result + 0.5));
    if (maxNTime > 0 && timeWindow > maxNTime) {
      timeWindow = maxNTime;
    }
    if (timeWindow % 2 == 0) {
      timeWindow--;
    }
    itsTimeWindowArr.push_back(timeWindow);
    itsTimeWindow = std::max(itsTimeWindow, (unsigned int)(timeWindow));
    // Evaluate threshold.
    node3.get(rec, result);
    itsThresholdArr.push_back(result);
    if (result > itsThreshold) {
      itsThreshold = result;
    }
  }
}

}  // namespace steps
}  // namespace dp3
