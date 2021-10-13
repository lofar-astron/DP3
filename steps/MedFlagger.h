// MedFlagger.h: DPPP step class to flag data based on median filtering
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to flag data based on median filtering
/// @author Ger van Diepen

#ifndef DPPP_MEDFLAGGER_H
#define DPPP_MEDFLAGGER_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/FlagCounter.h"

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DPPP step class to flag data based on median filtering

/// This class is a Step class flagging data points based on the median
/// of the absolute difference of the data and the median of the data.
/// Both medians are taken in a time/frequency window around the data point.
/// Only unflagged data points in the window are taken into account.
/// The size of the window is given in the parset file.
///
/// The window around data points at the edges is formed by mirroring the
/// data at the edge. For example, for channel 0 and a window size of 7
/// the data are mirrored, thus channels 3,2,1,0,1,2,3 will be used.
/// For channel 1 the channels 2,1,0,1,2,3,4 will be used.
/// The test program tMirror.cc can be used to check the correctness of
/// the alogorithm to determine the channels to use.
///
/// Taking the median is an O(N) operation, thus doing it for all data
/// points is an O(N^2) operation. The test program tMedian.cc can be
/// used to test the performance of the algorithms to determine the median.
/// It shows that casacore's kthLargest outperforms STL's nth_element.
/// <br>
/// Shuffling the data around to be able to determine the medians is also
/// an expensive operation, but takes less time than the medians themselves.
///
/// When a correlation is flagged, all correlations for that data point
/// are flagged. It is possible to specify which correlations have to be
/// taken into account when flagging. Using, say, only XX may boost
/// performance with a factor 4, but miss points to be flagged.
/// It is also possible to specify the order in which the correlations
/// have to be tested.
///
/// It is possible to flag specific baselines only using a selection on
/// baseline length.
/// <br>Furthermore it is possible to only flag the autocorrelations and
/// apply the resulting flags to the crosscorrelations, possibly selected
/// on baseline length.

class MedFlagger : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  MedFlagger(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~MedFlagger();

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  /// It is used to adjust the parms if needed.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the flagger counts.
  virtual void showCounts(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Flag for the entry at the given index.
  /// Use the given time entries for the medians.
  /// Process the result in the next step.
  void flag(unsigned int index, const std::vector<unsigned int>& timeEntries);

  /// Compute the median factors for given baseline, channel, and
  /// correlation.
  void computeFactors(const std::vector<unsigned int>& timeEntries,
                      unsigned int bl, int chan, int corr, int nchan, int ncorr,
                      float& Z1, float& Z2, float* tempBuf,
                      common::NSTimer& moveTimer, common::NSTimer& medianTimer);

  /// Get the values of the expressions for each baseline.
  void getExprValues(int maxNChan, int maxNTime);

 protected:
  InputStep* itsInput;
  string itsName;
  string itsThresholdStr;
  string itsFreqWindowStr;
  string itsTimeWindowStr;
  std::vector<float> itsThresholdArr;  ///< threshold per baseline
  std::vector<unsigned int>
      itsFreqWindowArr;  ///< freq window size per baseline
  std::vector<unsigned int>
      itsTimeWindowArr;  ///< time window size per baseline
  float itsThreshold;
  unsigned int itsFreqWindow;
  unsigned int itsTimeWindow;
  unsigned int itsNTimes;
  unsigned int itsNTimesDone;
  std::vector<unsigned int> itsFlagCorr;
  bool itsApplyAutoCorr;
  std::vector<int> itsAutoCorrIndex;  ///< baseline index of autocorrelations
  unsigned int itsNrAutoCorr;
  double itsMinBLength;            ///< minimum baseline length
  double itsMaxBLength;            ///< maximum baseline length
  std::vector<double> itsBLength;  ///< length of each baseline
  std::vector<base::DPBuffer> itsBuf;
  std::vector<casacore::Cube<float>> itsAmpl;  ///< amplitudes of the data
  base::FlagCounter itsFlagCounter;
  common::NSTimer itsTimer;
  common::NSTimer itsComputeTimer;  ///< move/median timer
  double itsMoveTime;               ///< data move timer (sum all threads)
  double itsMedianTime;             ///< median timer (sum of all threads)
};

}  // namespace steps
}  // namespace dp3

#endif
