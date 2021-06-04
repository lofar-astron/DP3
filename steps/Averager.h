// Averager.h: DPPP step class to average in time and/or freq
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to average in time and/or freq
/// @author Ger van Diepen

#ifndef DPPP_AVERAGER_H
#define DPPP_AVERAGER_H

#include "InputStep.h"

#include "../base/DPBuffer.h"

#include <casacore/casa/Arrays/Cube.h>

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

/// \brief DPPP step class to average in time and/or freq

/// This class is a Step class calculating the weighted average of
/// data in time and/or frequency.
/// <br>
/// Only unflagged data points are used. The average is calculated as
/// <tt>sum(data*weight) / sum(weight)</tt> and the sum of the weights
/// is the weight of the new data point. If all data point to use are
/// flagged, the resulting data point and weight are set to zero and flagged.
//
/// It keeps track of the FullResFlags. It sets them if the corresponding
/// data point is flagged. Note that multiple FullResFlags elements map to
/// a single data point if some averaging was done before.

class Averager : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Averager(InputStep&, const common::ParameterSet&, const string& prefix);

  /// Construct the object using the given parameters.
  Averager(InputStep&, const string& stepname, unsigned int nchanAvg,
           unsigned int ntimeAvg);

  virtual ~Averager();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  /// Average into itsBufOut.
  void average();

  /// Copy the fullRes flags in the input buffer to the correct
  /// time index in the output buffer.
  /// If a flag is set, set all flags in corresponding FullRes window.
  void copyFullResFlags(const casacore::Cube<bool>& fullResFlags,
                        const casacore::Cube<bool>& flags, int timeIndex);

  /// Get the value in Hertz of a string like "1000 MHz". If unit is
  /// omitted it defaults to Hertz
  double getFreqHz(const string& freqstr);

  InputStep& itsInput;
  std::string itsName;
  base::DPBuffer itsBuf;
  base::DPBuffer itsBufTmp;
  base::DPBuffer itsBufOut;
  casacore::Cube<int> itsNPoints;
  casacore::Cube<casacore::Complex> itsAvgAll;
  casacore::Cube<float> itsWeightAll;
  casacore::Cube<bool> itsFullResFlags;
  double itsFreqResolution;
  double itsTimeResolution;
  unsigned int itsNChanAvg;
  unsigned int itsNTimeAvg;
  unsigned int itsMinNPoint;
  float itsMinPerc;
  unsigned int itsNTimes;
  double itsTimeInterval;
  bool itsNoAvg;  ///< No averaging (i.e. both 1)?
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
