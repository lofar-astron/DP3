// Averager.h: DP3 step class to average in time and/or freq
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to average in time and/or freq
/// @author Ger van Diepen

#ifndef DP3_STEPS_AVERAGER_H_
#define DP3_STEPS_AVERAGER_H_

#include <memory>

#include <xtensor/xtensor.hpp>
#include <aocommon/staticfor.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/steps/Step.h>
#include "../common/Timer.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

/// \brief DP3 step class to average in time and/or freq

/// This class is a Step class calculating the weighted average of
/// data in time and/or frequency.
/// <br>
/// Only unflagged data points are used. The average is calculated as
/// <tt>sum(data*weight) / sum(weight)</tt> and the sum of the weights
/// is the weight of the new data point. If all data point to use are
/// flagged, the resulting data point and weight are set to zero and flagged.

class Averager : public Step {
 public:
  /// Required and provided fields are fixed. These constants simplify
  /// get*Fields implementations for steps that have Averagers as sub-step.
  static const common::Fields kRequiredFields;
  static const common::Fields kProvidedFields;

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Averager(const common::ParameterSet&, const string& prefix);

  /// Construct the object using the given parameters.
  Averager(const string& stepname, unsigned int nchanAvg,
           unsigned int ntimeAvg);
  Averager(const string& stepname, double freq_resolution,
           double time_resolution);

  ~Averager() override;

  common::Fields getRequiredFields() const override { return kRequiredFields; }

  common::Fields getProvidedFields() const override { return kProvidedFields; }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Get the value in Hertz of a string like "1000 MHz". If unit is
  /// omitted it defaults to Hertz
  static double getFreqHz(const string& freqstr);

 private:
  /// Update itsBuf so it contains averages.
  void average();

  std::string itsName;
  std::unique_ptr<base::DPBuffer> itsBuf;
  xt::xtensor<int, 3> itsNPoints;
  xt::xtensor<std::complex<float>, 3> itsAvgAll;
  xt::xtensor<float, 3> itsWeightAll;
  double itsFreqResolution;
  double itsTimeResolution;
  unsigned int itsNChanAvg;
  unsigned int itsNTimeAvg;
  unsigned int itsMinNPoint;
  float itsMinPerc;
  unsigned int itsNTimes;
  double itsOriginalTimeInterval;
  bool itsNoAvg;  ///< No averaging (i.e. both 1)?
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
