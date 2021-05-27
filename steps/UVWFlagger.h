// UVWFlagger.h: DPPP step class to flag data on UVW coordinates
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class to flag data on UVW coordinates
/// @author Ger van Diepen

#ifndef DPPP_UVWFLAGGER_H
#define DPPP_UVWFLAGGER_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/UVWCalculator.h"

namespace dp3 {
namespace common {
class ParameterSet;
class ParameterValue;
}  // namespace common

namespace steps {
/// @brief DPPP step class to flag data on UVW coordinates

/// This class is a Step class flagging data points based on data
/// selections given in the parset file.
/// The following selections can be given:
/// <ul>
///  <li> minimum and/or maximum UV distance
///  <li> minimum or maximum value for U, V, and/or W
///  <li> both can be used with a different phase center which can be
///       be given as a position or as a moving source like SUN or JUPITER.
/// </ul>
/// The UVW values can be given in meters or in wavelengths.

class UVWFlagger : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  /// The antenna names are used to find antenna numbers.
  /// The channel frequencies as they are in the input step must be given
  /// starting at the start-channel.
  UVWFlagger(InputStep*, const common::ParameterSet&, const string& prefix);

  virtual ~UVWFlagger();

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

  /// Show the flag counts.
  virtual void showCounts(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  /// Test if uvw matches a range in meters.
  bool testUVWm(double uvw, const std::vector<double>& ranges);

  /// Set flags for channels where uvw (in m) matches a range in wavelengths.
  void testUVWl(double uvw, const std::vector<double>& ranges, bool* flagPtr,
                unsigned int nrcorr);

  /// Return a vector with UVW ranges.
  /// It looks for the named parameter suffixed with 'range', 'min', and
  /// 'max'. The returned vector contains 2 subsequent values for each range
  /// (min and max are also turned into a range).
  /// Optionally the values are squared to avoid having to take a sqrt
  /// of the data's UVW coordinates.
  std::vector<double> fillUVW(const common::ParameterSet& parset,
                              const string& prefix, const string& name,
                              bool square);

  /// Handle the specification of a phase center.
  /// It setups the UVWCalculator.
  void handleCenter();

  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBuffer;
  unsigned int itsNTimes;

  std::vector<double> itsRecWavel;        ///< reciprokes of wavelengths
  const std::vector<double> itsRangeUVm;  ///< UV ranges (in m) to be flagged
  const std::vector<double> itsRangeUm;   ///< U  ranges (in m) to be flagged
  const std::vector<double> itsRangeVm;   ///< V  ranges (in m) to be flagged
  const std::vector<double> itsRangeWm;   ///< W  ranges (in m) to be flagged
  const std::vector<double>
      itsRangeUVl;                       ///< UV ranges (in wavel) to be flagged
  const std::vector<double> itsRangeUl;  ///< U  ranges (in wavel) to be flagged
  const std::vector<double> itsRangeVl;  ///< V  ranges (in wavel) to be flagged
  const std::vector<double> itsRangeWl;  ///< W  ranges (in wavel) to be flagged

  /// If nothing is filled in, this step just passes through data
  const bool itsIsDegenerate;

  std::unique_ptr<base::UVWCalculator> itsUVWCalc;
  const std::vector<string> itsCenter;
  common::NSTimer itsTimer;
  common::NSTimer itsUVWTimer;
  base::FlagCounter itsFlagCounter;
};

}  // namespace steps
}  // namespace dp3

#endif
