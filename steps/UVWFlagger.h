// UVWFlagger.h: DP3 step class to flag data on UVW coordinates
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to flag data on UVW coordinates
/// @author Ger van Diepen

#ifndef DP3_STEPS_UVWFLAGGER_H_
#define DP3_STEPS_UVWFLAGGER_H_

#include <dp3/base/BdaBuffer.h>
#include <dp3/base/DPBuffer.h>

#include "InputStep.h"

#include "../base/UVWCalculator.h"

namespace dp3 {
namespace common {
class ParameterSet;
class ParameterValue;
}  // namespace common

namespace steps {
/// @brief DP3 step class to flag data on UVW coordinates

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
  UVWFlagger(const common::ParameterSet&, const std::string& prefix,
             MsType inputType);

  ~UVWFlagger() override;

  common::Fields getRequiredFields() const override {
    if (itsIsDegenerate) return common::Fields();

    common::Fields fields = kFlagsField;

    if (itsCenter.empty()) fields |= kUvwField;

    return fields;
  }

  common::Fields getProvidedFields() const override {
    return itsIsDegenerate ? common::Fields() : kFlagsField;
  }

  /// Process the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  bool process(std::unique_ptr<base::BdaBuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  /// It is used to adjust the parms if needed.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the flag counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  bool accepts(MsType dt) const override { return dt == itsInputType; }

  MsType outputs() const override { return itsInputType; }

  bool isDegenerate() const { return itsIsDegenerate; }

 private:
  /// Test if uvw matches a range in meters.
  bool testUVWm(double uvw, const std::vector<double>& ranges);

  /// Set flags for channels where uvw (in m) matches a range in wavelengths.
  void testUVWl(double uvw, const std::vector<double>& ranges, bool* flagPtr,
                unsigned int nrcorr, unsigned int baseline_id);

  // Apply UVW flagging
  void doFlag(const std::array<double, 3>& uvw, bool* flagPtr,
              unsigned int n_correlations, unsigned int n_channels,
              unsigned int baseline_id = 0);

  /// Return a vector with UVW ranges.
  /// It looks for the named parameter suffixed with 'range', 'min', and
  /// 'max'. The returned vector contains 2 subsequent values for each range
  /// (min and max are also turned into a range).
  /// Optionally the values are squared to avoid having to take a sqrt
  /// of the data's UVW coordinates.
  std::vector<double> fillUVW(const common::ParameterSet& parset,
                              const std::string& prefix,
                              const std::string& name, bool square);

  /// Handle the specification of a phase center.
  /// It setups the UVWCalculator.
  void handleCenter();

  const MsType itsInputType;
  std::string itsName;
  unsigned int itsNTimes;

  /// Wavelengths may differ between baselines when using BDA.
  /// If all baselines have equal channels (i.e. when the input is a DPBuffer),
  /// the outer vector only holds one inner vector. When using BDA (i.e. when
  /// the input is a BdaBuffer), each baseline has its own inner vector.
  std::vector<std::vector<double>> itsRecWavel;  ///< reciprocals of wavelengths
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
  const std::vector<std::string> itsCenter;
  common::NSTimer itsTimer;
  common::NSTimer itsUVWTimer;
  base::FlagCounter itsFlagCounter;
};

}  // namespace steps
}  // namespace dp3

#endif
