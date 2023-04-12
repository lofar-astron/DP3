// StationAdder.h: DP3 step class to add stations as a superstation
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to add stations as a superstation
/// @author Ger van Diepen

#ifndef DP3_STEPS_STATIONADDER_H_
#define DP3_STEPS_STATIONADDER_H_

#include <casacore/measures/Measures/MPosition.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/steps/Step.h>

#include "../base/UVWCalculator.h"

#include "../common/ParameterRecord.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {
/// @brief DP3 step class to add stations as a superstation

/// This class is a Step class summing stations to a superstation.
///
/// It is possible to define one or more groups of stations to be summed.
/// Each group has a name which is the name of the new station.
/// The complex values of baselines are added for which one station occurs
/// in only one group. A baseline is not added if no or both stations are
/// member of a summing group.
/// <br>The summation is done in a weighted way, where the weight of a
/// new station is the sum of the original weights. Optionally weights 1
/// can be used instead of the original weights.
///
/// Only unflagged data points are used. If too few data points are
/// unflagged, the output data point is flagged.
///
/// Questions:
/// 1. check if phases do not differ too much? Flag if too much?
/// 2. must all stations exist or possible that some don't?

class StationAdder : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  StationAdder(const common::ParameterSet&, const std::string& prefix);

  ~StationAdder() override;

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  common::Fields getProvidedFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Add new meta info to the MS.
  void addToMS(const string& msName) override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Return the indices of the stations in antennaNames matching
  /// the pattern list.
  /// The patterns are processed from left to right. A pattern can start
  /// with ! or ^ meaning that the the matches are discarded. In this
  /// way first a broad pattern can be given, which can be narrowed down.
  /// A warning is given if a pattern does not match any station name.
  static std::vector<int> getMatchingStations(
      const std::vector<std::string>& antennaNames,
      const std::vector<string>& patterns);

 private:
  /// Update the beam info subtables.
  void updateBeamInfo(const string& msName, unsigned int origNant,
                      casacore::Table& antTab);

  std::string itsName;
  base::DPBuffer itsBuf;
  common::ParameterRecord itsStatRec;  ///< stations definitions
  std::vector<casacore::Vector<int>>
      itsParts;  ///< the stations in each superstation
  std::vector<std::vector<int>>
      itsBufRows;             ///< old baseline rows in each new baseline
  unsigned int itsMinNPoint;  ///< flag data if too few unflagged data
  bool itsMakeAutoCorr;       ///< also form new auto-correlations?
  bool itsSumAutoCorr;        ///< sum auto- or cross-correlations?
  bool itsDoAverage;          ///< average or sum?
  bool itsUseWeight;          ///< false = use weight 1 per station
  std::unique_ptr<base::UVWCalculator> itsUVWCalc;
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
