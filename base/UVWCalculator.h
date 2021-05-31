// UVWCalculator.h: Class to calculate UVW coordinates
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to calculate UVW coordinates
/// Note: this code is used by LOFAR and APERTIF software.
/// @author Ger van Diepen

#ifndef DPPP_UVWCALCULATOR_H
#define DPPP_UVWCALCULATOR_H

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MBaseline.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MCBaseline.h>
#include <casacore/casa/Arrays/Vector.h>

namespace dp3 {
namespace base {

/// @brief Class to calculate UVW coordinates

/// This class calculates the UVW coordinates for a given baseline and
/// time stamp in the same way as done in LofarStMan.
///
/// It calculates and caches the UVW coordinates per antenna and combines
/// them to get the baseline UVW coordinates. This is much faster than
/// calculating baseline UVW coordinates directly.

class UVWCalculator {
 public:
  /// Construct the object for the given phase direction, array position,
  /// and station positions.
  UVWCalculator(const casacore::MDirection& phaseDir,
                const casacore::MPosition& arrayPosition,
                const std::vector<casacore::MPosition>& stationPositions);

  /// get the UVW coordinates for the given baseline and time.
  std::array<double, 3> getUVW(unsigned int ant1, unsigned int ant2,
                               double time);

 private:
  casacore::MDirection itsPhaseDir;
  bool itsMovingPhaseDir;
  casacore::MDirection::Convert itsDirToJ2000;  ///< direction to J2000
  casacore::MBaseline::Convert itsBLToJ2000;    ///< convert ITRF to J2000
  casacore::MeasFrame itsFrame;
  std::vector<casacore::MBaseline> itsAntMB;
  std::vector<std::array<double, 3>> itsAntUvw;
  casacore::Block<bool> itsUvwFilled;
  double itsLastTime;
};

}  // namespace base
}  // namespace dp3

#endif
