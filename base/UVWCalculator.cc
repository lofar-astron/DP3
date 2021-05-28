// UVWCalculator.cc: Class to calculate UVW coordinates
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

// Note: this code is used by LOFAR and APERTIF software.

#include "UVWCalculator.h"

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/Muvw.h>

using namespace casacore;

namespace dp3 {
namespace base {

UVWCalculator::UVWCalculator(const MDirection& phaseDir,
                             const MPosition& arrayPos,
                             const vector<MPosition>& stationPositions) {
  // Convert the station positions to a baseline in ITRF.
  int nrant = stationPositions.size();
  Vector<Double> pos0;
  for (int i = 0; i < nrant; ++i) {
    // Get antenna positions and convert to ITRF.
    MPosition mpos = MPosition::Convert(stationPositions[i], MPosition::ITRF)();
    if (i == 0) {
      pos0 = mpos.getValue().getVector();
    }
    Vector<Double> pos = mpos.getValue().getVector();
    MVPosition mvpos((pos[0] - pos0[0]), (pos[1] - pos0[1]),
                     (pos[2] - pos0[2]));
    itsAntMB.push_back(MBaseline(MVBaseline(mvpos), MBaseline::ITRF));
  }
  // Initialize the converters.
  // Set up the frame for epoch and antenna position.
  itsFrame.set(arrayPos, MDirection());
  // Create converter for phase center direction to J2000.
  itsDirToJ2000.set(phaseDir, MDirection::Ref(MDirection::J2000, itsFrame));
  // We can already try to convert the phase reference direction to J2000.
  // If it fails, it is a time-dependent direction (like SUN).
  itsMovingPhaseDir = false;
  try {
    itsPhaseDir = itsDirToJ2000();
    itsFrame.resetDirection(itsPhaseDir);
  } catch (...) {
    itsMovingPhaseDir = true;
  }
  itsFrame.set(MEpoch());
  // Create converter for MBaseline ITRF to J2000.
  itsBLToJ2000.set(MBaseline(), MBaseline::Ref(MBaseline::J2000, itsFrame));
  // Initialize the rest which is used to cache the UVW per antenna.
  // The cache is only useful if the MS is accessed in time order, but that
  // is normally the case.
  itsLastTime = 0;
  itsAntUvw.resize(nrant);
  itsUvwFilled.resize(nrant);
  itsUvwFilled = false;
}

std::array<double, 3> UVWCalculator::getUVW(unsigned int ant1,
                                            unsigned int ant2, double time) {
  // If a different time, we have to calculate the UVWs.
  if (time != itsLastTime) {
    itsLastTime = time;
    Quantum<Double> tm(time, "s");
    itsFrame.resetEpoch(MEpoch(MVEpoch(tm.get("d").getValue()), MEpoch::UTC));
    itsUvwFilled = false;
    // If phase dir is moving, calculate it for this time.
    if (itsMovingPhaseDir) {
      itsPhaseDir = itsDirToJ2000();
      itsFrame.resetDirection(itsPhaseDir);
    }
  }
  // Calculate the UVWs for this timestamp if not done yet.
  int ant = ant1;
  for (int i = 0; i < 2; ++i) {
    if (!itsUvwFilled[ant]) {
      MBaseline& mbl = itsAntMB[ant];
      mbl.getRefPtr()->set(itsFrame);  // attach frame
      MBaseline::Convert mcvt(mbl, MBaseline::J2000);
      MVBaseline bas = mcvt().getValue();
      MVuvw jvguvw(bas, itsPhaseDir.getValue());
      const casacore::Vector<double>& uvw =
          Muvw(jvguvw, Muvw::J2000).getValue().getVector();
      std::copy_n(uvw.data(), itsAntUvw[ant].size(), itsAntUvw[ant].data());
      itsUvwFilled[ant] = true;
    }
    ant = ant2;
  }
  // The UVW of the baseline is the difference of the antennae.
  return {itsAntUvw[ant2][0] - itsAntUvw[ant1][0],
          itsAntUvw[ant2][1] - itsAntUvw[ant1][1],
          itsAntUvw[ant2][2] - itsAntUvw[ant1][2]};
}

}  // namespace base
}  // namespace dp3
