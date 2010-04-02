//# UVWCalculator.cc: Class to calculate UVW coordinates
//# Copyright (C) 2010
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$
//#
//# @author Ger van Diepen

#include <DPPP/UVWCalculator.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MBaseline.h>
#include <measures/Measures/MCBaseline.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/Muvw.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    UVWCalculator::UVWCalculator()
    {}

    UVWCalculator::UVWCalculator (const MDirection& phaseDir,
                                  const MPosition& arrayPos,
                                  const vector<MPosition>& stationPositions)
    {
      // Convert the station positions to a baseline in ITRF.
      int nrant = stationPositions.size();
      Vector<Double> pos0;
      for (int i=0; i<nrant; ++i) {
        // Get antenna positions and convert to ITRF.
        MPosition mpos = MPosition::Convert (stationPositions[i],
                                             MPosition::ITRF)();
        if (i == 0) {
          pos0 = mpos.getValue().getVector();
        }
        Vector<Double> pos = mpos.getValue().getVector();
        MVPosition mvpos((pos[0] - pos0[0]),
                         (pos[1] - pos0[1]),
                         (pos[2] - pos0[2]));
        itsAntMB.push_back (MBaseline (MVBaseline(mvpos), MBaseline::ITRF));
      }
      // Try to use the PHASE_DIR in J2000.
      // It fails if a moving source (e.g. SUN) is used.
      itsMovingPhaseDir = false;
      itsOrigPhaseDir   = phaseDir;
      try {
        itsPhaseDir = MDirection::Convert (phaseDir, MDirection::J2000)();
      } catch (...) {
        itsMovingPhaseDir = true;
      }
      // Create a reference frame. Use the middle antenna as array position.
      itsFrame.set (arrayPos);
      itsFrame.set (itsPhaseDir);
      // Initialize the rest which is used to cache the UVW per antenna.
      // The cache is only useful if the MS is accessed in time order, but that
      // is normally the case.
      itsLastTime = 0;
      itsAntUvw.resize (nrant);
      itsUvwFilled.resize (nrant);
      itsUvwFilled = false;
    }

    Vector<double> UVWCalculator::getUVW (uint ant1, uint ant2, double time)
    {
      // If a different time, we have to calculate the UVWs.
      if (time != itsLastTime) {
        itsLastTime = time;
        Quantum<Double> tm(time, "s");
        itsFrame.set (MEpoch(MVEpoch(tm.get("d").getValue()), MEpoch::UTC));
        itsUvwFilled = false;
        // If phase dir is moving, calculate it for this time.
        if (itsMovingPhaseDir) {
          itsPhaseDir = MDirection::Convert (itsOrigPhaseDir,
                                             MDirection::J2000)();
        }
      }
      // Calculate the UVWs for this timestamp if not done yet.
      int ant = ant1;
      for (int i=0; i<2; ++i) {
        if (!itsUvwFilled[ant]) {
          MBaseline& mbl = itsAntMB[ant];
          mbl.getRefPtr()->set(itsFrame);       // attach frame
          MBaseline::Convert mcvt(mbl, MBaseline::J2000);
          MVBaseline bas = mcvt().getValue();
          MVuvw jvguvw(bas, itsPhaseDir.getValue());
          itsAntUvw[ant] = Muvw(jvguvw, Muvw::J2000).getValue().getVector();
          itsUvwFilled[ant] = true;
        }
        ant = ant2;
      }
      // The UVW of the baseline is the difference of the antennae.
      return itsAntUvw[ant2] - itsAntUvw[ant1];
    }

  } //# end namespace
}
