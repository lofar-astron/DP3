//# PhaseShift.cc: DPPP step class to average in time and/or freq
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

#include <lofar_config.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/BasicSL/Constants.h>
#include <iostream>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    PhaseShift::PhaseShift (DPInput* input,
                            const ParSet& parset, const string& prefix)
      : itsInput   (input),
        itsName    (prefix),
        itsCenter  (parset.getStringVector(prefix+"phasecenter")),
        itsMachine (0)
    {}

    PhaseShift::~PhaseShift()
    {
      delete itsMachine;
    }

    void PhaseShift::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
      // Default phase center is the original one.
      MDirection newDir(itsInput->phaseCenter());
      ////      bool original = true;
      bool original = false;
      if (! itsCenter.empty()) {
        newDir = handleCenter();
        original = false;
      }
      itsMachine = new casa::UVWMachine(info.phaseCenter(), newDir);
      info.setPhaseCenter (newDir, original);
      // Calculate 2*pi*freq/C to get correct phase term (in wavelengths).
      const Vector<double>& freq = itsInput->chanFreqs(info.nchanAvg());
      itsFreqC.reserve (freq.size());
      for (uint i=0; i<freq.size(); ++i) {
        itsFreqC.push_back (2. * C::pi * freq[i] / C::c);
      }
    }

    void PhaseShift::show (std::ostream& os) const
    {
      os << "PhaseShift " << itsName << std::endl;
      os << "  phasecenter:    " << itsCenter << std::endl;
    }

    void PhaseShift::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " PhaseShift " << itsName << endl;
    }

    bool PhaseShift::process (const DPBuffer& buf)
    {
      itsTimer.start();
      DPBuffer newBuf(buf);
      RefRows rowNrs(newBuf.getRowNrs());
      if (newBuf.getUVW().empty()) {
        newBuf.getUVW().reference (itsInput->fetchUVW (newBuf, rowNrs));
      }
      newBuf.getData().unique();
      newBuf.getUVW().unique();
//       // Get the uvw-s per station, so we can do convertUVW per station
//       // instead of per baseline.
//       const Vector<Int>& ant1 = itsInput->getAnt1();
//       const Vector<Int>& ant2 = itsInput->getAnt2();
//       vector<Vector<double> > antUVW;
//       antUVW.resize (1 + std::max (max(ant1), max(ant2)));
//       // The first station is taken as reference (has value 0).
//       Vector<double> nullUVW(3, 0.);
//       uint ndone = 0;
//       bool moreTodo = true;
//       bool setFirst = true;      // set first unknown uvw to 0?
//       while (moreTodo) {
//         uint oldNdone = ndone;
//         moreTodo = false;
//         VectorIterator<double> iterUVW (newBuf.getUVW(), 0);
//         for (uint i=0; i<ant1.size(); ++i) {
//           int a1 = ant1[i];
//           int a2 = ant2[i];
//           cout<<a1<<' '<<a2<<' '<<iterUVW.vector()<<endl;
//           if (antUVW[a1].empty()) {
//             if (antUVW[a2].empty()) {
//               if (setFirst) {
//                 antUVW[a1] = nullUVW;
//                 antUVW[a2] = iterUVW.vector();
//                 setFirst = false;
//                 ndone++;
//               } else {
//                 // Both stations are still unknown, so need a next round.
//                 moreTodo = true;
//               }
//             } else {
//               // Get uvw of first station.
//               antUVW[a1] = antUVW[a2] - iterUVW.vector();
//               ndone++;
//             }
//           } else if (antUVW[a2].empty()) {
//             antUVW[a2] = antUVW[a1] + iterUVW.vector();
//             ndone++;
//           } else {
//             // UVW should be about the same.
//             ASSERT (allNear(antUVW[a2] - antUVW[a1], iterUVW.vector(), 1e-7));
//           }
//           cout << "  "<<antUVW[a1]<<antUVW[a2]<<endl;
//           iterUVW.next();
//         }
//         if (ndone == oldNdone) {
//           setFirst = moreTodo;
//           moreTodo = true;
//         }
//       }
//       vector<double> phases(nant);
//       for (uInt i=0; i<nant; ++i) {
//         if (! antUVW[i].empty()) {
//           itsMachine->convertUVW (phases[i], antUVW[i]);
//         }
//       }
      Complex* data = newBuf.getData().data();
      uint ncorr = newBuf.getData().shape()[0];
      uint nchan = newBuf.getData().shape()[1];
      uint nbl   = newBuf.getData().shape()[2];
      VectorIterator<double> uvwIter(newBuf.getUVW(), 0);
      double phase;
      //# If ever in the future a time dependent phase center is used,
      //# the machine must be reset for each new time, thus each new call
      //# to process.
      for (uint i=0; i<nbl; ++i) {
        // Convert the uvw coordinates and get the phase shift term.
        itsMachine->convertUVW (phase, uvwIter.vector());
        for (uint j=0; j<nchan; ++j) {
          // Shift the phase of the data of this baseline.
          // Convert the phase term to wavelengths (and apply 2*pi).
          // u_wvl = u_m / wvl = u_m * freq / c
          double phasewvl = phase * itsFreqC[j];
          Complex phasor(cos(phasewvl), sin(phasewvl));
          for (uint k=0; k<ncorr; ++k) {
            *data++ *= phasor;
          }
        }
        uvwIter.next();
     }
     itsTimer.stop();
     getNextStep()->process (newBuf);
     return true;
    }

    void PhaseShift::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

    MDirection PhaseShift::handleCenter()
    {
      // The phase center must be given in J2000 as two values (ra,dec).
      // In the future time dependent frames can be done as in UVWFlagger.
      ASSERTSTR (itsCenter.size() == 2,
                 "2 values must be given in PhaseShift phasecenter");
      ///ASSERTSTR (itsCenter.size() < 4,
      ///"Up to 3 values can be given in UVWFlagger phasecenter");
      MDirection phaseCenter;
      if (itsCenter.size() == 1) {
        string str = toUpper(itsCenter[0]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid source type"
                   " in UVWFlagger phasecenter");
        return MDirection(tp);
      }
      Quantity q0, q1;
      ASSERTSTR (MVAngle::read (q0, itsCenter[0]),
                 itsCenter[0] << " is an invalid RA or longitude"
                 " in UVWFlagger phasecenter");
      ASSERTSTR (MVAngle::read (q1, itsCenter[1]),
                 itsCenter[1] << " is an invalid DEC or latitude"
                 " in UVWFlagger phasecenter");
      MDirection::Types type = MDirection::J2000;
      if (itsCenter.size() > 2) {
        string str = toUpper(itsCenter[2]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid direction type in UVWFlagger"
                   " in UVWFlagger phasecenter");
      }
      return MDirection(q0, q1, type);
    }

  } //# end namespace
}
