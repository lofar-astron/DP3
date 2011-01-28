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
#include <casa/Arrays/ArrayMath.h>
#include <iostream>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    PhaseShift::PhaseShift (DPInput* input,
                            const ParSet& parset, const string& prefix)
      : itsInput  (input),
        itsName   (prefix),
        itsRaStr  (parset.getString(prefix+"ra", "")),
        itsDecStr (parset.getString(prefix+"dec", ""))
    {}

    PhaseShift::~PhaseShift()
    {}

    void PhaseShift::updateInfo (DPInfo& info)
    {
      // Default phase center is the original one.
      MDirection newDir(info.originalPhaseCenter());
      if (itsRaStr.empty() && itsDecStr.empty()) {
        itsRa  = MVTime::read (itsRaStr);
        itsDec = MVTime::read (itsDecStr);
        newDir = MDirection(itsRa, itsDec, MDirection::J2000);
      }
      new casa::UVWMachine(info.phaseCentre(), newDir, false, true);
      info.setNewPhaseCenter (newDir);
    }

    void PhaseShift::show (std::ostream& os) const
    {
      os << "PhaseShift " << itsName << std::endl;
      os << "  ra:             " << itsRaStr << std::endl;
      os << "  dec:            " << itsDecStr << std::endl;
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
      RefRows rowNrs(buf.getRowNrs());
      itsBuf.getUVW()  = itsInput->fetchUVW (buf, rowNrs);
      itsBuf.getData().resize (buf.getData().shape());
      const Complex* indata = buf.getData().data();
      Complex* outdata = itsBuf.getData().data();
      uint ncorr = buf.getData().shape()[0];
      uint nchan = buf.getData().shape()[1];
      uint ntime = buf.getData().shape()[2];
      VectorIterator uvwIter(itsBuf.getUVW());
      for (uint i=0; i<ntime; ++i) {
        // Convert the uvw coordinates and get the phase shift term.
        itsMachine.convertUVW (phase, uvwIter.vector());
        Complex phasor(cos(phase), sin(phase));
        // Shift the phase of the data.
        for (uint j=0; j<nchan; ++j) {
          for (uint k=0; j<ncorr; ++k) {
            *outdata = *indata * phasor;
          }
        }
        uvwIter.next();
     }
     itsTimer.stop();
     getNextStep()->process (buf);
     return true;
    }

    void PhaseShift::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
