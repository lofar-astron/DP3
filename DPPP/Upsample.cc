//# Upsample.cc: DPPP step class to Upsample visibilities
//# Copyright (C) 2013
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include "Upsample.h"

#include <iostream>

#include "../Common/ParameterSet.h"

#include <casacore/casa/BasicMath/Math.h> // nearAbs
#include <casacore/casa/Arrays/ArrayLogical.h> // anyTrue

#include <iomanip>
#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    Upsample::Upsample (DPInput*, 
                      const ParameterSet& parset,
                      const string& prefix)
    : itsOldTimeInterval(0),
      itsTimeStep(parset.getInt(prefix + "timestep")),
      itsFirstToFlush(0)
    {
      itsBuffers.resize(itsTimeStep);
    }

    Upsample::~Upsample()
    {}

    void Upsample::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();

      itsOldTimeInterval = info().timeInterval();
      info().setTimeInterval(itsOldTimeInterval / itsTimeStep);

      info().setMetaChanged();
    }

    void Upsample::show (std::ostream& os) const
    {
      os << "Upsample " << itsName << endl;
      os << "  time step : " << itsTimeStep <<endl;
    }

    bool Upsample::process (const DPBuffer& bufin)
    {

      double time0 = bufin.getTime() - 0.5 * itsOldTimeInterval;
      double exposure = bufin.getExposure() / itsTimeStep;

      // Duplicate the input buffer itsTimeStep times
      for (uint i=0; i<itsTimeStep; ++i) {
        itsBuffers[i].copy (bufin);
        // Update the time centroid and time exposure
        itsBuffers[i].setTime(time0 + info().timeInterval() * (i+0.5));
        itsBuffers[i].setExposure(exposure);
      }

      if (itsPrevBuffers.empty()) {
        // First time slot, ask for next time slot first
        itsPrevBuffers.resize(itsTimeStep);
        for (uint i=0; i<itsTimeStep; ++i) {
          itsPrevBuffers[i].copy(itsBuffers[i]); // No shallow copy
        }
        return false;
      }

      // Flush the itsPrevBuffers. Skip parts at the beginning as determined
      // in the previous call to process. Skip parts at the end as determined
      // now. Also, determine at which step the next buffer should start to
      // flush in the next call of process.
      uint curIndex = 0; // Index in the current buffers
      for (uint prevIndex=itsFirstToFlush; prevIndex<itsTimeStep; prevIndex++) {
        itsFirstToFlush = 0; // reset for next use
        // Advance curIndex until
        // buffers[curIndex].time >= prevBuffers[prevIndex].time
        while (itsPrevBuffers[prevIndex].getTime() > itsBuffers[curIndex].getTime() &&
               !nearAbs(itsPrevBuffers[prevIndex].getTime(),
                        itsBuffers[curIndex].getTime(),
                        0.4*info().timeInterval()) ) {
          curIndex++;
        }
        if (nearAbs(itsPrevBuffers[prevIndex].getTime(),
                    itsBuffers[curIndex].getTime(),
                    0.4*info().timeInterval())) {
           // Found double buffer, choose which one to use
           // If both totally flagged, prefer prevbuffer
           if (allTrue(itsBuffers[curIndex].getFlags())) {
             // Use prevBuffer
             itsFirstToFlush = curIndex+1;
             getNextStep()->process(itsPrevBuffers[prevIndex]);
           } else {
             // Use next buffer; stop processing prevbuffer.
             // This will uncorrectly give flagged if data has been flagged
             // and a time slot has been inserted and itsTimeStep > 2.
             break;
           }
        } else {
          // No double buffer, just flush the prevbuffer
          getNextStep()->process(itsPrevBuffers[prevIndex]);
        }
      }

      itsPrevBuffers.swap(itsBuffers); // itsBuffers will be overwritten later

      return false;
    }


    void Upsample::finish()
    {
      // Flush itsPrevBuffers
      for (uint i=itsFirstToFlush; i<itsTimeStep; ++i) {
        getNextStep()->process(itsPrevBuffers[i]);
      }

      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
