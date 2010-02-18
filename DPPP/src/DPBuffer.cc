//# DPBuffer.cc: Buffer holding the data of a timeslot/band
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

#include <DPPP/DPBuffer.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    DPBuffer::DPBuffer()
    {}

    DPBuffer::DPBuffer (const DPBuffer& that)
    {
      operator= (that);
    }

    DPBuffer& DPBuffer::operator= (const DPBuffer& that)
    {
      if (this != &that) {
        itsTime = that.itsTime;
        itsRowNrs.reference (that.itsRowNrs);
        itsData.reference (that.itsData);
        itsAmpl.reference (that.itsAmpl);
        itsFlags.reference (that.itsFlags);
        itsWeights.reference (that.itsWeights);
        itsPreAvgFlags.reference (that.itsPreAvgFlags);
      }
      return *this;
    }

    void DPBuffer::mergePreAvgFlags (Cube<bool>& preAvgFlags,
                                     const Cube<bool>& flags)
    {
      // Flag shape is [ncorr, newnchan, nbl].
      // PreAvg shape is [orignchan, navgtime, nbl]
      // where orignchan = navgchan * newnchan.
      const IPosition& preAvgShape = preAvgFlags.shape();
      const IPosition& flagShape   = flags.shape();
      int orignchan = preAvgShape[0];
      int newnchan  = flagShape[1];
      int navgchan  = orignchan / newnchan;
      int navgtime  = preAvgShape[1];
      int nbl       = preAvgShape[2];
      int ncorr     = flagShape[0];
      bool* preAvgPtr = preAvgFlags.data();
      const bool* flagPtr = flags.data();
      // Loop over all baselines and new channels.
      // Only use the first correlation in the loop.
      for (int j=0; j<nbl; ++j) {
        for (int i=0; i<newnchan; ++i) {
          // If ta data point is flagged, the flags in the corresponding
          // PreAvg window have to be set.
          // This is needed in case a data point is further averaged.
          if (*flagPtr) {
            for (int i=0; i<navgtime; ++i) {
              std::fill (preAvgPtr, preAvgPtr+navgchan, true);
              preAvgPtr += orignchan;
            }
            preAvgPtr -= navgtime*orignchan;
          }
          flagPtr   += ncorr;
          preAvgPtr += navgchan;
        }
        // Set pointer to next baseline.
        preAvgPtr += (navgtime-1)*orignchan;
      }
    }

  } //# end namespace
}
