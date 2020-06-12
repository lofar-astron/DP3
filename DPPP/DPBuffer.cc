// DPBuffer.cc: Buffer holding the data of a timeslot/band
// Copyright (C) 2010
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$
//
// @author Ger van Diepen

#include "DPBuffer.h"

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    DPBuffer::DPBuffer()
      : itsTime     (0),
        itsExposure (0)
    {}

    DPBuffer::DPBuffer (const DPBuffer& that)
    {
      operator= (that);
    }

    DPBuffer& DPBuffer::operator= (const DPBuffer& that)
    {
      if (this != &that) {
        itsTime     = that.itsTime;
        itsExposure = that.itsExposure;
        itsRowNrs.reference (that.itsRowNrs);
        itsData.reference (that.itsData);
        itsFlags.reference (that.itsFlags);
        itsWeights.reference (that.itsWeights);
        itsUVW.reference (that.itsUVW);
        itsFullResFlags.reference (that.itsFullResFlags);
      }
      return *this;
    }

    void DPBuffer::copy (const DPBuffer& that)
    {
      if (this != &that) {
        itsTime     = that.itsTime;
        itsExposure = that.itsExposure;
        itsRowNrs.assign (that.itsRowNrs);
        if (! that.itsData.empty()) {
          itsData.assign (that.itsData);
        }
        if (! that.itsFlags.empty()) {
          itsFlags.assign (that.itsFlags);
        }
        if (! that.itsWeights.empty()) {
          itsWeights.assign (that.itsWeights);
        }
        if (! that.itsUVW.empty()) {
          itsUVW.assign (that.itsUVW);
        }
        if (! that.itsFullResFlags.empty()) {
          itsFullResFlags.assign (that.itsFullResFlags);
        }
      }
    }

    void DPBuffer::referenceFilled (const DPBuffer& that)
    {
      if (this != &that) {
        itsTime     = that.itsTime;
        itsExposure = that.itsExposure;
        itsRowNrs.reference (that.itsRowNrs);
        if (! that.itsData.empty()) {
          itsData.reference (that.itsData);
        }
        if (! that.itsFlags.empty()) {
          itsFlags.reference (that.itsFlags);
        }
        if (! that.itsWeights.empty()) {
          itsWeights.reference (that.itsWeights);
        }
        if (! that.itsUVW.empty()) {
          itsUVW.reference (that.itsUVW);
        }
        if (! that.itsFullResFlags.empty()) {
          itsFullResFlags.reference (that.itsFullResFlags);
        }
      }
    }

    void DPBuffer::mergeFullResFlags (Cube<bool>& fullResFlags,
                                      const Cube<bool>& flags)
    {
      // Flag shape is [ncorr, newnchan, nbl].
      // FullRes shape is [orignchan, navgtime, nbl]
      // where orignchan = navgchan * newnchan.
      const IPosition& fullResShape = fullResFlags.shape();
      const IPosition& flagShape    = flags.shape();
      int orignchan = fullResShape[0];
      int newnchan  = flagShape[1];
      int navgchan  = orignchan / newnchan;
      int navgtime  = fullResShape[1];
      int nbl       = fullResShape[2];
      int ncorr     = flagShape[0];
      bool* fullResPtr = fullResFlags.data();
      const bool* flagPtr = flags.data();
      // Loop over all baselines and new channels.
      // Only use the first correlation in the loop.
      for (int j=0; j<nbl; ++j) {
        for (int i=0; i<newnchan; ++i) {
          // If ta data point is flagged, the flags in the corresponding
          // FullRes window have to be set.
          // This is needed in case a data point is further averaged.
          if (*flagPtr) {
            for (int i=0; i<navgtime; ++i) {
              std::fill (fullResPtr, fullResPtr+navgchan, true);
              fullResPtr += orignchan;
            }
            fullResPtr -= navgtime*orignchan;
          }
          flagPtr   += ncorr;
          fullResPtr += navgchan;
        }
        // Set pointer to next baseline.
        fullResPtr += (navgtime-1)*orignchan;
      }
    }

  } // end namespace
}
