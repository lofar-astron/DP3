//# DPBuffer.h: Buffer holding the data of a timeslot/band
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

#ifndef DPPP_DPBUFFER_H
#define DPPP_DPBUFFER_H

/// @file
/// @brief Buffer holding the data of a timeslot/band

#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Cube.h>
#include <casa/BasicSL/Complex.h>

namespace LOFAR {
  namespace DPPP {

    // @ingroup DPPP

    class DPBuffer
    {
    public:
      // Construct object with empty arrays.
      DPBuffer();

      // The copy constructor uses reference copies.
      DPBuffer (const DPBuffer&);

      // Assignment uses reference copies.
      DPBuffer& operator= (const DPBuffer&);

      void setData (const casa::Cube<casa::Complex>& data)
        { itsData.reference (data); }
      const casa::Cube<casa::Complex>& getData() const
        { return itsData; }
      casa::Cube<casa::Complex>& getData()
        { return itsData; }

      void setFlags (const casa::Cube<bool>& flags)
        { itsFlags.reference (flags); }
      const casa::Cube<bool>& getFlags() const
        { return itsFlags; }
      casa::Cube<bool>& getFlags()
        { return itsFlags; }

      void setAmplitudes (const casa::Cube<float>& ampl)
        { itsAmpl.reference (ampl); }
      const casa::Cube<float>& getAmplitudes() const
        { return itsAmpl; }
      casa::Cube<float>& getAmplitudes()
        { return itsAmpl; }

      void setWeights (const casa::Cube<float>& weights)
        { itsWeights.reference (weights); }
      const casa::Cube<float>& getWeights() const
        { return itsWeights; }
      casa::Cube<float>& getWeights()
        { return itsWeights; }

      void setPreAvgFlags (const casa::Cube<bool>& flags)
        { itsPreAvgFlags.reference (flags); }
      const casa::Cube<bool>& getPreAvgFlags() const
        { return itsPreAvgFlags; }
      casa::Cube<bool>& getPreAvgFlags()
        { return itsPreAvgFlags; }

      void setTime (double time)
        { itsTime = time; }
      double getTime() const
        { return itsTime; }

      void setRowNrs (const casa::Vector<uint>& rownrs)
        { itsRowNrs.reference (rownrs); }
      const casa::Vector<uint>& getRowNrs() const
        { return itsRowNrs; }

      void setUVW (const casa::Matrix<double>& uvw)
        { itsUVW.reference (uvw); }
      const casa::Matrix<double>& getUVW() const
        { return itsUVW; }

      bool hasNoFlags() const
        { return itsFlags.empty(); }

      // Resize the data and flag arrays to the other buffer sizes.
      void resize (const DPBuffer& other)
      {
        itsData.resize (other.getData().shape());
        itsFlags.resize (other.getFlags().shape());
      }

      // Take care that the arrays are unique (not referenced elsewhere).
      void makeUnique()
      {
        itsData.unique();
        itsFlags.unique();
      }

      // Merge the flags into the pre-average flags.
      // For each flagged point, the corresponding pre-average flags are set.
      static void mergePreAvgFlags (casa::Cube<bool>& preAvgFlags,
                                    const casa::Cube<bool>& flags);

    private:
      double                    itsTime;
      casa::Vector<uint>        itsRowNrs;
      casa::Cube<casa::Complex> itsData;        //# ncorr,nchan,nbasel
      casa::Cube<float>         itsAmpl;        //# amplitude of data
      casa::Cube<bool>          itsFlags;       //# ncorr,nchan,nbasel
      casa::Matrix<double>      itsUVW;         //# 3,nbasel
      casa::Cube<float>         itsWeights;     //# nchan,nbasel
      casa::Cube<bool>          itsPreAvgFlags; //# preavg_nchan,ntimeavg,nbasel
    };

  } //# end namespace
}

#endif
