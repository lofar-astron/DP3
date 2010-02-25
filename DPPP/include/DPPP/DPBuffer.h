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

    // @ingroup NDPPP

    // This class holds the data for one time slot in Array variables.
    // It makes heavy use of reference semantics to avoid data copying
    // when data are pushed from one step to another.
    // This means that a data array can be shared between DPStep objects. 
    // So if a DPStep object changes data in a buffer, it has to be sure
    // it can do it. If needed, Array::unique should be called to ensure
    // the array is not shared.
    //
    // The following data can be kept in a DPBuffer object.
    // <table>
    //  <tr>
    //   <td>TIME</td>
    //   <td>The time slot center of the current data (in MJD seconds).</td>
    //  </tr>
    //  <tr>
    //   <td>ROWNRS</td>
    //   <td>The row numbers of the current time slot. It can be empty
    //       when e.g. a time slot is inserted or if data are averaged.</td>
    //  </tr>
    //  <tr>
    //   <td>DATA</td>
    //   <td>The visibility data as [ncorr,nchan,nbaseline].</td>
    //  </tr>
    //  <tr>
    //   <td>FLAG</td>
    //   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
    //       Note that the ncorr axis is redundant because NDPPP will always
    //       have the same flag for all correlations. The reason all
    //       correlations are there is because the MS expects them.</td>
    //  </tr>
    //  <tr>
    //   <td>AMPLITUDE</td>
    //   <td>The amplitudes of the visibility data. It is used by the
    //       MedFlagger to avoid having to recalculate amplitudes.</td>
    //  </tr>
    //  <tr>
    //   <td>WEIGHT</td>
    //   <td>The data weights as [ncorr,nchan,nbaseline].
    //       Similarly to FLAG the ncorr axis is redundant because the
    //       same weight is used for all correlations.</td>
    //  </tr>
    //  <tr>
    //   <td>UVW</td>
    //   <td>The UVW coordinates in meters as [3,nbaseline].</td>
    //  </tr>
    //  <tr>
    //   <td>FULLRESFLAG</td>
    //   <td>The flags before any averaging was done. In the MS they are kept
    //       in column LOFAR_FULL_RES_FLAG. They are used to deal in BBS
    //       in a smart way with bandwidth and time smearing.
    //       The shape of the array is [nchan, ntimeavg, nbaseline], where
    //       ntimeavg is the number of time slots averaged to a single one.
    //       The number of channels averaged to a single one can be determined
    //       by dividing nchan by the number of channels in the data (or flags).
    //   </td>
    //  </tr>
    // </table>
    // The DATA and FLAG data members should always filled in, so each DPStep
    // should do that. Other data members do not need to be filled.
    // The DPInput::fetch functions should be used to get data for those
    // members. They take care that the buffer's data is used if available,
    // otherwise they get it from the DPInput object.
    // In that way as little memory as needed is used. Note that e.g. the
    // MedFlagger can use a lot of memory if a large time window is used.

    class DPBuffer
    {
    public:
      // Construct object with empty arrays.
      DPBuffer();

      // The copy constructor uses reference copies.
      DPBuffer (const DPBuffer&);

      // Assignment uses reference copies.
      DPBuffer& operator= (const DPBuffer&);

      // Set or get the visibility data per corr,chan,baseline.
      void setData (const casa::Cube<casa::Complex>& data)
        { itsData.reference (data); }
      const casa::Cube<casa::Complex>& getData() const
        { return itsData; }
      casa::Cube<casa::Complex>& getData()
        { return itsData; }

      // Set or get the flags per corr,chan,baseline.
      void setFlags (const casa::Cube<bool>& flags)
        { itsFlags.reference (flags); }
      const casa::Cube<bool>& getFlags() const
        { return itsFlags; }
      casa::Cube<bool>& getFlags()
        { return itsFlags; }

      // Set or get the amplitudes of the visibility data.
      // This is used by the MedFlagger to avoid calculating amplitudes
      // over and over again.
      void setAmplitudes (const casa::Cube<float>& ampl)
        { itsAmpl.reference (ampl); }
      const casa::Cube<float>& getAmplitudes() const
        { return itsAmpl; }
      casa::Cube<float>& getAmplitudes()
        { return itsAmpl; }

      // Set or get the weights per corr,chan,baseline.
      void setWeights (const casa::Cube<float>& weights)
        { itsWeights.reference (weights); }
      const casa::Cube<float>& getWeights() const
        { return itsWeights; }
      casa::Cube<float>& getWeights()
        { return itsWeights; }

      // Set or get the flags at the full resolution per chan,timeavg,baseline.
      void setFullResFlags (const casa::Cube<bool>& flags)
        { itsFullResFlags.reference (flags); }
      const casa::Cube<bool>& getFullResFlags() const
        { return itsFullResFlags; }
      casa::Cube<bool>& getFullResFlags()
        { return itsFullResFlags; }

      // Get or set the time.
      void setTime (double time)
        { itsTime = time; }
      double getTime() const
        { return itsTime; }

      // Get or set the row numbers used by the DPInput class.
      // It can be empty (e.g. when MSReader inserted an dummy time slot).
      void setRowNrs (const casa::Vector<uint>& rownrs)
        { itsRowNrs.reference (rownrs); }
      const casa::Vector<uint>& getRowNrs() const
        { return itsRowNrs; }

      // Get or set the UVW coordinates per baseline.
      void setUVW (const casa::Matrix<double>& uvw)
        { itsUVW.reference (uvw); }
      const casa::Matrix<double>& getUVW() const
        { return itsUVW; }
      casa::Matrix<double>& getUVW()
        { return itsUVW; }

      // Merge the flags into the pre-average flags.
      // For each flagged point, the corresponding pre-average flags are set.
      static void mergeFullResFlags (casa::Cube<bool>& fullResFlags,
                                     const casa::Cube<bool>& flags);

    private:
      double                    itsTime;
      casa::Vector<uint>        itsRowNrs;
      casa::Cube<casa::Complex> itsData;        //# ncorr,nchan,nbasel
      casa::Cube<float>         itsAmpl;        //# amplitude of data
      casa::Cube<bool>          itsFlags;       //# ncorr,nchan,nbasel
      casa::Matrix<double>      itsUVW;         //# 3,nbasel
      casa::Cube<float>         itsWeights;     //# nchan,nbasel
      casa::Cube<bool>          itsFullResFlags; //# fullres_nchan,ntimeavg,nbl
    };

  } //# end namespace
}

#endif
