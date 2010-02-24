//# DPInput.h: Abstract base class for a DPStep generating input
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

#ifndef DPPP_DPINPUT_H
#define DPPP_DPINPUT_H

// @file
// @brief Abstract base class for a DPStep generating input

#include <DPPP/DPStep.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/UVWCalculator.h>
#include <DPPP/FlagCounter.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/RefRows.h>
#include <casa/Arrays/Slicer.h>
#include <Common/lofar_vector.h>

namespace LOFAR {
  namespace DPPP {

    // @ingroup NDPPP

    // This class is the abstract base class for a DPStep object that
    // handles the input. A concrete example is MSReader that reads the
    // data from a MeasurementSet. However, it is also possible to have
    // input steps generating data on the fly as done in test programs
    // like tAverager.cc.
    //
    // A particular task of the class is to fetch the input for various
    // data items like weight, uvw, etc.. This is done by testing if the
    // item's data array is in the DPBuffer. If so, it will be returned.
    // Otherwise the appropriate 'get' function will be called to read the
    // data array from the input.
    // The derived classes should implement those 'get' functions, unless
    // they are sure the data arrays are always put in the buffer.

    class DPInput: public DPStep
    {
    public:
      virtual ~DPInput();

      // Read the UVW at the given row numbers.
      // The default implementation throws an exception.
      virtual casa::Matrix<double> getUVW (const casa::RefRows& rowNrs);

      // Read the weights at the given row numbers.
      // The default implementation throws an exception.
      virtual casa::Cube<float> getWeights (const casa::RefRows& rowNrs);

      // Read the preAvg flags (LOFAR_PREAVG_FLAG) at the given row numbers.
      // The default implementation throws an exception.
      virtual casa::Cube<bool> getPreAvgFlags (const casa::RefRows& rowNrs);

      // Read the given data column at the given row numbers.
      // The default implementation throws an exception.
      virtual casa::Cube<casa::Complex> getData
      (const casa::String& columnName, const casa::RefRows& rowNrs);

      // Get info.
      uint ncorr() const
        { return itsNrCorr; }
      uint nchan() const
        { return itsNrChan; }
      uint nbaselines() const
        { return itsNrBl; }

      // Fetch the PreAvg flags.
      // If defined in the buffer, they are taken from there.
      // Otherwise there are read from the input.
      // If not defined in the input, they are filled using the flags in the
      // buffer assuming that no averaging has been done so far.
      // If defined, they can be merged with the buffer's flags which means
      // that if an averaged channel is flagged, the corresponding PreAvg
      // flags are set.
      casa::Cube<bool> fetchPreAvgFlags (const DPBuffer& buf,
                                         const casa::RefRows& rowNrs,
                                         bool merge=false);

      // Fetch the weights.
      // If defined in the buffer, they are taken from there.
      // Otherwise there are read from the input.
      casa::Cube<float> fetchWeights (const DPBuffer& buf,
                                      const casa::RefRows& rowNrs);

      // Fetch the UVW.
      // If defined in the buffer, they are taken from there.
      // Otherwise there are read from the input.
      casa::Matrix<double> fetchUVW (const DPBuffer& buf,
                                     const casa::RefRows& rowNrs);

    protected:
      uint itsNrChan;
      uint itsNrCorr;
      uint itsNrBl;
    };

  } //# end namespace
}

#endif
