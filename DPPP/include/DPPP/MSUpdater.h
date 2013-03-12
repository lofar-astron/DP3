//# MSUpdater.h: DPPP step writing to an MS
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

#ifndef DPPP_MSUPDATER_H
#define DPPP_MSUPDATER_H

// @file
// @brief DPPP step writing to an MS

#include <DPPP/DPStep.h>
#include <Common/LofarTypes.h>
#include <tables/Tables/RefRows.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    //# Forward Declarations.
    class MSReader;

    // @ingroup NDPPP

    // This class updates the flags in an existing MeasurementSet.
    // Hardly anything is done in this class.
    // It uses function putFlags in MSReader to do the actual write.
    //
    // Like MSWriter it adds an entry to the HISTORY table of the MS
    // containing the parset values and DPPP version.

    class MSUpdater: public DPStep
    {
    public:
      MSUpdater (MSReader*, const ParameterSet& parset,
                 const std::string& prefix, int needWrite);

      virtual ~MSUpdater();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;

    private:
      // Write the flags at the given row numbers.
      void putFlags (const casa::RefRows& rowNrs,
                     const casa::Cube<bool>& flags);

      // Write the data at the given row numbers.
      void putData (const casa::RefRows& rowNrs,
                    const casa::Cube<casa::Complex>& data);

      //# Data members
      MSReader*   itsReader;
      bool        itsWriteData;
      uint        itsNrCorr;
      uint        itsNrChan;
      uint        itsNrBl;
      uint        itsNrTimesFlush; //# flush every N time slots (0=no flush)
      uint        itsNrDone;       //# nr of time slots written
      NSTimer     itsTimer;
    };

  } //# end namespace
}

#endif
