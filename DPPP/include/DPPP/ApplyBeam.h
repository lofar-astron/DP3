//# ApplyBeam.h: DPPP step class to ApplyBeam visibilities from a source model
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
//# $Id:
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_APPLYBEAM_H
#define DPPP_APPLYBEAM_H

// @file
// @brief DPPP step class to apply the beam model (optionally inverted)

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/Position.h>
#include <StationResponse/Station.h>
#include <StationResponse/Types.h>
#include <casacore/casa/Arrays/Cube.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to apply the beam model, optionally inverted.
    // The input MeasurementSet it operates on, must have the LOFAR subtables
    // defining the station layout and tiles/dipoles used.

    class ApplyBeam: public DPStep
    {
      public:
        // Modes for the beam: only array factor, only element beam, or both (default)
        enum BeamMode {DEFAULT=1, ARRAY_FACTOR=2, ELEMENT=3};

        // Construct the object.
        // Parameters are obtained from the parset using the given prefix.
        ApplyBeam (DPInput*, const ParameterSet&, const string& prefix, bool substep=false);

        ApplyBeam();

        virtual ~ApplyBeam();

        // Process the data.
        // It keeps the data.
        // When processed, it invokes the process function of the next step.
        virtual bool process(const DPBuffer&);

        // Finish the processing of this step and subsequent steps.
        virtual void finish();

        // Update the general info.
        virtual void updateInfo(const DPInfo&);

        // Show the step parameters.
        virtual void show(std::ostream&) const;

        // Show the timings.
        virtual void showTimings(std::ostream&, double duration) const;

        bool invert() {
          return itsInvert;
        }

        template<typename T>
        static void applyBeam(
            const DPInfo& info, double time, T* data0, float* weight0,
            const StationResponse::vector3r_t& srcdir,
            const StationResponse::vector3r_t& refdir,
            const StationResponse::vector3r_t& tiledir,
            const vector<StationResponse::Station::Ptr>& antBeamInfo,
            vector<StationResponse::matrix22c_t>& beamValues,
            bool useChannelFreq, bool invert, int mode,
            bool doUpdateWeights=false);

      private:
        StationResponse::vector3r_t dir2Itrf(
            const casacore::MDirection& dir,
            casacore::MDirection::Convert& measConverter);

        //# Data members.
        DPInput*             itsInput;
        string               itsName;
        DPBuffer             itsBuffer;
        bool                 itsInvert;
        bool                 itsUpdateWeights;
        bool                 itsUseChannelFreq;
        Position             itsPhaseRef;
        BeamMode             itsMode;

        uint                 itsDebugLevel;

        // The info needed to calculate the station beams.
        vector<vector<StationResponse::Station::Ptr> > itsAntBeamInfo;
        vector<casacore::MeasFrame> itsMeasFrames;
        vector<casacore::MDirection::Convert> itsMeasConverters;
        vector<vector<StationResponse::matrix22c_t> > itsBeamValues;

        NSTimer itsTimer;
    };

  } //# end namespace
}
#endif

#include "ApplyBeam.tcc"
