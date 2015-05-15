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

#ifndef DPPP_ApplyBeam_H
#define DPPP_ApplyBeam_H

// @file
// @brief DPPP step class to ApplyBeam visibilities from a source model

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/Patch.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/ModelComponent.h>
#include <StationResponse/Station.h>
#include <StationResponse/Types.h>
#include <casa/Arrays/Cube.h>
#include <casa/Quanta/MVEpoch.h>
#include <measures/Measures/MEpoch.h>
#include <casa/Arrays/ArrayMath.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to ApplyBeam visibilities with optionally beam

    typedef std::pair<size_t, size_t> Baseline;
    typedef std::pair<ModelComponent::ConstPtr, Patch::ConstPtr> Source;

    class ApplyBeam: public DPStep
    {
      public:
        // Construct the object.
        // Parameters are obtained from the parset using the given prefix.
        ApplyBeam (DPInput*, const ParameterSet&, const string& prefix);

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
            const DPInfo& info, double time, T* data0,
            const StationResponse::vector3r_t& srcdir,
            const StationResponse::vector3r_t& refdir,
            const StationResponse::vector3r_t& tiledir,
            const vector<StationResponse::Station::Ptr>& antBeamInfo,
            vector<StationResponse::matrix22c_t>& beamValues,
            bool useChannelFreq, bool invert);

      private:
        StationResponse::vector3r_t dir2Itrf(
            const casa::MDirection& dir,
            casa::MDirection::Convert& measConverter);

        //# Data members.
        DPInput*             itsInput;
        string               itsName;
        DPBuffer             itsBuffer;
        bool                 itsInvert;
        bool                 itsUseChannelFreq;
        Position             itsPhaseRef;

        uint                 itsDebugLevel;

        vector<Baseline>     itsBaselines;

        // Vector containing info on converting baseline uvw to station uvw
        vector<int>          itsUVWSplitIndex;

        // UVW coordinates per station (3 coordinates per station)
        casa::Matrix<double> itsUVW;

        // The info needed to calculate the station beams.
        vector<vector<StationResponse::Station::Ptr> > itsAntBeamInfo;
        vector<casa::MeasFrame> itsMeasFrames;
        vector<casa::MDirection::Convert> itsMeasConverters;
        vector<vector<StationResponse::matrix22c_t> > itsBeamValues;

        NSTimer itsTimer;
    };

  } //# end namespace
}
#endif

#include "ApplyBeam.tcc"
