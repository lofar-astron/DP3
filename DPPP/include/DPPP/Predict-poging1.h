//# Predict1.h: DPPP step class to Predict1 visibilities from a source model
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

#ifndef DPPP_Predict1_H
#define DPPP_Predict1_H

// @file
// @brief DPPP step class to Predict1 visibilities from a source model

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/Patch.h>
#include <DPPP/SourceDBUtil.h>
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

    // This class is a DPStep class to Predict1 visibilities and optionally beam

    typedef vector<Patch::ConstPtr> PatchList;
    typedef std::pair<size_t, size_t >Baseline;

    class Predict1: public DPStep
    {
    public:
      struct ThreadPrivateStorage
      {
        vector<double>    unknowns;
        vector<double>    uvw;
        vector<dcomplex>  model_patch; // Contains the model for only one patch
        vector<dcomplex>  model;
        vector<dcomplex>  model_subtr;
        vector<StationResponse::matrix22c_t> beamvalues; // [nst,nch]
        size_t            count_converged;

        //# Variables for conversion of directions to ITRF
        casa::MeasFrame                       measFrame;
        casa::MDirection::Convert             measConverter;
      };

      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      Predict1 (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~Predict1();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;


    private:
      void initThreadPrivateStorage(ThreadPrivateStorage &storage,
                                    size_t nDirection, size_t nStation,
                                    size_t nBaseline, size_t nChannel,
                                    size_t nChannelSubtr)
      {
        storage.unknowns.resize(nDirection * nStation * 8);
        storage.uvw.resize(nStation * 3);
        storage.model.resize(nBaseline * nChannel * 4);
        storage.model_patch.resize(nBaseline * nChannel * 4);
        storage.model_subtr.resize(nBaseline * nChannelSubtr * 4);
        storage.beamvalues.resize(nStation * nChannel);
        storage.count_converged = 0;

        // Create the Measure ITRF conversion info given the array position.
        // The time and direction are filled in later.
        info().arrayPosCopy();
        storage.measFrame.set (info().arrayPosCopy());
        storage.measFrame.set (casa::MEpoch(casa::MVEpoch(info().startTime()/86400), casa::MEpoch::UTC));
        storage.measConverter.set (casa::MDirection::J2000,
                              casa::MDirection::Ref(casa::MDirection::ITRF, storage.measFrame));
        // Do a dummy conversion, because Measure initialization does not
        // seem to be thread-safe.
        dir2Itrf(info().delayCenterCopy(),storage.measConverter);
      }

      // Convert a direction to ITRF.
      StationResponse::vector3r_t dir2Itrf (const casa::MDirection&, casa::MDirection::Convert&) const;

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      string           itsSourceDBName;
      bool             itsApplyBeam;
      bool             itsOneBeamPerPatch;
      bool             itsUseChannelFreq;
      Position         itsPhaseRef;

      uint             itsDebugLevel;

      vector<Baseline> itsBaselines;
      vector<ThreadPrivateStorage> itsThreadStorage;

      //# The info needed to calculate the station beams.
      vector<StationResponse::Station::Ptr> itsAntBeamInfo;

      PatchList        itsPatchList;

      string           itsOperation;

      NSTimer          itsTimer;
      NSTimer          itsTimerPredict1;
    };

  } //# end namespace
}

#endif
