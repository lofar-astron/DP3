//# GainCal.h: DPPP step class to calibrate (direction independent) gains
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
//# $Id: GainCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_GAINCAL_H
#define DPPP_GAINCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/Patch.h>
#include <ParmDB/ParmFacade.h>
#include <ParmDB/ParmSet.h>
#include <DPPP/SourceDBUtil.h>
#include <StationResponse/Station.h>
#include <StationResponse/Types.h>
#include <ParmDB/Parm.h>
#include <casa/Arrays/Cube.h>
#include <casa/Quanta/MVEpoch.h>
#include <measures/Measures/MEpoch.h>
#include <casa/Arrays/ArrayMath.h>

namespace LOFAR {

namespace {

}

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class to calibrate (direction independent) gains.

    typedef vector<Patch::ConstPtr> PatchList;
    typedef std::pair<size_t, size_t >Baseline;

    class GainCal: public DPStep
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
      GainCal (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~GainCal();

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
      struct StefVecs {
        casa::Matrix<casa::DComplex> g;
        casa::Matrix<casa::DComplex> gold;
        casa::Matrix<casa::DComplex> gx;
        casa::Matrix<casa::DComplex> gxx;
        casa::Matrix<casa::DComplex> h;
        casa::Matrix<casa::DComplex> z;
      };

      void exportToMatlab(uint ch);

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

      // Perform stefcal (polarized or unpolarized)
      void stefcal(string mode, uint solint=1);

      // Counts the number of antennas with non-flagged data, adds this to
      // dataPerAntenna
      void countAntUsedNotFlagged (const casa::Bool* flag);

      // Set a map for the used antennas
      void setAntennaMaps ();

      // Remove rows and colums corresponding to antennas with too much
      // flagged data from vis and mvis
      void removeDeadAntennas ();

      // Fills the matrices itsVis and itsMVis (single precision model)
      void fillMatrices (casa::Complex* model, casa::Complex* data, float* weight,
                         const casa::Bool* flag);

      // Fills the matrices itsVis and itsMVis (for double precision model)
      void fillMatrices (dcomplex* model, casa::Complex* data, float* weight,
                         const casa::Bool* flag);

      // Calculate the beam for the given sky direction and frequencies.
      // Apply it to the data.
      // If apply==False, nothing is done.
      void applyBeam (double time, const Position& pos, bool apply,
                      const casa::Vector<double>& chanFreqs, dcomplex* data,
                      StationResponse::vector3r_t& refdir,
                      StationResponse::vector3r_t& tiledir,
                      StationResponse::matrix22c_t* beamvalues,
                      casa::MDirection::Convert& converter);

      // Convert a direction to ITRF.
      StationResponse::vector3r_t dir2Itrf (const casa::MDirection&, casa::MDirection::Convert&) const;

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      string           itsSourceDBName;
      bool             itsUseModelColumn;
      string           itsParmDBName;
      bool             itsApplyBeam;
      bool             itsOneBeamPerPatch;
      bool             itsUseChannelFreq;
      shared_ptr<BBS::ParmDB> itsParmDB;
      Position         itsPhaseRef;

      string           itsMode;
      uint             itsTStep;

      uint             itsDebugLevel;
      bool             itsDetectStalling;
      string           itsStefcalVariant;

      vector<Baseline> itsBaselines;
      vector<ThreadPrivateStorage> itsThreadStorage;

      casa::Array<casa::DComplex> itsVis;
      casa::Array<casa::DComplex> itsMVis;

      vector<casa::Matrix<casa::DComplex> > itsSols; // for every timeslot, nSt gains with vector of length nCr values
      vector<vector<int> > itsAntUseds;
      vector<vector<int> > itsAntMaps;

      StefVecs         iS;

      casa::Vector<casa::String> itsAntennaUsedNames;
      casa::Vector<uint>     itsDataPerAntenna;
      map<string,int>  itsParmIdMap; //# -1 = new parm name

      uint             itsMaxIter;
      double           itsTolerance;
      bool             itsPropagateSolutions;
      uint             itsSolInt;
      uint             itsMinBLperAnt;      

      //# The info needed to calculate the station beams.
      vector<StationResponse::Station::Ptr> itsAntBeamInfo;

      PatchList        itsPatchList;

      string           itsOperation;

      uint             itsConverged;
      uint             itsNonconverged;
      uint             itsStalled;
      uint             itsNTimes;
      NSTimer          itsTimer;
      NSTimer          itsTimerPredict;
      NSTimer          itsTimerSolve;
      NSTimer          itsTimerWrite;
      NSTimer          itsTimerFill;
    };

  } //# end namespace
}

#endif
