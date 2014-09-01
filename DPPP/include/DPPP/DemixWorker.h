//# DemixWorker.h: Demixer helper class processing a time chunk
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
//# $Id: Demixer.h 23223 2012-12-07 14:09:42Z schoenmakers $
//#
//# @author Ger van Diepen

#ifndef DPPP_DEMIXWORKER_H
#define DPPP_DEMIXWORKER_H

// @file
// @brief DPPP step class to average in time and/or freq

#include <DPPP/DemixInfo.h>
#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/Patch.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/Filter.h>
#include <DPPP/EstimateNew.h>
#include <StationResponse/Station.h>
#include <ParmDB/ParmDB.h>

#include <casa/Arrays/Cube.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MeasureHolder.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>

namespace LOFAR {

  namespace DPPP {
    // @ingroup NDPPP

    // DemixWorker::process processes a single time window (say, 2 minutes).
    // It predicts the A-team and target sources to determine which sources
    // have to be taken into account and which antennae have to be solved for.
    // Multiple DemixWorker::process can be executed in parallel by the parent
    // class DemixerNew.
    //
    // Each DemixWorker object references a DemixInfo object containing the
    // general info and parameters.

    class DemixWorker
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      DemixWorker (DPInput*,
                   const string& prefix,
                   const DemixInfo& info,
                   const DPInfo& dpinfo,
		   int workernr);


      // Process the data in the input buffers and store the result in the
      // output buffers.
      void process (const DPBuffer* bufin, uint nbufin,
                    DPBuffer* bufout, vector<double>* solutions,
                    uint chunkNr);

      // Get the number of solves.
      uint nSolves() const
        { return itsNrSolves; }
      // Get the number of converged solves.
      uint nConverged() const
        { return itsNrConverged; }
      // Get the total nr of iterations used.
      uint nIterations() const
        { return itsNrIter; }
      // Get the number of times no demix was needed.
      uint nNoDemix() const
        { return itsNrNoDemix; }
      uint nIncludeStrongTarget() const
        { return itsNrIncludeStrongTarget; }
      uint nIncludeCloseTarget() const
        { return itsNrIncludeCloseTarget; }
      uint nIgnoreTarget() const
        { return itsNrIgnoreTarget; }
      uint nDeprojectTarget() const
        { return itsNrDeprojectTarget; }
      // Get nr of times a source was demixed.
      const casa::Vector<uint>& nsourcesDemixed() const
        { return itsNrSourcesDemixed; }
      // Get nr of times a station was demixed.
      const casa::Vector<uint>& nstationsDemixed() const
        { return itsNrStationsDemixed; }
      // Get nr of times a station/source was demixed.
      const casa::Matrix<uint>& statSourceDemixed() const
        { return itsStatSourceDemixed; }
      const casa::Matrix<double>& amplSubtrMean() const
        { return itsAmplSubtrMean; }
      const casa::Matrix<double>& amplSubtrM2() const
        { return itsAmplSubtrM2; }
      const casa::Matrix<size_t>& amplSubtrNr() const
        { return itsAmplSubtrNr; }

      // Get the timings of the various processing steps.
      // <group>
      double getTotalTime() const
        { return itsTimer.getElapsed(); }
      double getCoarseTime() const
        { return itsTimerCoarse.getElapsed(); }
      double getPhaseShiftTime() const
        { return itsTimerPhaseShift.getElapsed(); }
      double getDemixTime() const
        { return itsTimerDemix.getElapsed(); }
      double getPredictTime() const
        { return itsTimerPredict.getElapsed(); }
      double getSolveTime() const
        { return itsTimerSolve.getElapsed(); }
      double getSubtractTime() const
        { return itsTimerSubtract.getElapsed(); }
      // </group>

    private:
      // Setup the demix processing steps for this piece of data.
      // It fills itsFirstSteps, etc. for the sources to be demixed.
      // It also determines how to handle the target (include,deproject,ignore).
      void setupDemix (uint chunkNr);

      // Find the median ampltitude for the selected baselines.
      // It uses itsTmpAmpl as temporary buffer.
      float findMedian (const casa::Cube<float>& ampl, const bool* selbl);

      // Average the baseline UVWs in bufin and split them into UVW per station.
      // It returns the number of time averages.
      uint avgSplitUVW (const DPBuffer* bufin, uint nbufin,
                        uint ntimeAvg, const vector<uint>& selbl);

      // Predict the target StokesI amplitude.
      // It applies the beam at each target patch.
      void predictTarget (const vector<Patch::ConstPtr>& patchList,
                          uint ntime, double time, double timeStep);

      // Predict the StokesI amplitude of the Ateam patches and determine
      // which antennae and sources to use when demixing.
      // It applies the beam at each patch center.
      void predictAteam (const vector<Patch::ConstPtr>& patchList,
                         uint ntime, double time, double timeStep);

      // Add the StokesI of itsPredictVis to ampl.
      void addStokesI (casa::Matrix<float>& ampl);

      // Calculate the beam for demix resolution and apply to itsPredictVis.
      // If apply==False, nothing is done.
      void applyBeam (double time, const Position& pos, bool apply);

      // Calculate the beam for the given sky direction and frequencies.
      // Apply it to the data.
      // If apply==False, nothing is done.
      void applyBeam (double time, const Position& pos, bool apply,
                      const casa::Vector<double>& chanFreqs,
                      dcomplex* data);

      // Convert a direction to ITRF.
      StationResponse::vector3r_t dir2Itrf (const casa::MDirection&);

      // Calculate the StokesI amplitude from the predicted visibilities.
      // (0.5 * (XX+YY))
      void calcStokesI (casa::Matrix<float>& ampl);

      // Simply average the data if no demixing needs to bedone.
      void average (const DPBuffer* bufin, uint nbufin,
                    DPBuffer* bufout);

      // Add the decorrelation factor contribution for each time slot.
      void addFactors (const DPBuffer& newBuf,
                       casa::Array<casa::DComplex>& factorBuf);

      // Calculate the decorrelation factors by averaging them.
      // Apply the P matrix to deproject the sources without a model.
      void makeFactors (const casa::Array<casa::DComplex>& bufIn,
                        casa::Array<casa::DComplex>& bufOut,
                        const casa::Cube<float>& weightSums,
                        uint nChanOut,
                        uint nChanAvg);

      // Deproject the sources without a model.
      void deproject (casa::Array<casa::DComplex>& factors,
                      vector<MultiResultStep*> avgResults,
                      uint resultIndex);

      // Do the demixing.
      void handleDemix (DPBuffer* bufout, vector<double>* solutions,
                        double time, double timeStep);

      // Solve gains and subtract sources.
      void demix (vector<double>* solutions, double time, double timeStep);

      // Add amplitude subtracted to the arrays for mean and stddev.
      void addMeanM2 (const vector<float>& sourceAmpl, uint src);

      // Merge the data of the selected baselines from the subtract buffer
      // into the full buffer.
      void mergeSubtractResult();

      //# Data members.
      int                                   itsWorkerNr;
      const DemixInfo*                      itsMix;
      vector<PhaseShift*>                   itsOrigPhaseShifts;
      //# Phase shift and average steps for demix.
      vector<DPStep::ShPtr>                 itsOrigFirstSteps;
      //# Result of phase shifting and averaging the directions of interest
      //# at the demix resolution.
      vector<MultiResultStep*>              itsAvgResults;
      vector<PhaseShift*>                   itsPhaseShifts;
      vector<DPStep::ShPtr>                 itsFirstSteps;
      DPStep::ShPtr                         itsAvgStepSubtr;
      Filter                                itsFilter;
      Filter*                               itsFilterSubtr;
      //# Result of averaging the target at the subtract resolution.
      MultiResultStep*                      itsAvgResultFull;
      MultiResultStep*                      itsAvgResultSubtr;
      //# The sources to demix (excluding target).
      vector<Patch::ConstPtr>               itsDemixList;
      //# The info needed to calculate the station beams.
      vector<StationResponse::Station::Ptr> itsAntBeamInfo;
      //# Measure objects unique to this worker (thread).
      //# This is needed because they are not thread-safe.
      casa::MPosition                       itsArrayPos;
      casa::MDirection                      itsDelayCenter;
      casa::MDirection                      itsTileBeamDir;

      //# Variables set by setupDemix and used by handleDemix.
      uint                                  itsNDir;
      uint                                  itsNModel;
      uint                                  itsNSubtr;
      bool                                  itsIgnoreTarget;
      bool                                  itsIncludeTarget;
      //# Accumulator used for computing the demixing weights at the demix
      //# resolution. The shape of this buffer is #correlations x #channels
      //# x #baselines x #directions x #directions (fastest axis first).
      casa::Array<casa::DComplex>           itsFactorBuf;
      //# Buffer of demixing weights at the demix resolution. Each Array is a
      //# cube of shape #correlations x #channels x #baselines of matrices of
      //# shape #directions x #directions.
      vector<casa::Array<casa::DComplex> >  itsFactors;
      //# Accumulator used for computing the demixing weights. The shape of this
      //# buffer is #correlations x #channels x #baselines x #directions
      //# x #directions (fastest axis first).
      casa::Array<casa::DComplex>           itsFactorBufSubtr;
      //# Buffer of demixing weights at the subtract resolution. Each Array is a
      //# cube of shape #correlations x #channels x #baselines of matrices of
      //# shape #directions x #directions.
      vector<casa::Array<casa::DComplex> >  itsFactorsSubtr;

      //# Variables for conversion of directions to ITRF.
      casa::MeasFrame                       itsMeasFrame;
      casa::MDirection::Convert             itsMeasConverter;
      vector<StationResponse::matrix22c_t>  itsBeamValues;  //# [nst,nch]

      //# Indices telling which Ateam sources to use.
      vector<uint>                          itsSrcSet;
      //# UVW per station per demix time slot
      casa::Cube<double>                    itsStationUVW;  //# UVW per station
      casa::Matrix<double>                  itsAvgUVW;      //# temp buffer
      casa::Cube<dcomplex>                  itsPredictVis;  //# temp buffer
      //# #nfreq x #bl x #time StokesI amplitude per A-source.
      vector<casa::Cube<float> >            itsAteamAmpl;
      //# #bl x #src telling if baseline has sufficient Ateam flux.
      casa::Matrix<bool>                    itsAteamAmplSel;
      //# #nfreq x #bl x #time StokesI amplitude of target.
      casa::Cube<float>                     itsTargetAmpl;
      //# Temporary buffer to determine medians.
      vector<float>                         itsTmpAmpl;
      //# Per A-source and for target the min and max amplitude.
      vector<double>                        itsAteamMinAmpl;
      vector<double>                        itsAteamMaxAmpl;
      double                                itsTargetMinAmpl;
      double                                itsTargetMaxAmpl;
      //# Per A-source the stations to use (matching the minimum amplitude).
      vector<vector<uint> >                 itsStationsToUse;
      casa::Block<bool>                     itsSolveStation; //# solve station i?
      //# Per station and source the index in the unknowns vector.
      //# Note there are 8 unknowns (4 pol, ampl/phase) per source/station.
      vector<vector<int> >                  itsUnknownsIndex;
      //# The estimater (solver).
      EstimateNew                           itsEstimate;
      //# Variables for the predict.
      casa::Matrix<double>                  itsUVW;
      vector<dcomplex>                      itsModelVis;
      uint                                  itsNTimeOut;
      uint                                  itsNTimeOutSubtr;
      uint                                  itsTimeIndex;
      vector<float>                         itsObservedAmpl;
      vector<float>                         itsSourceAmpl;
      vector<float>                         itsSumSourceAmpl;
      //# Statistics
      uint                                  itsNrSolves;
      uint                                  itsNrConverged;
      uint                                  itsNrIter;
      uint                                  itsNrNoDemix;
      uint                                  itsNrIncludeStrongTarget;
      uint                                  itsNrIncludeCloseTarget;
      uint                                  itsNrIgnoreTarget;
      uint                                  itsNrDeprojectTarget;
      //# Nr of times a source is demixed.
      casa::Vector<uint>                    itsNrSourcesDemixed;
      //# Nr of times a station is demixed.
      casa::Vector<uint>                    itsNrStationsDemixed;
      //# Nr of times a source/station is demixed.
      casa::Matrix<uint>                    itsStatSourceDemixed;
      //# Average amplitude subtracted for middle channel [nbl,nsrc]
      casa::Matrix<double>                  itsAmplSubtrMean;
      //# M2n to calculate stddev online in stable way (see Wikipedia)
      casa::Matrix<double>                  itsAmplSubtrM2;
      //# N for mean/stddev amplitude calculations.
      casa::Matrix<size_t>                  itsAmplSubtrNr;
      //# Timers.
      NSTimer                               itsTimer;
      NSTimer                               itsTimerCoarse;
      NSTimer                               itsTimerPhaseShift;
      NSTimer                               itsTimerDemix;
      NSTimer                               itsTimerPredict;
      NSTimer                               itsTimerSolve;
      NSTimer                               itsTimerSubtract;
    };

  } //# end namespace
} //# end namespace

#endif
