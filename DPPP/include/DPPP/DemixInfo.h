//# DemixInfo.h: DPPP struct to hold the common demix variables
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

#ifndef DPPP_DEMIXINFO_H
#define DPPP_DEMIXINFO_H

// @file
// @brief DPPP struct to hold the common demix variables

#include <DPPP/DPInfo.h>
#include <DPPP/BaselineSelection.h>
#include <DPPP/Baseline.h>
#include <DPPP/Patch.h>
#include <Common/ParameterSet.h>

#include <casa/Arrays/Vector.h>

namespace LOFAR {
  namespace DPPP {
    // @ingroup NDPPP

    // This struct holds the common demix variables.
    // It can be shared between the parallel DemixWorker objects.

    class DemixInfo
    {
    public:
      // Constructor to read and initialize the values.
      DemixInfo (const ParameterSet&, const string& prefix);

      // Update the info.
      void update (const DPInfo& infoSel, DPInfo& info);

      // Show parameters.
      void show (ostream&) const;

      // Get the DPInfo object.
      const DPInfo& getInfo() const
        { return itsInfoSel; }

      // Get settings.
      //#    0=test 1=include 2=deproject 3=ignore
      uint   targetHandling() const              {return itsTargetHandling;}
      uint   verbose() const                     {return itsVerbose;}
      uint   maxIter() const                     {return itsMaxIter;}
      uint   minNBaseline() const                {return itsMinNBaseline;}
      uint   minNStation() const                 {return itsMinNStation;}
      uint   nstation() const                    {return itsNStation;}
      uint   nbl() const                         {return itsNBl;}
      uint   ncorr() const                       {return itsNCorr;}
      uint   nchanIn() const                     {return itsNChanIn;}
      uint   nchanAvg() const                    {return itsNChanAvg;}
      uint   nchanAvgSubtr() const               {return itsNChanAvgSubtr;}
      uint   nchanOut() const                    {return itsNChanOut;}
      uint   nchanOutSubtr() const               {return itsNChanOutSubtr;}
      uint   ntimeAvg() const                    {return itsNTimeAvg;}
      uint   ntimeAvgSubtr() const               {return itsNTimeAvgSubtr;}
      uint   ntimeOut() const                    {return itsNTimeOut;}
      uint   ntimeOutSubtr() const               {return itsNTimeOutSubtr;}
      uint   ntimeChunk() const                  {return itsNTimeChunk;}
      uint   chunkSize() const                   {return itsChunkSize;}
      double timeIntervalAvg() const             {return itsTimeIntervalAvg;}
      double ratio1() const                      {return itsRatio1;}
      double ratio2() const                      {return itsRatio2;}
      double ateamAmplThreshold() const          {return itsAteamAmplThreshold;}
      double targetAmplThreshold() const         {return itsTargetAmplThreshold;}
      double defaultGain() const                 {return itsDefaultGain;}
      bool   isAteamNearby() const               {return itsIsAteamNearby;}
      bool   propagateSolution() const           {return itsPropagateSolution;}
      bool   applyBeam() const                   {return itsApplyBeam;}
      bool   solveBoth() const                   {return itsSolveBoth;}
      bool   doSubtract() const                  {return itsDoSubtract;}
      const BaselineSelection& selBL() const     {return itsSelBL;}
      const vector<int>& uvwSplitIndex() const   {return itsUVWSplitIndex;}
      const string& predictModelName() const     {return itsPredictModelName;}
      const string& demixModelName() const       {return itsDemixModelName;}
      const string& targetModelName() const      {return itsTargetModelName;}
      const vector<string>& sourceNames() const  {return itsSourceNames;}
      const Position& phaseRef() const           {return itsPhaseRef;}
      const vector<Baseline>& baselines() const  {return itsBaselines;}
      const casa::Vector<bool> selTarget() const {return itsSelTarget;}
      const casa::Vector<double>& freqDemix() const      {return itsFreqDemix;}
      const casa::Vector<double>& freqSubtr() const      {return itsFreqSubtr;}
      const vector<Patch::ConstPtr>& ateamList() const   {return itsAteamList;}
      const vector<Patch::ConstPtr>& targetList() const  {return itsTargetList;}
      const vector<Patch::ConstPtr>& ateamDemixList() const
        {return itsAteamDemixList;}
      const vector<Patch::ConstPtr>& targetDemixList() const
        {return itsTargetDemixList;}

      // Get the baselines.
      const casa::Vector<casa::Int>& getAnt1() const
        { return itsInfoSel.getAnt1(); }
      const casa::Vector<casa::Int>& getAnt2() const
        { return itsInfoSel.getAnt2(); }

      // Get the antenna names and used antennas.
      const casa::Vector<casa::String>& antennaNames() const
        { return itsInfoSel.antennaNames(); }

      // Get cosine of the angular distance between two sky positions.
      static double getCosAngDist (double ra1, double dec1,
                                   double ra2, double dec2)
      {
        return sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2);
      }

      // Test if two positions in the sky are within delta radians.
      static bool testAngDist (double ra1, double dec1,
                               double ra2, double dec2, double cosDelta)
      {
        return getCosAngDist (ra1, dec1, ra2, dec2) >= cosDelta;
      }

    private:
      // Create a list of patches (and components).
      vector<Patch::ConstPtr> makePatchList (const string& sdbName,
                                             const vector<string>& patchNames);

      // Make the target list for demixing with a detailed model for the
      // possible Ateam sources in it.
      void makeTargetDemixList();

      //# Data members.
      DPInfo                  itsInfoSel;
      BaselineSelection       itsSelBL;
      BaselineSelection       itsSelBLTarget;
      vector<int>             itsUVWSplitIndex;
      string                  itsPredictModelName;
      string                  itsDemixModelName;
      string                  itsTargetModelName;
      vector<string>          itsSourceNames;
      double                  itsRatio1;
      double                  itsRatio2;
      double                  itsAteamAmplThreshold;
      double                  itsTargetAmplThreshold;
      double                  itsCosTargetDelta;
      double                  itsAngdistThreshold;
      double                  itsAngdistRefFreq;
      double                  itsDefaultGain;
      bool                    itsIsAteamNearby;
      bool                    itsPropagateSolution;
      bool                    itsApplyBeam;
      bool                    itsSolveBoth;    //# solve if both stat solvable
      bool                    itsDoSubtract;
      uint                    itsTargetHandling;
      uint                    itsVerbose;           //# trace verbosity level
      uint                    itsMaxIter;           //# max #iter in solve
      uint                    itsMinNBaseline;      //# min #baselines for solve
      uint                    itsMinNStation;       //# min #stations for solve
      uint                    itsNStation;
      uint                    itsNBl;
      uint                    itsNCorr;
      uint                    itsNChanIn;
      uint                    itsNChanAvgSubtr;     //# subtract averaging
      uint                    itsNChanAvg;          //# demix averaging
      uint                    itsNChanOutSubtr;
      uint                    itsNChanOut;
      uint                    itsNTimeAvgSubtr;     //# subtract averaging
      uint                    itsNTimeAvg;          //# demix averaging
      uint                    itsNTimeOutSubtr;     //# #output times per chunk
      uint                    itsNTimeOut;          //# #demix times per chunk
      uint                    itsChunkSize;         //# predict time step
      uint                    itsNTimeChunk;        //# nr chunks in parallel
      double                  itsTimeIntervalAvg;
      Position                itsPhaseRef;          //# original phaseref
      vector<Baseline>        itsBaselines;
      casa::Vector<bool>      itsSelTarget;     //# baselines in target estimate
      casa::Vector<double>    itsFreqDemix;
      casa::Vector<double>    itsFreqSubtr;
      vector<Patch::ConstPtr> itsAteamList;
      vector<Patch::ConstPtr> itsTargetList;
      vector<Patch::ConstPtr> itsAteamDemixList;
      vector<Patch::ConstPtr> itsTargetDemixList;
      vector<string>          itsAteamRemoved;
      vector<string>          itsTargetReplaced;
    };

  } //# end namespace
} //# end namespace

#endif
