//# DemixInfo.cc: Struct to hold the common demix variables
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

#include <lofar_config.h>
#include <DPPP/DemixInfo.h>
#include <DPPP/PointSource.h>
#include <DPPP/GaussianSource.h>
#include <DPPP/Stokes.h>
#include <DPPP/Simulate.h>
#include <ParmDB/SourceDB.h>
#include <Common/ParameterSet.h>
#include <Common/StreamUtil.h>
#include <Common/OpenMP.h>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MeasConvert.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    DemixInfo::DemixInfo (const ParameterSet& parset, const string& prefix)
      : itsSelBL            (parset, prefix, false, "cross"),
        itsSelBLTarget      (parset, prefix+"target.", false, "cross", "CS*&"),
        itsPredictModelName (parset.getString(prefix+"estimate.skymodel", "")),
        itsDemixModelName   (parset.getString(prefix+"ateam.skymodel")),
        itsTargetModelName  (parset.getString(prefix+"target.skymodel")),
        itsSourceNames      (parset.getStringVector (prefix+"sources")),
        itsRatio1           (parset.getDouble (prefix+"ratio1", 5.)),
        itsRatio2           (parset.getDouble (prefix+"ratio2", 0.25)),
        // Default thresholds depend on freq, so filled by function update.
        itsAteamAmplThreshold  (parset.getDouble (prefix+"ateam.threshold",
                                                  0.)),
        itsTargetAmplThreshold (parset.getDouble (prefix+"target.threshold",
                                                  0.)),
        itsAngdistThreshold (parset.getDouble (prefix+"distance.threshold", 60)),
        itsAngdistRefFreq   (parset.getDouble (prefix+"distance.reffreq", 60e6)),
        itsDefaultGain      (parset.getDouble (prefix+"defaultgain", 1e-3)),
        itsPropagateSolution(parset.getBool   (prefix+"propagatesolutions",
                                               false)),
        itsApplyBeam        (parset.getBool   (prefix+"applybeam", true)),
        itsSolveBoth        (parset.getBool   (prefix+"solveboth", false)),
        itsDoSubtract       (parset.getBool   (prefix+"subtract", true)),
        itsTargetHandling   (parset.getUint   (prefix+"targethandling", 0)),
        itsVerbose          (parset.getUint   (prefix+"verbose", 0)),
        itsMaxIter          (parset.getUint   (prefix+"maxiter", 50)),
        itsMinNBaseline     (parset.getUint   (prefix+"minnbaseline", 6)),
        itsMinNStation      (parset.getUint   (prefix+"minnstation", 5)),
        itsNStation         (0),
        itsNBl              (0),
        itsNCorr            (0),
        itsNChanIn          (0),
        itsNChanAvgSubtr    (parset.getUint  (prefix+"freqstep", 1)),
        itsNChanAvg         (parset.getUint  (prefix+"demixfreqstep",
                                              itsNChanAvgSubtr)),
        itsNChanOutSubtr    (0),
        itsNChanOut         (0),
        itsNTimeAvgSubtr    (parset.getUint  (prefix+"timestep", 1)),
        itsNTimeAvg         (parset.getUint  (prefix+"demixtimestep",
                                              itsNTimeAvgSubtr)),
        itsChunkSize        (parset.getUint  (prefix+"chunksize",
                                              itsNTimeAvg)),
        itsNTimeChunk       (parset.getUint  (prefix+"ntimechunk",
                                              OpenMP::maxThreads())),
        itsTimeIntervalAvg  (0)
    {
      // Get delta in arcsec and take cosine of it (convert to radians first).
      double delta = parset.getDouble (prefix+"target.delta", 60.);
      itsCosTargetDelta = cos (delta / 3600. * casa::C::pi / 180.);
      ASSERTSTR (!(itsTargetModelName.empty() || itsDemixModelName.empty()),
                 "An empty name is given for a sky model");
      // If the estimate source model is given, read it.
      if (! itsPredictModelName.empty()) {
        itsAteamList = makePatchList (itsPredictModelName, itsSourceNames);
        // Use all predict patch names if no source names given.
        // In this way we're sure both A-team lists have the same sources
        // in the same order.
        if (itsSourceNames.empty()) {
          itsSourceNames.reserve (itsAteamList.size());
          for (size_t i=0; i<itsAteamList.size(); ++i) {
            itsSourceNames.push_back (itsAteamList[i]->name());
          }
        }
      }
      itsAteamDemixList = makePatchList (itsDemixModelName, itsSourceNames);
      itsTargetList     = makePatchList (itsTargetModelName, vector<string>());
      // If no estimate model is given, use the demix model.
      if (itsAteamList.empty()) {
        itsAteamList = itsAteamDemixList;
      }
      if (itsSourceNames.empty()) {
        itsSourceNames.reserve (itsAteamList.size());
        for (size_t i=0; i<itsAteamList.size(); ++i) {
          itsSourceNames.push_back (itsAteamList[i]->name());
        }
      }
      // Note that the A-team models are in the same order of name.
      // Check they have matching positions.
      ASSERT (itsAteamList.size() == itsAteamDemixList.size());
      for (size_t i=0; i<itsAteamList.size(); ++i) {
        ASSERT (itsAteamList[i]->name() == itsAteamDemixList[i]->name());
        ASSERTSTR (testAngDist (itsAteamDemixList[i]->position()[0],
                                itsAteamDemixList[i]->position()[1],
                                itsAteamList[i]->position()[0],
                                itsAteamList[i]->position()[1],
                                itsCosTargetDelta),
                   "Position mismatch of source " << itsAteamList[i]->name()
                   << " in A-team SourceDBs (["
                   << itsAteamDemixList[i]->position()[0] << ", "
                   << itsAteamDemixList[i]->position()[1] << "] and ["
                   << itsAteamList[i]->position()[0] << ", "
                   << itsAteamList[i]->position()[1] << "])");
      }
      makeTargetDemixList();
    }

    void DemixInfo::makeTargetDemixList()
    {
      // Get all A-team models for demixing.
      // Note that in the constructor only some sources were read.
      // First open the SourceDB.
      BBS::SourceDB sdb(BBS::ParmDBMeta(string(), itsDemixModelName));
      sdb.lock();
      vector<Patch::ConstPtr> patchList = makePatchList (itsDemixModelName,
                                                         vector<string>());
      // The demix target list is the same as the predict list, but A-team
      // sources must be replaced with their demix model.
      // Also these sources must be removed from the A-team model.
      vector<Patch::ConstPtr> targetDemixList;
      uint ncomponent = 0;
      itsTargetDemixList.reserve (itsTargetList.size());
      for (size_t i=0; i<itsTargetList.size(); ++i) {
        itsTargetDemixList.push_back (itsTargetList[i]);
        ncomponent += itsTargetList[i]->nComponents();
        // Look if an A-team source matches this target source.
        for (size_t j=0; j<patchList.size(); ++j) {
          if (testAngDist (itsTargetList[i]->position()[0],
                           itsTargetList[i]->position()[1],
                           patchList[j]->position()[0],
                           patchList[j]->position()[1],
                           itsCosTargetDelta)) {
            // Match, so use the detailed A-team model.
            itsTargetDemixList[i] = patchList[j];
            ncomponent += (patchList[j]->nComponents() -
                           itsTargetList[i]->nComponents());
            itsTargetReplaced.push_back (patchList[j]->name());
            // A-source is in target, so remove from A-team models (if in there).
            for (size_t k=0; k<itsAteamList.size(); ++k) {
              if (testAngDist (itsTargetDemixList[i]->position()[0],
                               itsTargetDemixList[i]->position()[1],
                               itsAteamList[k]->position()[0],
                               itsAteamList[k]->position()[1],
                               itsCosTargetDelta)) {
                itsAteamRemoved.push_back (itsAteamList[k]->name());
                itsAteamList.erase (itsAteamList.begin() + k);
                itsAteamDemixList.erase (itsAteamDemixList.begin() + k);
                break;
              }
            }
            break;
          }
        }
      }
    }
 
    void DemixInfo::update (const DPInfo& infoSel, DPInfo& info)
    {
      // Remove unused antennae and renumber remaining ones.
      itsInfoSel = infoSel;
      itsInfoSel.removeUnusedAnt();
      // Get size info.
      itsNChanIn = infoSel.nchan();
      itsNCorr   = infoSel.ncorr();
      ASSERTSTR (itsNCorr==4, "Demixing requires data with 4 polarizations");
      // NB. The number of baselines and stations refer to the number of
      // selected baselines and the number of unique stations participating
      // in the selected baselines.
      itsNBl = infoSel.nbaselines();
      itsNStation = infoSel.antennaUsed().size();

      // The default thresholds depend on frequency.
      if (itsAteamAmplThreshold <= 0) {
        if (info.refFreq() < 100e6) {
          itsAteamAmplThreshold = 50;
        } else {
          itsAteamAmplThreshold = 5;
        }
      }
      if (itsTargetAmplThreshold <= 0) {
        if (info.refFreq() < 100e6) {
          itsTargetAmplThreshold = 200;
        } else {
          itsTargetAmplThreshold = 100;
        }
      }

      // Setup the baseline index vector used to split the UVWs.
      itsUVWSplitIndex = nsetupSplitUVW (itsInfoSel.nantenna(),
                                         itsInfoSel.getAnt1(),
                                         itsInfoSel.getAnt2());
      if (itsVerbose > 1) {
        cout << "splitindex="<<itsUVWSplitIndex<<endl;
      }

      // Determine which baselines to use when estimating A-team and target.
      itsSelTarget = itsSelBLTarget.applyVec (infoSel);

      // Form the baselines.
      // the numbering due to unused stations.
      /// Why is that needed for predict/solve?
      for (uint i=0; i<itsNBl; ++i) {
        itsBaselines.push_back (Baseline(itsInfoSel.getAnt1()[i],
                                         itsInfoSel.getAnt2()[i]));
      }

      // Adapt averaging to available nr of channels and times.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infoDemix(infoSel);
      itsNTimeAvg = std::min (itsNTimeAvg, infoSel.ntime());
      itsNChanAvg = infoDemix.update (itsNChanAvg, itsNTimeAvg);
      itsNChanOut = infoDemix.nchan();
      itsTimeIntervalAvg = infoDemix.timeInterval();
      ///      itsNTimeDemix      = infoDemix.ntime();

      // Update the overall Demixer DPInfo object.
      itsNTimeAvgSubtr = std::min (itsNTimeAvgSubtr, infoSel.ntime());
      itsNChanAvgSubtr = info.update (itsNChanAvgSubtr, itsNTimeAvgSubtr);
      itsNChanOutSubtr = info.nchan();
      ASSERTSTR (itsNChanAvg % itsNChanAvgSubtr == 0,
        "Demix frequency averaging " << itsNChanAvg
        << " must be a multiple of output averaging "
        << itsNChanAvgSubtr);
      ASSERTSTR (itsNTimeAvg % itsNTimeAvgSubtr == 0,
        "Demix time averaging " << itsNTimeAvg
        << " must be a multiple of output averaging "
        << itsNTimeAvgSubtr);
      ASSERTSTR (itsChunkSize % itsNTimeAvg == 0,
        "Demix predict time chunk size " << itsChunkSize
        << " must be a multiple of averaging time step "
        << itsNTimeAvg);
      itsNTimeOut = itsChunkSize / itsNTimeAvg;
      itsNTimeOutSubtr = itsChunkSize / itsNTimeAvgSubtr;
      // Store channel frequencies for the demix and subtract resolutions.
      itsFreqDemix = infoDemix.chanFreqs();
      itsFreqSubtr = info.chanFreqs();

      // Store phase center position in J2000.
      MDirection dirJ2000(MDirection::Convert(infoSel.phaseCenter(),
                                              MDirection::J2000)());
      Quantum<Vector<Double> > angles = dirJ2000.getAngle();
      itsPhaseRef = Position(angles.getBaseValue()[0],
                             angles.getBaseValue()[1]);

      // Determine if the minimum distance (scaled with freq) of A-sources
      // to target is within the threshold.
      // First get the target position (average of its patches).
      BBS::PatchSumInfo sumInfo(0);
      for (size_t i=0; i<itsTargetList.size(); ++i) {
        sumInfo.add (itsTargetList[i]->position()[0],
                     itsTargetList[i]->position()[1],
                     1.);
      }
      double targetRa  = sumInfo.getRa();
      double targetDec = sumInfo.getDec();
      // Determine the minimum distance.
      double minDist   = 1e30;
      double freqRatio = info.refFreq() / itsAngdistRefFreq;
      for (size_t i=0; i<itsAteamList.size(); ++i) {
        double dist = acos (getCosAngDist (itsAteamList[i]->position()[0],
                                           itsAteamList[i]->position()[1],
                                           targetRa, targetDec));
        dist *= freqRatio;
        if (verbose() > 10) {
          cout << "Target distance to " << itsAteamList[i]->name()
               << " = " << dist*180./C::pi << " deg" << endl;
        }
        if (dist < minDist) minDist = dist;
      }
      itsIsAteamNearby = cos(minDist) > cos(itsAngdistThreshold*C::pi/180.);
    }

    void DemixInfo::show (ostream& os) const
    {
      os << "  estimate.skymodel:  " << itsPredictModelName << endl;
      os << "  ateam.skymodel:     " << itsDemixModelName << endl;
      os << "  target.skymodel:    " << itsTargetModelName << endl;
      os << "  sources:            " << itsSourceNames << endl;
      os << "                        " << itsAteamRemoved
         << " removed from A-team model (in target)" << endl;
      os << "                        " << itsTargetReplaced
         << " replaced in target model (better A-team model)" << endl;
      os << "  ratio1:             " << itsRatio1 << endl;
      os << "  ratio2:             " << itsRatio2 << endl;
      os << "  ateam.threshold:    " << itsAteamAmplThreshold << endl;
      os << "  target.threshold:   " << itsTargetAmplThreshold << endl;
      os << "  target.delta:       "
         << acos(itsCosTargetDelta) * 3600. / casa::C::pi * 180.
         << " arcsec" << endl;
      os << "  distance.threshold: " << itsAngdistThreshold << " deg" << endl;
      os << "  distance.reffreq:   " << itsAngdistRefFreq << " Hz" << endl;
      os << "  minnbaseline:       " << itsMinNBaseline << endl;
      os << "  minnstation:        " << itsMinNStation << endl;
      os << "  maxiter:            " << itsMaxIter << endl;
      os << "  defaultgain:        " << itsDefaultGain << endl;
      os << "  propagatesolutions: " << (itsPropagateSolution ? "True":"False")
         << endl;
      os << "  applybeam:          " << (itsApplyBeam ? "True":"False") << endl;
      os << "  solveboth:          " << (itsSolveBoth ? "True":"False") << endl;
      os << "  subtract:           " << (itsDoSubtract ? "True":"False")
         << endl;
      os << "  freqstep:           " << itsNChanAvgSubtr << endl;
      os << "  timestep:           " << itsNTimeAvgSubtr << endl;
      os << "  demixfreqstep:      " << itsNChanAvg << endl;
      os << "  demixtimestep:      " << itsNTimeAvg << endl;
      os << "  chunksize:          " << itsChunkSize << endl;
      os << "  ntimechunk:         " << itsNTimeChunk << endl;
      os << "  target estimate";
      itsSelBLTarget.show (os, "    ");
      os << "  demix";
      itsSelBL.show (os, "    ");
    }

    vector<Patch::ConstPtr>
    DemixInfo::makePatchList (const string& sdbName,
                              const vector<string>& patchNames)
    {
      // Open the SourceDB.
      BBS::SourceDB sdb(BBS::ParmDBMeta(string(), sdbName));
      sdb.lock();
      // Get all patches from it.
      vector<string> names(sdb.getPatches());
      vector<string>::const_iterator pnamesIter = patchNames.begin();
      vector<string>::const_iterator pnamesEnd  = patchNames.end();
      if (patchNames.empty()) {
        pnamesIter = names.begin();
        pnamesEnd  = names.end();
      }
      // Create a patch component list for each matching patch name.
      vector<Patch::ConstPtr> patchList;
      patchList.reserve (patchNames.size());
      for (; pnamesIter != pnamesEnd; ++pnamesIter) {
        ASSERTSTR (std::find (names.begin(), names.end(), *pnamesIter)
                   != names.end(),
                   "Demixer: sourcename " << *pnamesIter
                   << " not found in SourceDB " << sdbName);
        // Use this patch; get all its sources.
        vector<BBS::SourceData> patch = sdb.getPatchSourceData (*pnamesIter);
        vector<ModelComponent::Ptr> componentList;
        componentList.reserve (patch.size());
        for (vector<BBS::SourceData>::const_iterator iter=patch.begin();
             iter!=patch.end(); ++iter) {
          const BBS::SourceData& src = *iter;
          // Fetch position.
          ASSERT (src.getInfo().getRefType() == "J2000");
          Position position;
          position[0] = src.getRa();
          position[1] = src.getDec();

          // Fetch stokes vector.
          Stokes stokes;
          stokes.I = src.getI();
          stokes.V = src.getV();
          if (!src.getInfo().getUseRotationMeasure()) {
            stokes.Q = src.getQ();
            stokes.U = src.getU();
          }

          PointSource::Ptr source;
          switch (src.getInfo().getType()) {
          case BBS::SourceInfo::POINT:
            {
              source = PointSource::Ptr(new PointSource(position, stokes));
            }
            break;
          case BBS::SourceInfo::GAUSSIAN:
            {
              GaussianSource::Ptr gauss(new GaussianSource(position, stokes));
              const double deg2rad = (casa::C::pi / 180.0);
              gauss->setPositionAngle(src.getOrientation() * deg2rad);
              const double arcsec2rad = (casa::C::pi / 3600.0) / 180.0;
              gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
              gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
              source = gauss;
            }
            break;
          default:
            {
              ASSERTSTR(false, "Only point sources and Gaussian sources are"
                        " supported at this time.");
            }
          }

          // Fetch spectral index attributes (if applicable).
          if (src.getSpectralIndex().size() > 0) {
            source->setSpectralIndex(src.getInfo().getSpectralIndexRefFreq(),
                                     src.getSpectralIndex().begin(),
                                     src.getSpectralIndex().end());
          }

          // Fetch rotation measure attributes (if applicable).
          if (src.getInfo().getUseRotationMeasure()) {
            source->setRotationMeasure(src.getPolarizedFraction(),
                                       src.getPolarizationAngle(),
                                       src.getRotationMeasure());
          }

          // Add the source definition.
          componentList.push_back(source);
        }

        // Add the component list as a patch to the list of patches.
        Patch::Ptr ppatch(new Patch(*pnamesIter,
                                    componentList.begin(),
                                    componentList.end()));
        vector<BBS::PatchInfo> patchInfo(sdb.getPatchInfo(-1, *pnamesIter));
        ASSERT (patchInfo.size() == 1);
        // Set the position and apparent flux of the patch.
        Position patchPosition;
        patchPosition[0] = patchInfo[0].getRa();
        patchPosition[1] = patchInfo[0].getDec();
        ppatch->setPosition (patchPosition);
        ppatch->setBrightness (patchInfo[0].apparentBrightness());
        patchList.push_back (ppatch);
      }
      return patchList;
    }


  } //# end namespace
} //# end namespace
