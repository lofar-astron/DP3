// DemixInfo.cc: Struct to hold the common demix variables
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DemixInfo.h"
#include "Exceptions.h"
#include "PointSource.h"
#include "GaussianSource.h"
#include "Stokes.h"
#include "Simulate.h"

#include "../parmdb/SourceDB.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>

using dp3::common::operator<<;

using namespace casacore;

namespace dp3 {
namespace base {

DemixInfo::DemixInfo(const common::ParameterSet& parset, const string& prefix)
    : itsSelBL(parset, prefix, false, "cross"),
      itsSelBLTarget(parset, prefix + "target.", false, "cross", "CS*&"),
      itsPredictModelName(parset.getString(prefix + "estimate.skymodel", "")),
      itsDemixModelName(parset.getString(prefix + "ateam.skymodel")),
      itsTargetModelName(parset.getString(prefix + "target.skymodel")),
      itsSourceNames(parset.getStringVector(prefix + "sources")),
      itsRatio1(parset.getDouble(prefix + "ratio1", 5.)),
      itsRatio2(parset.getDouble(prefix + "ratio2", 0.25)),
      // Default thresholds depend on freq, so filled by function update.
      itsAteamAmplThreshold(parset.getDouble(prefix + "ateam.threshold", 0.)),
      itsTargetAmplThreshold(parset.getDouble(prefix + "target.threshold", 0.)),
      itsAngdistThreshold(parset.getDouble(prefix + "distance.threshold", 60)),
      itsAngdistRefFreq(parset.getDouble(prefix + "distance.reffreq", 60e6)),
      itsDefaultGain(parset.getDouble(prefix + "defaultgain", 1e-3)),
      itsPropagateSolution(
          parset.getBool(prefix + "propagatesolutions", false)),
      itsApplyBeam(parset.getBool(prefix + "applybeam", true)),
      itsSolveBoth(parset.getBool(prefix + "solveboth", false)),
      itsDoSubtract(parset.getBool(prefix + "subtract", true)),
      itsTargetHandling(parset.getUint(prefix + "targethandling", 0)),
      itsVerbose(parset.getUint(prefix + "verbose", 0)),
      itsMaxIter(parset.getUint(prefix + "maxiter", 50)),
      itsMinNBaseline(parset.getUint(prefix + "minnbaseline", 6)),
      itsMinNStation(parset.getUint(prefix + "minnstation", 5)),
      itsNStation(0),
      itsNBl(0),
      itsNCorr(0),
      itsNChanIn(0),
      itsNChanAvgSubtr(parset.getUint(prefix + "freqstep", 1)),
      itsNChanAvg(parset.getUint(prefix + "demixfreqstep", itsNChanAvgSubtr)),
      itsNChanOutSubtr(0),
      itsNChanOut(0),
      itsNTimeAvgSubtr(parset.getUint(prefix + "timestep", 1)),
      itsNTimeAvg(parset.getUint(prefix + "demixtimestep", itsNTimeAvgSubtr)),
      itsChunkSize(parset.getUint(prefix + "chunksize", itsNTimeAvg)),
      itsNTimeChunk(parset.getUint(prefix + "ntimechunk", 0)),
      itsTimeIntervalAvg(0) {
  // Get delta in arcsec and take cosine of it (convert to radians first).
  double delta = parset.getDouble(prefix + "target.delta", 60.);
  itsCosTargetDelta = cos(delta / 3600. * casacore::C::pi / 180.);
  if (itsTargetModelName.empty() || itsDemixModelName.empty())
    throw Exception("An empty name is given for a sky model");
  // If the estimate source model is given, read it.
  if (!itsPredictModelName.empty()) {
    itsAteamList = makePatchList(itsPredictModelName, itsSourceNames);
    // Use all predict patch names if no source names given.
    // In this way we're sure both A-team lists have the same sources
    // in the same order.
    if (itsSourceNames.empty()) {
      itsSourceNames.reserve(itsAteamList.size());
      for (size_t i = 0; i < itsAteamList.size(); ++i) {
        itsSourceNames.push_back(itsAteamList[i]->name());
      }
    }
  }
  itsAteamDemixList = makePatchList(itsDemixModelName, itsSourceNames);
  if (itsTargetHandling != 3) {
    itsTargetList = makePatchList(itsTargetModelName, vector<string>());
  }
  // If no estimate model is given, use the demix model.
  if (itsAteamList.empty()) {
    itsAteamList = itsAteamDemixList;
  }
  if (itsSourceNames.empty()) {
    itsSourceNames.reserve(itsAteamList.size());
    for (size_t i = 0; i < itsAteamList.size(); ++i) {
      itsSourceNames.push_back(itsAteamList[i]->name());
    }
  }
  // Note that the A-team models are in the same order of name.
  // Check they have matching positions.
  assert(itsAteamList.size() == itsAteamDemixList.size());
  for (size_t i = 0; i < itsAteamList.size(); ++i) {
    assert(itsAteamList[i]->name() == itsAteamDemixList[i]->name());
    if (!testAngDist(itsAteamDemixList[i]->direction().ra,
                     itsAteamDemixList[i]->direction().dec,
                     itsAteamList[i]->direction().ra,
                     itsAteamList[i]->direction().dec, itsCosTargetDelta))
      throw Exception("Direction mismatch of source " +
                      itsAteamList[i]->name() + " in A-team SourceDBs ([" +
                      itsAteamDemixList[i]->direction().ra + ", " +
                      itsAteamDemixList[i]->direction().dec + "] and [" +
                      itsAteamList[i]->direction().ra + ", " +
                      itsAteamList[i]->direction().dec + "])");
  }
  makeTargetDemixList();
}

void DemixInfo::makeTargetDemixList() {
  // Get all A-team models for demixing.
  // Note that in the constructor only some sources were read.
  // First open the SourceDB.
  parmdb::SourceDB sdb(parmdb::ParmDBMeta(string(), itsDemixModelName), false,
                       false);
  sdb.lock();
  vector<Patch::ConstPtr> patchList =
      makePatchList(itsDemixModelName, vector<string>());
  // The demix target list is the same as the predict list, but A-team
  // sources must be replaced with their demix model.
  // Also these sources must be removed from the A-team model.
  vector<Patch::ConstPtr> targetDemixList;
  unsigned int ncomponent = 0;
  itsTargetDemixList.reserve(itsTargetList.size());
  for (size_t i = 0; i < itsTargetList.size(); ++i) {
    itsTargetDemixList.push_back(itsTargetList[i]);
    ncomponent += itsTargetList[i]->nComponents();
    // Look if an A-team source matches this target source.
    for (size_t j = 0; j < patchList.size(); ++j) {
      if (testAngDist(itsTargetList[i]->direction().ra,
                      itsTargetList[i]->direction().dec,
                      patchList[j]->direction().ra,
                      patchList[j]->direction().dec, itsCosTargetDelta)) {
        // Match, so use the detailed A-team model.
        itsTargetDemixList[i] = patchList[j];
        ncomponent +=
            (patchList[j]->nComponents() - itsTargetList[i]->nComponents());
        itsTargetReplaced.push_back(patchList[j]->name());
        // A-source is in target, so remove from A-team models (if in there).
        for (size_t k = 0; k < itsAteamList.size(); ++k) {
          if (testAngDist(itsTargetDemixList[i]->direction().ra,
                          itsTargetDemixList[i]->direction().dec,
                          itsAteamList[k]->direction().ra,
                          itsAteamList[k]->direction().dec,
                          itsCosTargetDelta)) {
            itsAteamRemoved.push_back(itsAteamList[k]->name());
            itsAteamList.erase(itsAteamList.begin() + k);
            itsAteamDemixList.erase(itsAteamDemixList.begin() + k);
            break;
          }
        }
        break;
      }
    }
  }
}

void DemixInfo::update(const DPInfo& infoSel, DPInfo& info, size_t nThreads) {
  if (itsNTimeChunk == 0) {
    itsNTimeChunk = nThreads;
  }
  // Remove unused antennae and renumber remaining ones.
  itsInfoSel = infoSel;
  itsInfoSel.removeUnusedAnt();
  // Get size info.
  itsNChanIn = infoSel.nchan();
  itsNCorr = infoSel.ncorr();
  if (itsNCorr != 4)
    throw Exception("Demixing requires data with 4 polarizations");
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
  itsUVWSplitIndex = nsetupSplitUVW(itsInfoSel.nantenna(), itsInfoSel.getAnt1(),
                                    itsInfoSel.getAnt2());
  if (itsVerbose > 1) {
    cout << "splitindex=" << itsUVWSplitIndex << '\n';
  }

  // Determine which baselines to use when estimating A-team and target.
  itsSelTarget = itsSelBLTarget.applyVec(infoSel);

  // Form the baselines.
  // the numbering due to unused stations.
  /// Why is that needed for predict/solve?
  for (unsigned int i = 0; i < itsNBl; ++i) {
    itsBaselines.push_back(
        Baseline(itsInfoSel.getAnt1()[i], itsInfoSel.getAnt2()[i]));
  }

  // Adapt averaging to available nr of channels and times.
  // Use a copy of the DPInfo, otherwise it is updated multiple times.
  DPInfo infoDemix(infoSel);
  itsNTimeAvg = std::min(itsNTimeAvg, infoSel.ntime());
  itsNChanAvg = infoDemix.update(itsNChanAvg, itsNTimeAvg);
  itsNChanOut = infoDemix.nchan();
  itsTimeIntervalAvg = infoDemix.timeInterval();
  ///      itsNTimeDemix      = infoDemix.ntime();

  // Update the overall Demixer DPInfo object.
  itsNTimeAvgSubtr = std::min(itsNTimeAvgSubtr, infoSel.ntime());
  itsNChanAvgSubtr = info.update(itsNChanAvgSubtr, itsNTimeAvgSubtr);
  itsNChanOutSubtr = info.nchan();
  if (itsNChanAvg % itsNChanAvgSubtr != 0)
    throw Exception("Demix frequency averaging " + std::to_string(itsNChanAvg) +
                    " must be a multiple of output averaging " +
                    std::to_string(itsNChanAvgSubtr));
  if (itsNTimeAvg % itsNTimeAvgSubtr != 0)
    throw Exception("Demix time averaging " + std::to_string(itsNTimeAvg) +
                    " must be a multiple of output averaging " +
                    std::to_string(itsNTimeAvgSubtr));
  if (itsChunkSize % itsNTimeAvg != 0)
    throw Exception("Demix predict time chunk size " +
                    std::to_string(itsChunkSize) +
                    " must be a multiple of averaging time step " +
                    std::to_string(itsNTimeAvg));
  itsNTimeOut = itsChunkSize / itsNTimeAvg;
  itsNTimeOutSubtr = itsChunkSize / itsNTimeAvgSubtr;
  // Store channel frequencies for the demix and subtract resolutions.
  itsFreqDemix = infoDemix.chanFreqs();
  itsFreqSubtr = info.chanFreqs();

  // Store phase center direction in J2000.
  MDirection dirJ2000(
      MDirection::Convert(infoSel.phaseCenter(), MDirection::J2000)());
  Quantum<Vector<Double>> angles = dirJ2000.getAngle();
  itsPhaseRef = Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);

  // Determine if the minimum distance (scaled with freq) of A-sources
  // to target is within the threshold.
  // First get the target direction (average of its patches).
  parmdb::PatchSumInfo sumInfo(0);
  for (size_t i = 0; i < itsTargetList.size(); ++i) {
    sumInfo.add(itsTargetList[i]->direction().ra,
                itsTargetList[i]->direction().dec, 1.);
  }
  double targetRa = sumInfo.getRa();
  double targetDec = sumInfo.getDec();
  // Determine the minimum distance.
  double minDist = 1e30;
  double freqRatio = info.refFreq() / itsAngdistRefFreq;
  for (size_t i = 0; i < itsAteamList.size(); ++i) {
    double dist = acos(getCosAngDist(itsAteamList[i]->direction().ra,
                                     itsAteamList[i]->direction().dec, targetRa,
                                     targetDec));
    dist *= freqRatio;
    if (verbose() > 10) {
      cout << "Target distance to " << itsAteamList[i]->name() << " = "
           << dist * 180. / C::pi << " deg" << '\n';
    }
    if (dist < minDist) minDist = dist;
  }
  itsIsAteamNearby = cos(minDist) > cos(itsAngdistThreshold * C::pi / 180.);
}

void DemixInfo::show(ostream& os) const {
  os << "  estimate.skymodel:  " << itsPredictModelName << '\n';
  os << "  ateam.skymodel:     " << itsDemixModelName << '\n';
  os << "  target.skymodel:    " << itsTargetModelName << '\n';
  os << "  sources:            " << itsSourceNames << '\n';
  os << "                        " << itsAteamRemoved
     << " removed from A-team model (in target)" << '\n';
  os << "                        " << itsTargetReplaced
     << " replaced in target model (better A-team model)" << '\n';
  os << "  ratio1:             " << itsRatio1 << '\n';
  os << "  ratio2:             " << itsRatio2 << '\n';
  os << "  ateam.threshold:    " << itsAteamAmplThreshold << '\n';
  os << "  target.threshold:   " << itsTargetAmplThreshold << '\n';
  os << "  target.delta:       "
     << acos(itsCosTargetDelta) * 3600. / casacore::C::pi * 180. << " arcsec"
     << '\n';
  os << "  distance.threshold: " << itsAngdistThreshold << " deg" << '\n';
  os << "  distance.reffreq:   " << itsAngdistRefFreq << " Hz" << '\n';
  os << "  minnbaseline:       " << itsMinNBaseline << '\n';
  os << "  minnstation:        " << itsMinNStation << '\n';
  os << "  maxiter:            " << itsMaxIter << '\n';
  os << "  defaultgain:        " << itsDefaultGain << '\n';
  os << "  propagatesolutions: " << (itsPropagateSolution ? "True" : "False")
     << '\n';
  os << "  applybeam:          " << (itsApplyBeam ? "True" : "False") << '\n';
  os << "  solveboth:          " << (itsSolveBoth ? "True" : "False") << '\n';
  os << "  subtract:           " << (itsDoSubtract ? "True" : "False") << '\n';
  os << "  freqstep:           " << itsNChanAvgSubtr << '\n';
  os << "  timestep:           " << itsNTimeAvgSubtr << '\n';
  os << "  demixfreqstep:      " << itsNChanAvg << '\n';
  os << "  demixtimestep:      " << itsNTimeAvg << '\n';
  os << "  chunksize:          " << itsChunkSize << '\n';
  os << "  ntimechunk:         " << itsNTimeChunk << '\n';
  os << "  target estimate";
  itsSelBLTarget.show(os, "    ");
  os << "  demix";
  itsSelBL.show(os, "    ");
}

vector<Patch::ConstPtr> DemixInfo::makePatchList(
    const string& sdbName, const vector<string>& patchNames) {
  // Open the SourceDB.
  parmdb::SourceDB sdb(parmdb::ParmDBMeta(string(), sdbName), false, false);
  sdb.lock();
  // Get all patches from it.
  vector<string> names(sdb.getPatches());
  vector<string>::const_iterator pnamesIter = patchNames.begin();
  vector<string>::const_iterator pnamesEnd = patchNames.end();
  if (patchNames.empty()) {
    pnamesIter = names.begin();
    pnamesEnd = names.end();
  }
  // Create a patch component list for each matching patch name.
  vector<Patch::ConstPtr> patchList;
  patchList.reserve(patchNames.size());
  for (; pnamesIter != pnamesEnd; ++pnamesIter) {
    if (std::find(names.begin(), names.end(), *pnamesIter) == names.end())
      throw Exception("Demixer: sourcename " + *pnamesIter +
                      " not found in SourceDB " + sdbName);
    // Use this patch; get all its sources.
    vector<parmdb::SourceData> patch = sdb.getPatchSourceData(*pnamesIter);
    vector<ModelComponent::Ptr> componentList;
    componentList.reserve(patch.size());
    for (vector<parmdb::SourceData>::const_iterator iter = patch.begin();
         iter != patch.end(); ++iter) {
      const parmdb::SourceData& src = *iter;
      // Fetch direction.
      assert(src.getInfo().getRefType() == "J2000");
      const Direction direction(src.getRa(), src.getDec());

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
        case parmdb::SourceInfo::POINT: {
          source = PointSource::Ptr(new PointSource(direction, stokes));
        } break;
        case parmdb::SourceInfo::GAUSSIAN: {
          GaussianSource::Ptr gauss(new GaussianSource(direction, stokes));
          const double deg2rad = (casacore::C::pi / 180.0);
          gauss->setPositionAngle(src.getOrientation() * deg2rad);
          const double arcsec2rad = (casacore::C::pi / 3600.0) / 180.0;
          gauss->setMajorAxis(src.getMajorAxis() * arcsec2rad);
          gauss->setMinorAxis(src.getMinorAxis() * arcsec2rad);
          source = gauss;
        } break;
        default: {
          throw Exception(
              "Only point sources and Gaussian sources are"
              " supported at this time.");
        }
      }

      // Fetch spectral index attributes (if applicable).
      if (src.getSpectralTerms().size() > 0) {
        bool isLogarithmic = true;
        source->setSpectralTerms(src.getInfo().getSpectralTermsRefFreq(),
                                 isLogarithmic, src.getSpectralTerms().begin(),
                                 src.getSpectralTerms().end());
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
    auto ppatch = std::make_shared<Patch>(*pnamesIter, componentList.begin(),
                                          componentList.end());
    std::vector<parmdb::PatchInfo> patchInfo(sdb.getPatchInfo(-1, *pnamesIter));
    assert(patchInfo.size() == 1);
    // Set the direction and apparent flux of the patch.
    const Direction patchDirection(patchInfo[0].getRa(), patchInfo[0].getDec());
    ppatch->setDirection(patchDirection);
    ppatch->setBrightness(patchInfo[0].apparentBrightness());
    patchList.push_back(std::move(ppatch));
  }
  return patchList;
}

}  // namespace base
}  // namespace dp3
