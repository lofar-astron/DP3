// GainCal.cc: DPPP step class to predict visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Predict.h"

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"

#include "../ParmDB/ParmDBMeta.h"
#include "../ParmDB/PatchInfo.h"

#include "DPInfo.h"
#include "Exceptions.h"
#include "FlagCounter.h"
#include "Position.h"
#include "ApplyBeam.h"
#include "Simulate.h"
#include "Simulator.h"
#include "Stokes.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include "../ParmDB/SourceDB.h"

#include <aocommon/threadpool.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/RefRows.h>

#include <boost/make_unique.hpp>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>

using namespace casacore;
using namespace DP3::BBS;

namespace DP3 {
namespace DPPP {

Predict::Predict(DPInput* input, const ParameterSet& parset,
                 const string& prefix)
    : itsThreadPool(nullptr), itsMeasuresMutex(nullptr) {
  init(input, parset, prefix,
       parset.getStringVector(prefix + "sources", vector<string>()));
}

Predict::Predict(DPInput* input, const ParameterSet& parset,
                 const string& prefix, const vector<string>& sourcePatterns)
    : itsThreadPool(nullptr), itsMeasuresMutex(nullptr) {
  init(input, parset, prefix, sourcePatterns);
}

void Predict::init(DPInput* input, const ParameterSet& parset,
                   const string& prefix, const vector<string>& sourcePatterns) {
  itsInput = input;
  itsName = prefix;
  itsSourceDBName = parset.getString(prefix + "sourcedb");
  setOperation(parset.getString(prefix + "operation", "replace"));
  itsApplyBeam = parset.getBool(prefix + "usebeammodel", false);
  itsDebugLevel = parset.getInt(prefix + "debuglevel", 0);
  itsPatchList = vector<Patch::ConstPtr>();

  if (!File(itsSourceDBName).exists())
    throw std::runtime_error("Specified source DB name does not exist");
  BBS::SourceDB sourceDB(BBS::ParmDBMeta("", itsSourceDBName), false);

  // Save directions specifications to pass to applycal
  stringstream ss;
  ss << sourcePatterns;
  itsDirectionsStr = ss.str();

  vector<string> patchNames;
  try {
    patchNames = makePatchList(sourceDB, sourcePatterns);
    itsPatchList = makePatches(sourceDB, patchNames, patchNames.size());
    if (itsPatchList.empty()) {
      throw Exception("Couldn't find patch for direction " + itsDirectionsStr);
    }
  } catch (std::exception& exception) {
    throw std::runtime_error(std::string("Something went wrong while reading "
                                         "the source model. The error was: ") +
                             exception.what());
  }

  if (itsApplyBeam) {
    itsUseChannelFreq = parset.getBool(prefix + "usechannelfreq", true);
    itsOneBeamPerPatch = parset.getBool(prefix + "onebeamperpatch", false);
    itsBeamProximityLimit =
        parset.getDouble(prefix + "beamproximitylimit", 60.0) *
        (M_PI / (180.0 * 60.0 * 60.0));

    string mode =
        boost::to_lower_copy(parset.getString(prefix + "beammode", "default"));
    if (mode == "default") {
      itsBeamMode = FullBeamCorrection;
    } else if (mode == "array_factor") {
      itsBeamMode = ArrayFactorBeamCorrection;
    } else if (mode == "element") {
      itsBeamMode = ElementBeamCorrection;
    } else {
      throw Exception("Beammode should be DEFAULT, ARRAY_FACTOR or ELEMENT");
    }

    string element_model = boost::to_lower_copy(
        parset.getString(prefix + "elementmodel", "hamaker"));
    if (element_model == "hamaker") {
      itsElementResponseModel = everybeam::ElementResponseModel::kHamaker;
    } else if (element_model == "lobes") {
      itsElementResponseModel = everybeam::ElementResponseModel::kLOBES;
    } else if (element_model == "oskar") {
      itsElementResponseModel =
          everybeam::ElementResponseModel::kOSKARSphericalWave;
    } else if (element_model == "oskardipole") {
      itsElementResponseModel = everybeam::ElementResponseModel::kOSKARDipole;
    } else {
      throw Exception(
          "Elementmodel should be HAMAKER, LOBES, OSKAR or OSKARDIPOLE");
    }

    // By default, a source model has each direction in one patch. Therefore,
    // if one-beam-per-patch is requested, we don't have to do anything.
    if (!itsOneBeamPerPatch) {
      if (itsBeamProximityLimit > 0.0) {
        // Rework patch list to cluster proximate sources
        itsPatchList =
            clusterProximateSources(itsPatchList, itsBeamProximityLimit);
      } else {
        // Rework patch list to contain a patch for every source
        itsPatchList = makeOnePatchPerComponent(itsPatchList);
      }
    }
  }

  // If called from h5parmpredict, applycal gets set by that step,
  // so must not be read from parset
  if (parset.isDefined(prefix + "applycal.parmdb") ||
      parset.isDefined(prefix + "applycal.steps")) {
    setApplyCal(input, parset, prefix + "applycal.");
  } else {
    itsDoApplyCal = false;
  }

  itsSourceList = makeSourceList(itsPatchList);

  // Determine whether any sources are polarized. If not, enable Stokes-I-
  // only mode (note that this mode cannot be used with itsApplyBeam)
  if (itsApplyBeam && itsBeamMode != ArrayFactorBeamCorrection) {
    itsStokesIOnly = false;
  } else {
    itsStokesIOnly = !checkPolarized(sourceDB, patchNames, patchNames.size());
  }
}

void Predict::setApplyCal(DPInput* input, const ParameterSet& parset,
                          const string& prefix) {
  itsDoApplyCal = true;
  itsApplyCalStep = ApplyCal(input, parset, prefix, true, itsDirectionsStr);
  if (itsOperation != "replace" &&
      parset.getBool(prefix + "applycal.updateweights", false))
    throw std::invalid_argument(
        "Weights cannot be updated when operation is not replace");
  itsResultStep = std::make_shared<ResultStep>();
  itsApplyCalStep.setNextStep(itsResultStep);
}

Predict::~Predict() {}

void Predict::initializeThreadData() {
  const size_t nBl = info().nbaselines();
  const size_t nSt = info().nantenna();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();
  const size_t nThreads = getInfo().nThreads();

  itsUVW.resize(3, nSt);
  itsUVWSplitIndex =
      nsetupSplitUVW(info().nantenna(), info().getAnt1(), info().getAnt2());

  itsModelVis.resize(nThreads);
  itsModelVisPatch.resize(nThreads);
  itsBeamValues.resize(nThreads);
  itsBeamValuesSingle.resize(nThreads);
  itsAntBeamInfo.resize(nThreads);
  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  itsMeasConverters.resize(nThreads);
  itsMeasFrames.resize(nThreads);

  for (unsigned int thread = 0; thread < nThreads; ++thread) {
    if (itsStokesIOnly) {
      itsModelVis[thread].resize(1, nCh, nBl);
    } else {
      itsModelVis[thread].resize(nCr, nCh, nBl);
    }
    bool needMeasConverters = itsMovingPhaseRef;
    needMeasConverters = needMeasConverters || itsApplyBeam;
    if (needMeasConverters) {
      // Prepare measures converters
      itsMeasFrames[thread].set(info().arrayPosCopy());
      itsMeasFrames[thread].set(
          MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
      itsMeasConverters[thread].set(
          MDirection::J2000,
          MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
    }
    if (itsApplyBeam) {
      if (itsStokesIOnly) {
        itsModelVisPatch[thread].resize(1, nCh, nBl);
      } else {
        itsModelVisPatch[thread].resize(nCr, nCh, nBl);
      }
      itsBeamValues[thread].resize(nSt * nCh);
      itsBeamValuesSingle[thread].resize(nSt * nCh);
      itsInput->fillBeamInfo(itsAntBeamInfo[thread], info().antennaNames(),
                             itsElementResponseModel);
    }
  }
}

void Predict::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  info().setBeamCorrectionMode(NoBeamCorrection);

  const size_t nBl = info().nbaselines();
  for (size_t i = 0; i != nBl; ++i) {
    itsBaselines.push_back(Baseline(info().getAnt1()[i], info().getAnt2()[i]));
  }

  try {
    MDirection dirJ2000(
        MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
    Quantum<Vector<Double> > angles = dirJ2000.getAngle();
    itsMovingPhaseRef = false;
    itsPhaseRef = Position(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  } catch (AipsError&) {
    // Phase direction (in J2000) is time dependent
    itsMovingPhaseRef = true;
  }

  initializeThreadData();

  if (itsDoApplyCal) {
    info() = itsApplyCalStep.setInfo(info());
  }
}

std::pair<double, double> Predict::getFirstDirection() const {
  std::pair<double, double> res;

  res.first = itsPatchList[0]->position()[0];
  res.second = itsPatchList[0]->position()[1];
  return res;
}

void Predict::setOperation(const std::string& operation) {
  itsOperation = operation;
  if (itsOperation != "replace" && itsOperation != "add" &&
      itsOperation != "subtract")
    throw std::invalid_argument(
        "Operation must be 'replace', 'add' or 'subtract'.");
}

void Predict::show(std::ostream& os) const {
  os << "Predict " << itsName << '\n';
  os << "  sourcedb:           " << itsSourceDBName << '\n';
  os << "   number of patches: " << itsPatchList.size() << '\n';
  os << "   number of sources: " << itsSourceList.size() << '\n';
  os << "   all unpolarized:   " << boolalpha << itsStokesIOnly << '\n';
  os << "  apply beam:         " << boolalpha << itsApplyBeam << '\n';
  if (itsApplyBeam) {
    os << "   mode:              ";
    if (itsBeamMode == FullBeamCorrection)
      os << "default";
    else if (itsBeamMode == ArrayFactorBeamCorrection)
      os << "array_factor";
    else
      os << "element";
    os << '\n';
    os << "   use channelfreq:   " << boolalpha << itsUseChannelFreq << '\n';
    os << "   one beam per patch:" << boolalpha << itsOneBeamPerPatch << '\n';
    os << "   beam proximity lim:" << itsBeamProximityLimit << '\n';
  }
  os << "  operation:          " << itsOperation << '\n';
  os << "  threads:            " << getInfo().nThreads() << '\n';
  if (itsDoApplyCal) {
    itsApplyCalStep.show(os);
  }
}

void Predict::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " Predict " << itsName << '\n';
}

bool Predict::process(const DPBuffer& bufin) {
  itsTimer.start();
  itsTempBuffer.copy(bufin);
  itsInput->fetchUVW(bufin, itsTempBuffer, itsTimer);
  itsInput->fetchWeights(bufin, itsTempBuffer, itsTimer);

  // Determine the various sizes.
  // const size_t nDr = itsPatchList.size();
  const size_t nSt = info().nantenna();
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();
  const size_t nSamples = nBl * nCh * nCr;

  itsTimerPredict.start();

  nsplitUVW(itsUVWSplitIndex, itsBaselines, itsTempBuffer.getUVW(), itsUVW);

  double time = itsTempBuffer.getTime();
  // Set up directions for beam evaluation
  everybeam::vector3r_t refdir, tiledir;

  bool needMeasConverters = itsMovingPhaseRef;
  needMeasConverters = needMeasConverters || itsApplyBeam;
  if (needMeasConverters) {
    // Because multiple predict steps might be predicting simultaneously, and
    // Casacore is not thread safe, this needs synchronization.
    std::unique_lock<std::mutex> lock;
    if (itsMeasuresMutex != nullptr)
      lock = std::unique_lock<std::mutex>(*itsMeasuresMutex);
    for (unsigned int thread = 0; thread != getInfo().nThreads(); ++thread) {
      itsMeasFrames[thread].resetEpoch(
          MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
      // Do a conversion on all threads
      refdir = dir2Itrf(info().delayCenter(), itsMeasConverters[thread]);
      tiledir = dir2Itrf(info().tileBeamDir(), itsMeasConverters[thread]);
    }
  }

  if (itsMovingPhaseRef) {
    // Convert phase reference to J2000
    MDirection dirJ2000(MDirection::Convert(
        info().phaseCenter(),
        MDirection::Ref(MDirection::J2000, itsMeasFrames[0]))());
    Quantum<Vector<Double> > angles = dirJ2000.getAngle();
    itsPhaseRef = Position(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  }

  std::unique_ptr<aocommon::ThreadPool> localThreadPool;
  aocommon::ThreadPool* pool = itsThreadPool;
  if (pool == nullptr) {
    // If no ThreadPool was specified, we create a temporary one just
    // for executation of this part.
    localThreadPool =
        boost::make_unique<aocommon::ThreadPool>(info().nThreads());
    pool = localThreadPool.get();
  } else {
    if (pool->NThreads() != info().nThreads())
      throw std::runtime_error(
          "Thread pool has inconsistent number of threads!");
  }
  std::vector<Simulator> simulators;
  simulators.reserve(pool->NThreads());
  for (size_t thread = 0; thread != pool->NThreads(); ++thread) {
    itsModelVis[thread] = dcomplex();
    itsModelVisPatch[thread] = dcomplex();

    // When applying beam, simulate into patch vector
    Cube<dcomplex>& simulatedest =
        (itsApplyBeam ? itsModelVisPatch[thread] : itsModelVis[thread]);
    simulators.emplace_back(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
                            info().chanFreqs(), itsUVW, simulatedest,
                            itsStokesIOnly);
  }
  std::vector<Patch::ConstPtr> curPatches(pool->NThreads());

  pool->For(0, itsSourceList.size(), [&](size_t iter, size_t thread) {
    // Keep on predicting, only apply beam when an entire patch is done
    Patch::ConstPtr& curPatch = curPatches[thread];
    if (itsApplyBeam && curPatch != itsSourceList[iter].second &&
        curPatch != nullptr) {
      if (itsStokesIOnly) {
        addBeamToData(curPatch, time, refdir, tiledir, thread, nSamples / nCr,
                      itsModelVisPatch[thread].data(), itsStokesIOnly);
      } else {
        addBeamToData(curPatch, time, refdir, tiledir, thread, nSamples,
                      itsModelVisPatch[thread].data(), itsStokesIOnly);
      }
    }
    simulators[thread].simulate(itsSourceList[iter].first);
    curPatch = itsSourceList[iter].second;
  });
  // Apply beam to the last patch
  pool->For(0, pool->NThreads(), [&](size_t thread, size_t) {
    if (itsApplyBeam && curPatches[thread] != nullptr) {
      if (itsStokesIOnly) {
        addBeamToData(curPatches[thread], time, refdir, tiledir, thread,
                      nSamples / nCr, itsModelVisPatch[thread].data(),
                      itsStokesIOnly);
      } else {
        addBeamToData(curPatches[thread], time, refdir, tiledir, thread,
                      nSamples, itsModelVisPatch[thread].data(),
                      itsStokesIOnly);
      }
    }
  });

  // Add all thread model data to one buffer
  itsTempBuffer.getData() = Complex();
  Complex* tdata = itsTempBuffer.getData().data();
  for (unsigned int thread = 0; thread < pool->NThreads(); ++thread) {
    if (itsStokesIOnly) {
      for (unsigned int i = 0, j = 0; i < nSamples; i += nCr, j++) {
        tdata[i] += itsModelVis[thread].data()[j];
        tdata[i + nCr - 1] += itsModelVis[thread].data()[j];
      }
    } else {
      std::transform(tdata, tdata + nSamples, itsModelVis[thread].data(), tdata,
                     std::plus<dcomplex>());
    }
  }

  // Call ApplyCal step
  if (itsDoApplyCal) {
    itsApplyCalStep.process(itsTempBuffer);
    itsTempBuffer = itsResultStep->get();
    tdata = itsTempBuffer.getData().data();
  }

  // Put predict result from temp buffer into the 'real' buffer
  if (itsOperation == "replace") {
    itsBuffer = itsTempBuffer;
  } else {
    itsBuffer.copy(bufin);
    Complex* data = itsBuffer.getData().data();
    if (itsOperation == "add") {
      std::transform(data, data + nSamples, tdata, data, std::plus<dcomplex>());
    } else if (itsOperation == "subtract") {
      std::transform(data, data + nSamples, tdata, data,
                     std::minus<dcomplex>());
    }
  }

  itsTimerPredict.stop();

  itsTimer.stop();
  getNextStep()->process(itsBuffer);
  return false;
}

everybeam::vector3r_t Predict::dir2Itrf(const MDirection& dir,
                                        MDirection::Convert& measConverter) {
  const MDirection& itrfDir = measConverter(dir);
  const Vector<Double>& itrf = itrfDir.getValue().getValue();
  everybeam::vector3r_t vec;
  vec[0] = itrf[0];
  vec[1] = itrf[1];
  vec[2] = itrf[2];
  return vec;
}

void Predict::addBeamToData(Patch::ConstPtr patch, double time,
                            const everybeam::vector3r_t& refdir,
                            const everybeam::vector3r_t& tiledir,
                            unsigned int thread, unsigned int nSamples,
                            dcomplex* data0, bool stokesIOnly) {
  // Apply beam for a patch, add result to itsModelVis
  MDirection dir(MVDirection(patch->position()[0], patch->position()[1]),
                 MDirection::J2000);
  everybeam::vector3r_t srcdir = dir2Itrf(dir, itsMeasConverters[thread]);

  float* dummyweight = 0;

  if (stokesIOnly) {
    ApplyBeam::applyBeamStokesIArrayFactor(
        info(), time, data0, dummyweight, srcdir, refdir, tiledir,
        itsAntBeamInfo[thread], itsBeamValuesSingle[thread], itsUseChannelFreq,
        false, itsBeamMode, false);
  } else {
    ApplyBeam::applyBeam(info(), time, data0, dummyweight, srcdir, refdir,
                         tiledir, itsAntBeamInfo[thread], itsBeamValues[thread],
                         itsUseChannelFreq, false, itsBeamMode, false);
  }

  // Add temporary buffer to itsModelVis
  std::transform(itsModelVis[thread].data(),
                 itsModelVis[thread].data() + nSamples, data0,
                 itsModelVis[thread].data(), std::plus<dcomplex>());
  itsModelVisPatch[thread] = dcomplex();
}

void Predict::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace DPPP
}  // namespace DP3
