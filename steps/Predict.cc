// GainCal.cc: DPPP step class to predict visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Predict.h"
#include "ApplyBeam.h"

#include <iostream>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"

#include "../parmdb/ParmDBMeta.h"
#include "../parmdb/PatchInfo.h"

#include "../base/DPInfo.h"
#include "../base/Exceptions.h"
#include "../base/FlagCounter.h"
#include "../base/Position.h"
#include "../base/Simulate.h"
#include "../base/Simulator.h"
#include "../base/Stokes.h"
#include "../base/PointSource.h"
#include "../base/GaussianSource.h"

#include "../parmdb/SourceDB.h"

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

using casacore::Cube;
using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVDirection;
using casacore::MVEpoch;
using casacore::Quantum;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

Predict::Predict(InputStep* input, const common::ParameterSet& parset,
                 const string& prefix)
    : itsThreadPool(nullptr), itsMeasuresMutex(nullptr) {
  init(input, parset, prefix,
       parset.getStringVector(prefix + "sources", std::vector<string>()));
}

Predict::Predict(InputStep* input, const common::ParameterSet& parset,
                 const string& prefix,
                 const std::vector<string>& sourcePatterns)
    : itsThreadPool(nullptr), itsMeasuresMutex(nullptr) {
  init(input, parset, prefix, sourcePatterns);
}

void Predict::init(InputStep* input, const common::ParameterSet& parset,
                   const string& prefix,
                   const std::vector<string>& sourcePatterns) {
  itsInput = input;
  itsName = prefix;
  itsSourceDBName = parset.getString(prefix + "sourcedb");
  setOperation(parset.getString(prefix + "operation", "replace"));
  itsApplyBeam = parset.getBool(prefix + "usebeammodel", false);
  itsDebugLevel = parset.getInt(prefix + "debuglevel", 0);
  itsPatchList = std::vector<base::Patch::ConstPtr>();

  if (!casacore::File(itsSourceDBName).exists())
    throw std::runtime_error("Specified source DB name does not exist");
  parmdb::SourceDB sourceDB(parmdb::ParmDBMeta("", itsSourceDBName), false);

  // Save directions specifications to pass to applycal
  std::stringstream ss;
  ss << sourcePatterns;
  itsDirectionsStr = ss.str();

  std::vector<string> patchNames;
  try {
    patchNames = base::makePatchList(sourceDB, sourcePatterns);
    itsPatchList = base::makePatches(sourceDB, patchNames, patchNames.size());
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
      itsBeamMode = base::FullBeamCorrection;
    } else if (mode == "array_factor") {
      itsBeamMode = base::ArrayFactorBeamCorrection;
    } else if (mode == "element") {
      itsBeamMode = base::ElementBeamCorrection;
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
  if (itsApplyBeam && itsBeamMode != base::ArrayFactorBeamCorrection) {
    itsStokesIOnly = false;
  } else {
    itsStokesIOnly =
        !base::checkPolarized(sourceDB, patchNames, patchNames.size());
  }
}

void Predict::setApplyCal(InputStep* input, const common::ParameterSet& parset,
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
  const size_t nCr = itsStokesIOnly ? 1 : info().ncorr();
  const size_t nThreads = getInfo().nThreads();

  itsUVW.resize(3, nSt);
  itsUVWSplitIndex = base::nsetupSplitUVW(info().nantenna(), info().getAnt1(),
                                          info().getAnt2());

  if (!itsPredictBuffer) {
    itsPredictBuffer = std::make_shared<base::PredictBuffer>();
  }
  if (itsApplyBeam && itsPredictBuffer->GetStationList().empty()) {
    itsInput->fillBeamInfo(itsPredictBuffer->GetStationList(),
                           info().antennaNames(), itsElementResponseModel);
  }
  itsPredictBuffer->resize(nThreads, nCr, nCh, nBl, nSt, itsApplyBeam);
  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  itsMeasConverters.resize(nThreads);
  itsMeasFrames.resize(nThreads);

  for (size_t thread = 0; thread < nThreads; ++thread) {
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
  }
}

void Predict::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  info().setBeamCorrectionMode(base::NoBeamCorrection);

  const size_t nBl = info().nbaselines();
  for (size_t i = 0; i != nBl; ++i) {
    itsBaselines.push_back(Baseline(info().getAnt1()[i], info().getAnt2()[i]));
  }

  try {
    MDirection dirJ2000(
        MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
    Quantum<casacore::Vector<double> > angles = dirJ2000.getAngle();
    itsMovingPhaseRef = false;
    itsPhaseRef =
        base::Position(angles.getBaseValue()[0], angles.getBaseValue()[1]);
  } catch (casacore::AipsError&) {
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
  os << "   all unpolarized:   " << std::boolalpha << itsStokesIOnly << '\n';
  os << "  apply beam:         " << std::boolalpha << itsApplyBeam << '\n';
  if (itsApplyBeam) {
    os << "   mode:              ";
    if (itsBeamMode == base::FullBeamCorrection)
      os << "default";
    else if (itsBeamMode == base::ArrayFactorBeamCorrection)
      os << "array_factor";
    else
      os << "element";
    os << '\n';
    os << "   use channelfreq:   " << std::boolalpha << itsUseChannelFreq
       << '\n';
    os << "   one beam per patch:" << std::boolalpha << itsOneBeamPerPatch
       << '\n';
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
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " Predict " << itsName << '\n';
}

bool Predict::process(const DPBuffer& bufin) {
  itsTimer.start();
  DPBuffer tempBuffer;
  tempBuffer.copy(bufin);
  itsInput->fetchUVW(bufin, tempBuffer, itsTimer);
  itsInput->fetchWeights(bufin, tempBuffer, itsTimer);

  // Determine the various sizes.
  // const size_t nDr = itsPatchList.size();
  const size_t nSt = info().nantenna();
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();
  const size_t nBeamValues = itsStokesIOnly ? nBl * nCh : nBl * nCh * nCr;

  itsTimerPredict.start();

  base::nsplitUVW(itsUVWSplitIndex, itsBaselines, tempBuffer.getUVW(), itsUVW);

  double time = tempBuffer.getTime();
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
    for (size_t thread = 0; thread != getInfo().nThreads(); ++thread) {
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
    Quantum<casacore::Vector<double> > angles = dirJ2000.getAngle();
    itsPhaseRef =
        base::Position(angles.getBaseValue()[0], angles.getBaseValue()[1]);
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
  std::vector<base::Simulator> simulators;
  simulators.reserve(pool->NThreads());
  for (size_t thread = 0; thread != pool->NThreads(); ++thread) {
    itsPredictBuffer->GetModel(thread) = dcomplex();
    if (itsApplyBeam) itsPredictBuffer->GetPatchModel(thread) = dcomplex();

    // When applying beam, simulate into patch vector
    Cube<dcomplex>& simulatedest =
        (itsApplyBeam ? itsPredictBuffer->GetPatchModel(thread)
                      : itsPredictBuffer->GetModel(thread));
    simulators.emplace_back(itsPhaseRef, nSt, nBl, nCh, itsBaselines,
                            info().chanFreqs(), itsUVW, simulatedest,
                            itsStokesIOnly);
  }
  std::vector<base::Patch::ConstPtr> curPatches(pool->NThreads());

  pool->For(0, itsSourceList.size(), [&](size_t sourceIndex, size_t thread) {
    // Predict the source model and apply beam when an entire patch is done
    base::Patch::ConstPtr& curPatch = curPatches[thread];
    const bool patchIsFinished =
        curPatch != itsSourceList[sourceIndex].second && curPatch != nullptr;
    if (itsApplyBeam && patchIsFinished) {
      // Apply the beam and add PatchModel to Model
      addBeamToData(curPatch, time, refdir, tiledir, thread, nBeamValues,
                    itsPredictBuffer->GetPatchModel(thread).data(),
                    itsStokesIOnly);
      // Initialize patchmodel to zero for the next patch
      itsPredictBuffer->GetPatchModel(thread) = dcomplex();
    }
    // Depending on itsApplyBeam, the following call will add to either
    // the Model or the PatchModel of the predict buffer
    simulators[thread].simulate(itsSourceList[sourceIndex].first);

    curPatch = itsSourceList[sourceIndex].second;
  });
  // Apply beam to the last patch
  if (itsApplyBeam) {
    pool->For(0, pool->NThreads(), [&](size_t thread, size_t) {
      if (curPatches[thread] != nullptr) {
        addBeamToData(
            curPatches[thread], time, refdir, tiledir, thread, nBeamValues,
            itsPredictBuffer->GetPatchModel(thread).data(), itsStokesIOnly);
      }
    });
  }

  // Add all thread model data to one buffer
  tempBuffer.getData() = casacore::Complex();
  casacore::Complex* tdata = tempBuffer.getData().data();
  const size_t nVisibilities = nBl * nCh * nCr;
  for (size_t thread = 0; thread < pool->NThreads(); ++thread) {
    if (itsStokesIOnly) {
      for (size_t i = 0, j = 0; i < nVisibilities; i += nCr, j++) {
        tdata[i] += itsPredictBuffer->GetModel(thread).data()[j];
        tdata[i + nCr - 1] += itsPredictBuffer->GetModel(thread).data()[j];
      }
    } else {
      std::transform(tdata, tdata + nVisibilities,
                     itsPredictBuffer->GetModel(thread).data(), tdata,
                     std::plus<dcomplex>());
    }
  }

  // Call ApplyCal step
  if (itsDoApplyCal) {
    itsApplyCalStep.process(tempBuffer);
    tempBuffer = itsResultStep->get();
    tdata = tempBuffer.getData().data();
  }

  // Put predict result from temp buffer into the 'real' buffer
  if (itsOperation == "replace") {
    itsBuffer = tempBuffer;
  } else {
    itsBuffer.copy(bufin);
    casacore::Complex* data = itsBuffer.getData().data();
    if (itsOperation == "add") {
      std::transform(data, data + nVisibilities, tdata, data,
                     std::plus<dcomplex>());
    } else if (itsOperation == "subtract") {
      std::transform(data, data + nVisibilities, tdata, data,
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
  const casacore::Vector<double>& itrf = itrfDir.getValue().getValue();
  everybeam::vector3r_t vec;
  vec[0] = itrf[0];
  vec[1] = itrf[1];
  vec[2] = itrf[2];
  return vec;
}

void Predict::addBeamToData(base::Patch::ConstPtr patch, double time,
                            const everybeam::vector3r_t& refdir,
                            const everybeam::vector3r_t& tiledir, size_t thread,
                            size_t nBeamValues, dcomplex* data0,
                            bool stokesIOnly) {
  // Apply beam for a patch, add result to Model
  MDirection dir(MVDirection(patch->position()[0], patch->position()[1]),
                 MDirection::J2000);
  everybeam::vector3r_t srcdir = dir2Itrf(dir, itsMeasConverters[thread]);

  float* dummyweight = 0;

  if (stokesIOnly) {
    ApplyBeam::applyBeamStokesIArrayFactor(
        info(), time, data0, dummyweight, srcdir, refdir, tiledir,
        itsPredictBuffer->GetStationList(),
        itsPredictBuffer->GetScalarBeamValues(thread), itsUseChannelFreq, false,
        itsBeamMode, false);
  } else {
    ApplyBeam::applyBeam(info(), time, data0, dummyweight, srcdir, refdir,
                         tiledir, itsPredictBuffer->GetStationList(),
                         itsPredictBuffer->GetFullBeamValues(thread),
                         itsUseChannelFreq, false, itsBeamMode, false);
  }

  // Add temporary buffer to Model
  std::transform(itsPredictBuffer->GetModel(thread).data(),
                 itsPredictBuffer->GetModel(thread).data() + nBeamValues, data0,
                 itsPredictBuffer->GetModel(thread).data(),
                 std::plus<dcomplex>());
}

void Predict::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
