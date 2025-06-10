// Demixer.cc: DP3 step class to subtract A-team sources
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "Demixer.h"

#include <iomanip>

#include <aocommon/dynamicfor.h>
#include <aocommon/logger.h>
#include <aocommon/xt/utensor.h>

#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include <dp3/base/DP3.h>

#include "../base/Apply.h"
#include "../base/CursorUtilCasa.h"
#include "../base/EstimateMixed.h"
#include "../base/SubtractMixed.h"
#include "../base/Simulate.h"
#include "../base/Simulator.h"

#include "../model/SkyModelCache.h"
#include "../model/SourceDBUtil.h"

#include "../parmdb/Axis.h"
#include "../parmdb/SourceDB.h"
#include "../parmdb/ParmDB.h"
#include "../parmdb/ParmSet.h"
#include "../parmdb/ParmCache.h"
#include "../parmdb/Parm.h"

#include "../common/StreamUtil.h"

#include <aocommon/logger.h>

#include "Averager.h"
#include "MsReader.h"
#include "NullStep.h"
#include "PhaseShift.h"

using casacore::IPosition;
using casacore::Matrix;
using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVEpoch;
using casacore::Quantum;

using aocommon::Logger;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

using dp3::common::operator<<;

namespace {
std::string toString(double value);
}  // end unnamed namespace

Demixer::Demixer(const common::ParameterSet& parset, const std::string& prefix)
    : itsName(prefix),
      itsSkyName(parset.getString(prefix + "skymodel", "sky")),
      itsInstrumentName(
          parset.getString(prefix + "instrumentmodel", "instrument")),
      itsDefaultGain(parset.getDouble(prefix + "defaultgain", 1.0)),
      itsMaxIter(parset.getInt(prefix + "maxiter", 50)),
      itsSelBL(parset, prefix, false, "cross"),
      itsAvgResultSubtr(nullptr),
      itsIgnoreTarget(parset.getBool(prefix + "ignoretarget", false)),
      itsTargetSource(parset.getString(prefix + "targetsource", std::string())),
      itsSubtrSources(parset.getStringVector(prefix + "subtractsources",
                                             std::vector<std::string>())),
      itsModelSources(parset.getStringVector(prefix + "modelsources",
                                             std::vector<std::string>())),
      itsExtraSources(parset.getStringVector(prefix + "othersources",
                                             std::vector<std::string>())),
      itsPropagateSolutions(
          parset.getBool(prefix + "propagatesolutions", false)),
      itsNBl(0),
      itsNChanAvgSubtr(parset.getUint(prefix + "freqstep", 1)),
      itsNChanIn(0),
      itsNChanOut(0),
      itsNChanOutSubtr(0),
      itsNCorr(0),
      itsNDir(0),
      itsNModel(0),
      itsNStation(0),
      itsNTimeAvgSubtr(parset.getUint(prefix + "timestep", 1)),
      itsNTimeChunk(parset.getUint(prefix + "ntimechunk", 0)),
      itsNTimeChunkSubtr(0),
      itsNTimeDemix(0),
      itsNTimeIn(0),
      itsNTimeOut(0),
      itsNTimeOutSubtr(0),
      itsTimeIntervalAvg(0),
      itsUseLBFGS(parset.getBool(prefix + "uselbfgssolver", false)),
      itsLBFGShistory(parset.getUint(prefix + "lbfgs.historysize", 10)),
      itsLBFGSrobustdof(parset.getDouble(prefix + "lbfgs.robustdof", 2.0)),
      itsRangeLBFGSsol(parset.getDoubleVector(prefix + "lbfgs.solution.range",
                                              std::vector<double>())),
      itsTimeIndex(0),
      itsNConverged(0) {
  if (itsSkyName.empty() || itsInstrumentName.empty())
    throw std::runtime_error(
        "An empty name is given for the sky and/or instrument model");
  if (itsIgnoreTarget && !itsTargetSource.empty())
    throw std::runtime_error(
        "Target source name cannot be given if ignoretarget=true");

  // Try parsing for demix{time,freq} resolution keys first,
  // if not demix{time,freq}step parset keys
  itsFreqResolution = Averager::getFreqHz(
      parset.getString(prefix + "demixfreqresolution", "0"));
  if (itsFreqResolution <= 0) {
    itsNChanAvg = parset.getUint(prefix + "demixfreqstep", itsNChanAvgSubtr);
  }
  itsTimeResolution = parset.getFloat(prefix + "demixtimeresolution", 0.);
  if (itsTimeResolution <= 0) {
    itsNTimeAvg = parset.getUint(prefix + "demixtimestep", itsNTimeAvgSubtr);
  }
  if ((itsFreqResolution > 0 && itsTimeResolution == 0) ||
      (itsFreqResolution == 0 && itsTimeResolution > 0))
    throw std::runtime_error(
        "Both time and frequency resolutions should be given");
  const bool use_resolution = (itsFreqResolution > 0 && itsTimeResolution > 0);

  if (itsFreqResolution > 0) {
    itsNChanAvg = 1;  // will be updated in updateInfo()
  }
  if (itsTimeResolution > 0) {
    itsNTimeAvg = 1;  // will be updated in updateInfo()
  }

#ifndef HAVE_LIBDIRAC
  if (itsUseLBFGS)
    throw std::runtime_error(
        "uselbfgssolver=true but libdirac is not available");
#endif /* ! HAVE_LIBDIRAC */

#ifdef HAVE_LIBDIRAC
  if (itsUseLBFGS) {
    if (itsRangeLBFGSsol.size() == 0) {
      // set up default range to [0,0]
      itsRangeLBFGSsol.push_back(0.);
      itsRangeLBFGSsol.push_back(0.);
    } else if (itsRangeLBFGSsol.size() == 2) {
      // user has already specified a range, check its validity
      if (itsRangeLBFGSsol[0] >= itsRangeLBFGSsol[1]) {
        throw std::runtime_error(
            "Invalid range for lbfgs.solution.range=[sol_min,sol_max], sol_min "
            "< sol_max");
      }
    } else {
      throw std::runtime_error(
          "Invalid range for lbfgs.solution.range=[sol_min,sol_max], sol_min < "
          "sol_max");
    }
  }
#endif /* HAVE_LIBDIRAC */

  // Add a result step as last step in the filter.
  itsFilter = std::make_shared<Filter>(itsSelBL);
  itsFilterResult = std::make_shared<ResultStep>();
  itsFilter->setNextStep(itsFilterResult);
  // Default nr of time chunks is maximum number of threads.
  if (itsNTimeChunk == 0) {
    itsNTimeChunk = aocommon::ThreadPool::GetInstance().NThreads();
  }
  // Check that time windows fit integrally.
  // This check will be done later when using time/freq resolution
  if (!use_resolution && (itsNTimeChunk * itsNTimeAvg) % itsNTimeAvgSubtr != 0)
    throw std::runtime_error(
        "time window should fit final averaging integrally");
  itsNTimeChunkSubtr = (itsNTimeChunk * itsNTimeAvg) / itsNTimeAvgSubtr;
  // Collect all source names.
  itsNModel = itsSubtrSources.size() + itsModelSources.size();
  itsNDir = itsNModel + itsExtraSources.size() + 1;
  itsAllSources.reserve(itsNDir);
  itsAllSources.insert(itsAllSources.end(), itsSubtrSources.begin(),
                       itsSubtrSources.end());
  itsAllSources.insert(itsAllSources.end(), itsModelSources.begin(),
                       itsModelSources.end());
  itsAllSources.insert(itsAllSources.end(), itsExtraSources.begin(),
                       itsExtraSources.end());
  if (!itsTargetSource.empty()) {
    itsAllSources.push_back(itsTargetSource);
  }

  // Get the source info of all patches from the SourceDB table.
  std::vector<std::string> patchNames(itsAllSources);
  // If the target source is given, add it to the model.
  // Because the target source has to be the last direction, it means
  // that (for the time being) no extra sources can be given.
  if (!itsTargetSource.empty()) {
    patchNames[itsNModel++] = itsTargetSource;
    // The target has to be the last demix direction.
    // If it has a source model, there cannot be any extra source
    // because the sources to be predicted have to be a consecutive vector.
    if (!itsExtraSources.empty())
      throw std::runtime_error(
          "Currently no extrasources can "
          "be given if the targetsource is given");
  }
  model::SkyModelCache& cache = model::SkyModelCache::GetInstance();
  itsPatchList =
      cache.GetSkyModel(itsSkyName)
          .Filter(patchNames, model::SourceDBWrapper::FilterMode::kValue)
          .MakePatchList();
  assert(itsPatchList.size() == itsNModel);

  // Size buffers.
  itsFactors.resize(itsNTimeChunk);
  itsFactorsSubtr.resize(itsNTimeChunkSubtr);
  itsPhaseShifts.reserve(itsNDir - 1);  // not needed for target direction
  itsFirstSteps.reserve(itsNDir);
  itsAvgResults.reserve(itsNDir);

  // Create the solve and subtract steps for the sources to be removed.
  // Solving consists of the following steps:
  // - select the requested baselines (longer baselines may need no demix)
  // - phaseshift selected data to each demix source
  // - average selected data in each direction, also original phasecenter.
  // - determine and average demix factors for all directions
  // - predict and solve in each direction. It is possible to predict
  //   more directions than to solve (for strong sources in field).
  // Subtract consists of the following steps:
  // - average all data (possibly different averaging than used in solve)
  // - determine and average demix factors (using select output in solve)
  // - select the requested baselines
  // - subtract sources for selected data
  // - merge subtract result into averaged data. This is not needed if
  //   no selection is done.
  // Note that multiple time chunks are handled jointly, so a
  // MultiResultStep is used to catch the results of all time chunks.
  for (unsigned int i = 0; i < itsNDir - 1; ++i) {
    // First make the phaseshift and average steps for each demix source.
    // The resultstep gets the result.
    // The phasecenter can be given in a parameter. Its name is the default.
    // Look up the source direction in the patch table.
    // If found, turn it into a vector of strings.
    std::vector<std::string> sourceVec(1, itsAllSources[i]);
    if (i < itsNModel) {
      sourceVec[0] = toString(itsPatchList[i]->Direction().ra);
      sourceVec.push_back(toString(itsPatchList[i]->Direction().dec));
    }
    auto step1 = std::make_shared<PhaseShift>(
        parset, prefix + itsAllSources[i] + '.', sourceVec);
    itsFirstSteps.push_back(step1);
    itsPhaseShifts.push_back(step1);
    auto step2 =
        (use_resolution
             ? std::make_shared<Averager>(prefix, itsFreqResolution,
                                          itsTimeResolution)
             : std::make_shared<Averager>(prefix, itsNChanAvg, itsNTimeAvg));
    step1->setNextStep(step2);
    auto step3 = std::make_shared<MultiResultStep>(itsNTimeChunk);
    step2->setNextStep(step3);
    // There is a single demix factor step which needs to get all results.
    itsAvgResults.push_back(step3);
  }

  // Now create the step to average the data themselves.
  auto targetAvg =
      (use_resolution
           ? std::make_shared<Averager>(prefix, itsFreqResolution,
                                        itsTimeResolution)
           : std::make_shared<Averager>(prefix, itsNChanAvg, itsNTimeAvg));
  itsFirstSteps.push_back(targetAvg);
  auto targetAvgRes = std::make_shared<MultiResultStep>(itsNTimeChunk);
  targetAvg->setNextStep(targetAvgRes);
  itsAvgResults.push_back(targetAvgRes);

  // Create the data average step for the subtract.
  // The entire average result is needed for the next NDPPP step.
  // Only the selected baselines need to be subtracted, so add a
  // filter step as the last one.
  itsAvgStepSubtr =
      std::make_shared<Averager>(prefix, itsNChanAvgSubtr, itsNTimeAvgSubtr);
  itsAvgResultFull = std::make_shared<MultiResultStep>(itsNTimeChunkSubtr);
  itsFilterSubtr = std::make_shared<Filter>(itsSelBL);
  itsAvgResultSubtr = std::make_shared<MultiResultStep>(itsNTimeChunkSubtr);
  itsAvgStepSubtr->setNextStep(itsAvgResultFull);
  itsAvgResultFull->setNextStep(itsFilterSubtr);
  itsFilterSubtr->setNextStep(itsAvgResultSubtr);
}

common::Fields Demixer::getRequiredFields() const {
  // Demixer always runs at least one Averager. There are no substeps with
  // additional requirements.
  return Averager::kRequiredFields;
}

common::Fields Demixer::getProvidedFields() const {
  // Demixer always runs at least one Averager. There are no substeps that
  // provide additional fields.
  return Averager::kProvidedFields;
}

void Demixer::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  // Get size info.
  itsNChanIn = infoIn.nchan();
  itsNCorr = infoIn.ncorr();
  if (itsNCorr != 4)
    throw std::runtime_error("Demixing requires data with 4 polarizations");

  // Handle possible data selection.
  itsFilter->setInfo(infoIn);
  const DPInfo& infoSel = itsFilter->getInfo();
  // NB. The number of baselines and stations refer to the number of
  // selected baselines and the number of unique stations that participate
  // in the selected baselines.
  itsNBl = infoSel.nbaselines();
  itsNStation = infoSel.antennaUsed().size();

  // Re-number the station IDs in the selected baselines, removing gaps in
  // the numbering due to unused stations.
  const std::vector<int>& antennaMap = infoSel.antennaMap();
  for (unsigned int i = 0; i < itsNBl; ++i) {
    itsBaselines.push_back(base::Baseline(antennaMap[infoSel.getAnt1()[i]],
                                          antennaMap[infoSel.getAnt2()[i]]));
  }

  // Prepare conversion from relative to absolute UVW
  std::vector<int> newAnt1(itsNBl);
  std::vector<int> newAnt2(itsNBl);
  for (unsigned int i = 0; i < itsNBl; ++i) {
    newAnt1[i] = antennaMap[infoSel.getAnt1()[i]];
    newAnt2[i] = antennaMap[infoSel.getAnt2()[i]];
  }
  itsUVWSplitIndex = base::nsetupSplitUVW(itsNStation, newAnt1, newAnt2);

  // Allocate buffers used to compute the smearing factors.
  itsFactorBuf.resize(
      {(itsNDir * (itsNDir - 1) / 2), itsNBl, itsNChanIn, itsNCorr});
  itsFactorBufSubtr.resize(
      {(itsNDir * (itsNDir - 1) / 2), itsNBl, itsNChanIn, itsNCorr});
  itsFactorBuf.fill(std::complex<double>(0.0, 0.0));
  itsFactorBufSubtr.fill(std::complex<double>(0.0, 0.0));

  // Adapt averaging to available nr of channels and times.
  // Use a copy of the DPInfo, otherwise it is updated multiple times.
  DPInfo infoDemix(infoSel);
  if (itsTimeResolution > 0) {
    double time_interval = infoDemix.timeInterval();
    itsNTimeAvg = std::max(1, (int)(itsTimeResolution / time_interval + 0.5));
    if ((itsNTimeChunk * itsNTimeAvg) % itsNTimeAvgSubtr != 0)
      throw std::runtime_error(
          "time window should fit final averaging integrally");
  }

  itsNTimeAvg = std::min(itsNTimeAvg, infoSel.ntime());
  if (itsFreqResolution > 0) {
    double chan_width = infoDemix.chanWidths()[0];
    itsNChanAvg = std::max(1, (int)(itsFreqResolution / chan_width + 0.5));
  }
  itsNChanAvg = infoDemix.update(itsNChanAvg, itsNTimeAvg);
  itsNChanOut = infoDemix.nchan();
  itsTimeIntervalAvg = infoDemix.timeInterval();
  itsNTimeDemix = infoDemix.ntime();

  unsigned updated_itsNTimeChunkSubtr =
      (itsNTimeChunk * itsNTimeAvg) / itsNTimeAvgSubtr;

  if (updated_itsNTimeChunkSubtr != itsNTimeChunkSubtr) {
    itsNTimeChunkSubtr = updated_itsNTimeChunkSubtr;
    itsAvgResultFull = std::make_shared<MultiResultStep>(itsNTimeChunkSubtr);
    itsAvgResultSubtr = std::make_shared<MultiResultStep>(itsNTimeChunkSubtr);
    itsAvgStepSubtr->setNextStep(itsAvgResultFull);
    itsAvgResultFull->setNextStep(itsFilterSubtr);
    itsFilterSubtr->setNextStep(itsAvgResultSubtr);

    itsFactorsSubtr.resize(itsNTimeChunkSubtr);
  }

  itsNTimeDemix = infoDemix.ntime();

  // Let the internal steps update their data.
  for (unsigned int i = 0; i < itsFirstSteps.size(); ++i) {
    itsFirstSteps[i]->setInfo(infoSel);
  }
  itsAvgStepSubtr->setInfo(infoIn);
  // Update the info of this object.
  itsNTimeAvgSubtr = std::min(itsNTimeAvgSubtr, infoSel.ntime());
  itsNChanAvgSubtr = info().update(itsNChanAvgSubtr, itsNTimeAvgSubtr);
  itsNChanOutSubtr = info().nchan();
  if (itsNChanAvg % itsNChanAvgSubtr != 0)
    throw std::runtime_error("Demix averaging " + std::to_string(itsNChanAvg) +
                             " must be multiple of output averaging " +
                             std::to_string(itsNChanAvgSubtr));
  if (itsNTimeAvg % itsNTimeAvgSubtr != 0)
    throw std::runtime_error("Demix averaging " + std::to_string(itsNTimeAvg) +
                             " must be multiple of output averaging " +
                             std::to_string(itsNTimeAvgSubtr));
  // Store channel frequencies for the demix and subtract resolutions.
  itsFreqDemix = infoDemix.chanFreqs();
  itsFreqSubtr = getInfo().chanFreqs();

  // Store phase center direction in J2000.
  try {
    MDirection dirJ2000(
        MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
    Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
    itsPhaseRef =
        base::Direction(angles.getBaseValue()[0], angles.getBaseValue()[1]);
    itsMovingPhaseRef = false;
  } catch (casacore::AipsError&) {
    // Phase direction (in J2000) is time dependent
    itsMovingPhaseRef = true;
    Logger::Warn
        << "WARNING: Demixing with moving phase reference is not tested.\n";
  }

  // Initialize the unknowns.
  itsUnknowns.resize(itsNTimeDemix * itsNModel * itsNStation * 8);
  itsPrevSolution.resize(itsNModel * itsNStation * 8);
  std::vector<double>::iterator it = itsPrevSolution.begin();
  std::vector<double>::iterator it_end = itsPrevSolution.end();
  while (it != it_end) {
    *it++ = itsDefaultGain;
    *it++ = 0.0;
    *it++ = 0.0;
    *it++ = 0.0;
    *it++ = 0.0;
    *it++ = 0.0;
    *it++ = itsDefaultGain;
    *it++ = 0.0;
  }
  // Initialize the flag counters.
  itsFlagCounter.init(getInfo());
}

void Demixer::show(std::ostream& os) const {
  os << "Demixer " << itsName << '\n';
  os << "  skymodel:           " << itsSkyName << '\n';
  os << "  instrumentmodel:    " << itsInstrumentName << '\n';
  os << "  default gain:       " << itsDefaultGain << '\n';
  os << "  max iterations:     " << itsMaxIter << '\n';
  itsSelBL.show(os);
  if (itsSelBL.hasSelection()) {
    os << "    demixing " << itsFilter->getInfo().nbaselines() << " out of "
       << getInfo().nbaselines() << " baselines   ("
       << itsFilter->getInfo().antennaUsed().size() << " out of "
       << getInfo().antennaUsed().size() << " stations)" << '\n';
  }
  os << "  targetsource:       " << itsTargetSource << '\n';
  os << "  subtractsources:    " << itsSubtrSources << '\n';
  unsigned int inx = 0;
  for (unsigned int i = 0; i < itsSubtrSources.size(); ++i) {
    os << "                        " << itsPhaseShifts[inx++]->getPhaseCenter()
       << '\n';
  }
  os << "  modelsources:       " << itsModelSources << '\n';
  for (unsigned int i = 0; i < itsModelSources.size(); ++i) {
    os << "                        " << itsPhaseShifts[inx++]->getPhaseCenter()
       << '\n';
  }
  os << "  extrasources:       " << itsExtraSources << '\n';
  for (unsigned int i = 0; i < itsExtraSources.size(); ++i) {
    os << "                        " << itsPhaseShifts[inx++]->getPhaseCenter()
       << '\n';
  }
  os << "  propagatesolutions: " << std::boolalpha << itsPropagateSolutions
     << std::noboolalpha << '\n';
  os << "  freqstep:           " << itsNChanAvgSubtr << '\n';
  os << "  timestep:           " << itsNTimeAvgSubtr << '\n';
  os << "  demixfreqstep:      " << itsNChanAvg << '\n';
  os << "  demixtimestep:      " << itsNTimeAvg << '\n';
  os << "  demixfreqresolution (Hz):      " << itsFreqResolution << '\n';
  os << "  demixtimeresolution (s):      " << itsTimeResolution << '\n';
  os << "  ntimechunk:         " << itsNTimeChunk << '\n';
}

void Demixer::showCounts(std::ostream& os) const {
  os << '\n' << "Statistics for Demixer " << itsName;
  os << '\n' << "======================" << '\n';
  os << '\n'
     << "Converged: " << itsNConverged << "/" << itsNTimeDemix << " cells"
     << '\n';
}

void Demixer::showTimings(std::ostream& os, double duration) const {
  const double self = itsTimer.getElapsed();

  os << "  ";
  FlagCounter::showPerc1(os, self, duration);
  os << " Demixer " << itsName << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerPhaseShift.getElapsed(), self);
  os << " of it spent in phase shifting/averaging data" << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, itsTimerDemix.getElapsed(), self);
  os << " of it spent in calculating decorrelation factors" << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, itsTimerSolve.getElapsed(), self);
  os << " of it spent in estimating gains and computing residuals" << '\n';
  os << "          ";
  FlagCounter::showPerc1(os, itsTimerDump.getElapsed(), self);
  os << " of it spent in writing gain solutions to disk" << '\n';
}

bool Demixer::process(std::unique_ptr<DPBuffer> buffer) {
  itsTimer.start();
  // Update the count.
  itsNTimeIn++;

  // Do the filter step first.
  common::Fields subChainReqFields = base::GetChainRequiredFields(itsFilter);
  itsFilter->process(std::make_unique<DPBuffer>(*buffer, subChainReqFields));
  std::unique_ptr<DPBuffer> selection_buffer = itsFilterResult->take();

  // Do the next steps (phaseshift and average) on the filter output.
  itsTimerPhaseShift.start();
  subChainReqFields = base::GetChainRequiredFields(itsFirstSteps[0]);
  for (int i = 0; i < int(itsFirstSteps.size()); ++i) {
    itsFirstSteps[i]->process(
        std::make_unique<DPBuffer>(*selection_buffer, subChainReqFields));
  }
  // Do the average and filter step for the output for all data.
  itsAvgStepSubtr->process(std::move(buffer));
  itsTimerPhaseShift.stop();

  // Per direction pair, calculate the phase rotation factors for the selected
  // data and accumulate them over time. Since the solving and subtract parts
  // apply different averaging parameters, the results are stored separately
  // in itsFactorBuf and itsFactorBufSubtr, respectively.
  itsTimerDemix.start();
  addFactors(std::move(selection_buffer));
  // The solving part
  if (itsNTimeIn % itsNTimeAvg == 0) {
    makeFactors(itsFactorBuf, itsFactors[itsNTimeOut],
                itsAvgResults[0]->get()[itsNTimeOut]->GetWeights(), itsNChanOut,
                itsNChanAvg);
    // Deproject sources without a model.
    deproject(itsFactors[itsNTimeOut], itsNTimeOut);
    itsFactorBuf.fill(std::complex<double>(0.0, 0.0));  // Clear buffer
    itsNTimeOut++;
  }
  // The subtract part
  if (itsNTimeIn % itsNTimeAvgSubtr == 0) {
    makeFactors(itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                itsAvgResultSubtr->get()[itsNTimeOutSubtr]->GetWeights(),
                itsNChanOutSubtr, itsNChanAvgSubtr);
    itsFactorBufSubtr.fill(std::complex<double>(0.0, 0.0));  // Clear buffer
    itsNTimeOutSubtr++;
  }
  itsTimerDemix.stop();

  // Estimate gains and subtract source contributions when sufficient time
  // slots have been collected.
  if (itsNTimeOut >= itsNTimeChunk) {
    handleDemix();
  }
  itsTimer.stop();
  return true;
}

void Demixer::finish() {
  Logger::Info << "  " << itsNTimeIn << " time slots to finish in Demixer ..."
               << '\n';
  itsTimer.start();

  // Process remaining entries.
  if (itsNTimeIn > 0) {
    // Finish the initial steps (phaseshift and average).
    itsTimerPhaseShift.start();
    for (int i = 0; i < int(itsFirstSteps.size()); ++i) {
      itsFirstSteps[i]->finish();
    }
    itsAvgStepSubtr->finish();
    itsTimerPhaseShift.stop();
    // Only average if there is some unaveraged data.
    itsTimerDemix.start();
    if (itsNTimeIn % itsNTimeAvg != 0) {
      makeFactors(itsFactorBuf, itsFactors[itsNTimeOut],
                  itsAvgResults[0]->get()[itsNTimeOut]->GetWeights(),
                  itsNChanOut, itsNChanAvg);
      // Deproject sources without a model.
      deproject(itsFactors[itsNTimeOut], itsNTimeOut);
      itsNTimeOut++;
    }
    if (itsNTimeIn % itsNTimeAvgSubtr != 0) {
      makeFactors(itsFactorBufSubtr, itsFactorsSubtr[itsNTimeOutSubtr],
                  itsAvgResultSubtr->get()[itsNTimeOutSubtr]->GetWeights(),
                  itsNChanOutSubtr, itsNChanAvgSubtr);
      itsNTimeOutSubtr++;
    }
    itsTimerDemix.stop();
    // Resize lists of mixing factors to the number of valid entries.
    itsFactors.resize(itsNTimeOut);
    itsFactorsSubtr.resize(itsNTimeOutSubtr);

    // Demix the source directions.
    handleDemix();
  }

  // Write solutions to disk in ParmDB format.
  itsTimerDump.start();
  dumpSolutions();
  itsTimerDump.stop();

  // Normalize variance ratio over total runs
  itsVarianceRatio /= float(itsTotalDemixRuns);

  itsTimer.stop();

  // Let the next steps finish.
  getNextStep()->finish();
}

void Demixer::handleDemix() {
  if (itsNModel > 0) {
    itsTimerSolve.start();
    demix();
    itsTimerSolve.stop();
    // If selection was done, merge the subtract results back into the
    // buffer.
  }
  // If needed, merge in the deselected baselines.
  if (itsSelBL.hasSelection()) {
    mergeSubtractResult();
  }

  // Clear the input buffers (MultiResultSteps).
  for (size_t i = 0; i < itsAvgResults.size(); ++i) {
    itsAvgResults[i]->clear();
  }
  // Let the next step process the data.
  for (unsigned int i = 0; i < itsNTimeOutSubtr; ++i) {
    itsTimer.stop();
    std::unique_ptr<DPBuffer> buffer_out;
    if (itsSelBL.hasSelection()) {
      buffer_out = std::move(itsAvgResultFull->get()[i]);
    } else {
      buffer_out = std::move(itsAvgResultSubtr->get()[i]);
    }

    MsReader::FlagInfinityNan(*buffer_out, itsFlagCounter);
    getNextStep()->process(std::move(buffer_out));
    itsTimer.start();
  }

  // Clear the output buffer (MultiResultSteps).
  itsAvgResultFull->clear();
  itsAvgResultSubtr->clear();

  // Reset counters.
  itsNTimeIn = 0;
  itsNTimeOut = 0;
  itsNTimeOutSubtr = 0;
  itsTimeIndex += itsNTimeChunk;
  itsTotalDemixRuns += 1;
}

void Demixer::mergeSubtractResult() {
  // Merge the selected baselines from the subtract buffer (itsAvgResultSubtr)
  // into the full buffer (itsAvgResultFull). Do it for all timestamps.
  for (unsigned int i = 0; i < itsNTimeOutSubtr; ++i) {
    const DPBuffer::DataType& arr = itsAvgResultSubtr->get()[i]->GetData();
    for (size_t baseline_in = 0; baseline_in < itsFilter->getIndicesBL().size();
         ++baseline_in) {
      const size_t baseline_out = itsFilter->getIndicesBL()[baseline_in];
      xt::view(itsAvgResultFull->get()[i]->GetData(), baseline_out, xt::all(),
               xt::all()) = xt::view(arr, baseline_in, xt::all(), xt::all());
    }
  }
}

void Demixer::addFactors(std::unique_ptr<DPBuffer> newBuf) {
  if (itsNDir <= 1) return;  // Nothing to do if only target direction.

  const size_t nbl = newBuf->GetData().shape(0);
  const size_t nchan = newBuf->GetData().shape(1);
  const size_t ncorr = newBuf->GetData().shape(2);

  const DPBuffer::FlagsType& flags = newBuf->GetFlags();
  const DPBuffer::WeightsType& weights = newBuf->GetWeights();

  // If ever in the future a time dependent phase center is used,
  // the machine must be reset for each new time, thus each new call
  // to process.
  // Add the weighted factors for each pair of directions.
  // The input factor is the phaseshift from target direction to
  // source direction. By combining them you get the shift from one
  // source direction to another.
  int dirnr = 0;  // direction pair number
  aocommon::StaticFor<size_t> loop;
  for (unsigned int dir0 = 0; dir0 < itsNDir - 1; ++dir0) {
    const xt::xtensor<std::complex<double>, 2>& phasors0 =
        itsPhaseShifts[dir0]->getPhasors();
    for (unsigned int dir1 = dir0 + 1; dir1 < itsNDir; ++dir1) {
      if (dir1 == itsNDir - 1) {
        // The last direction is the target direction, so no need to
        // combine the factors. Take conj to get shift source to target.
        loop.Run(0, nbl, [&](size_t start_baseline, size_t end_baseline) {
          for (size_t bl = start_baseline; bl < end_baseline; ++bl) {
            for (size_t chan = 0; chan < nchan; ++chan) {
              const std::complex<double> factor = std::conj(phasors0(bl, chan));
              for (size_t corr = 0; corr < ncorr; ++corr) {
                const std::complex<double> weighted_factor =
                    factor *
                    double(!flags(bl, chan, corr) * weights(bl, chan, corr));
                itsFactorBuf(dirnr, bl, chan, corr) += weighted_factor;
                itsFactorBufSubtr(dirnr, bl, chan, corr) += weighted_factor;
              }
            }
          }
        });  // end parallel for
      } else {
        // Different source directions; take both phase terms into account.
        const xt::xtensor<std::complex<double>, 2>& phasors1 =
            itsPhaseShifts[dir1]->getPhasors();

        loop.Run(0, nbl, [&](size_t start_baseline, size_t end_baseline) {
          for (size_t bl = start_baseline; bl < end_baseline; ++bl) {
            for (size_t chan = 0; chan < nchan; ++chan) {
              const std::complex<double> factor =
                  phasors1(bl, chan) * std::conj(phasors0(bl, chan));
              for (size_t corr = 0; corr < ncorr; ++corr) {
                const std::complex<double> weighted_factor =
                    factor *
                    double(!flags(bl, chan, corr) * weights(bl, chan, corr));
                itsFactorBuf(dirnr, bl, chan, corr) += weighted_factor;
                itsFactorBufSubtr(dirnr, bl, chan, corr) += weighted_factor;
              }
            }
          }
        });  // end parallel for
      }

      // Next direction pair.
      ++dirnr;
    }
  }
}

void Demixer::makeFactors(
    const aocommon::xt::UTensor<std::complex<double>, 4>& bufIn,
    aocommon::xt::UTensor<std::complex<double>, 5>& bufOut,
    const DPBuffer::WeightsType& weightSums, size_t nChanOut, size_t nChanAvg) {
  if (itsNDir <= 1) return;  // Nothing to do if only target direction.
  assert(weightSums.size() != 0);

  const size_t kMaxNrCorrelations = 4;
  bufOut.resize({itsNBl, nChanOut, itsNCorr, itsNDir, itsNDir});
  bufOut.fill(std::complex<double>(1.0, 0.0));
  // Fill the factors for each combination of different directions.
  unsigned int dirnr = 0;  // direction pair number of the input buffer
  aocommon::StaticFor<size_t> loop;
  for (unsigned int dir0 = 0; dir0 < itsNDir - 1; ++dir0) {
    for (unsigned int dir1 = dir0 + 1; dir1 < itsNDir; ++dir1) {
      // Average factors by summing channels.
      // Note that summing in time is done in addFactors.
      // The sum per output channel is divided by the summed weight.
      // Note there is a summed weight per baseline,outchan,corr.
      loop.Run(0, itsNBl, [&](size_t start_baseline, size_t end_baseline) {
        for (size_t bl = start_baseline; bl < end_baseline; ++bl) {
          size_t ch_in = 0;
          for (size_t ch_out = 0; ch_out < nChanOut; ++ch_out) {
            // Sum the factors for the input channels to average.
            std::array<std::complex<double>, kMaxNrCorrelations> sum = {0.0};
            // In theory, the last output channel could consist of fewer
            // input channels, so take care of that.
            const size_t n_input_channels =
                std::min(nChanAvg, itsNChanIn - ch_out * nChanAvg);
            const size_t last_input_channel = ch_in + n_input_channels;
            for (; ch_in < last_input_channel; ++ch_in) {
              for (unsigned int corr = 0; corr < itsNCorr; ++corr) {
                sum[corr] += bufIn(dirnr, bl, ch_in, corr);
              }
            }
            for (unsigned int corr = 0; corr < itsNCorr; ++corr) {
              const std::complex<double> factor_out =
                  sum[corr] / double(weightSums(bl, ch_out, corr));
              bufOut(bl, ch_out, corr, dir0, dir1) = factor_out;
              bufOut(bl, ch_out, corr, dir1, dir0) = std::conj(factor_out);
            }
          }
        }
      });  // end parallel for
      // Next input direction pair.
      ++dirnr;
    }
  }
}

void Demixer::deproject(aocommon::xt::UTensor<std::complex<double>, 5>& factors,
                        unsigned int resultIndex) {
  // Sources without a model have to be deprojected.
  // Optionally no deprojection of target direction.
  unsigned int nrDeproject = itsNDir - itsNModel;
  if (itsIgnoreTarget) {
    nrDeproject--;
  }
  // Nothing to do if only target direction or nothing to deproject.
  if (itsNDir <= 1 || nrDeproject == 0) return;
  // Get pointers to the data for the various directions.
  std::vector<std::complex<float>*> resultPtr(itsNDir);
  for (unsigned int j = 0; j < itsNDir; ++j) {
    resultPtr[j] = itsAvgResults[j]->get()[resultIndex]->GetData().data();
  }
  // The projection matrix is given by
  //     P = I - A * inv(A.T.conj * A) * A.T.conj
  // where A is the last column of the demixing matrix M.
  // The BBS equations get:
  //     P * M' * v_predict = P * v_averaged
  // where M' is obtained by removing the last column of demixing matrix M.
  // The dimensions of the matrices/vectors are:
  //     P : NxN
  //     M' : Nx(N-1)
  //     v_predict : (N-1) x 1
  //     v_averaged: N x 1
  // where N is the number of modeled sources to use in demixing.
  // In the general case, S sources might not have a source model.
  // In that case, A is the NxS matrix containing all these columns
  // from M and M' is the Nx(N-S) matrix without all these columns.

  // Calculate P for all baselines,channels,correlations.
  std::array<size_t, 5> shape = factors.shape();
  int nvis = shape[0] * shape[1] * shape[2];  // bl * chan * corr
  shape[3] = itsNModel;
  aocommon::xt::UTensor<std::complex<double>, 5> newFactors(shape);
  IPosition inShape(2, itsNDir, itsNDir);
  IPosition outShape(2, itsNDir, itsNModel);
  {
    Matrix<casacore::DComplex> a(itsNDir, nrDeproject);
    Matrix<casacore::DComplex> ma(itsNDir, itsNModel);
    std::vector<std::complex<double>> vec(itsNDir);
    for (int i = 0; i < nvis; ++i) {
      // Split the matrix into the modeled and deprojected sources.
      // Copy the columns to the individual matrices.
      const std::complex<double>* inptr =
          factors.data() + i * itsNDir * itsNDir;
      std::complex<double>* outptr =
          newFactors.data() + i * itsNDir * itsNModel;
      Matrix<casacore::DComplex> out(outShape, outptr, casacore::SHARE);
      // Copying a bit of data is probably faster than taking a matrix subset.
      casacore::objcopy(ma.data(), inptr, itsNDir * itsNModel);
      casacore::objcopy(a.data(), inptr + itsNDir * itsNModel,
                        itsNDir * nrDeproject);
      // Calculate conjugated transpose of A, multiply with A, and invert.
      Matrix<casacore::DComplex> at(adjoint(a));
      Matrix<casacore::DComplex> ata(invert(product(at, a)));
      if (ata.empty()) {
        ata.resize(nrDeproject, nrDeproject);
      }
      assert(ata.ncolumn() == nrDeproject && ata.nrow() == nrDeproject);
      // Calculate P = I - A * ata * A.T.conj
      Matrix<casacore::DComplex> aata(product(a, ata));
      Matrix<casacore::DComplex> p(-product(product(a, ata), at));
      casacore::Vector<casacore::DComplex> diag(p.diagonal());
      diag += casacore::DComplex(1, 0);
      // Multiply the demixing factors with P (get stored in newFactors).
      out = product(p, ma);
      // Multiply the averaged data point with P.
      std::fill(vec.begin(), vec.end(), casacore::DComplex());
      for (unsigned int j = 0; j < itsNDir; ++j) {
        for (unsigned int k = 0; k < itsNDir; ++k) {
          vec[k] += casacore::DComplex(resultPtr[j][i]) * p(k, j);
        }
      }
      // Put result back in averaged data for those sources.
      for (unsigned int j = 0; j < itsNDir; ++j) {
        resultPtr[j][i] = vec[j];
      }
    }
  }
  // Set the new demixing factors.
  factors = std::move(newFactors);
}

namespace {
struct ThreadPrivateStorage {
  std::vector<double> unknowns;
  xt::xtensor<double, 2> uvw;
  std::vector<casacore::Cube<std::complex<double>>> model;
  casacore::Cube<std::complex<double>> model_subtr;
  size_t count_converged;
  float variance_before;
  float variance_after;
};

void initThreadPrivateStorage(ThreadPrivateStorage& storage, size_t nDirection,
                              size_t nStation, size_t nBaseline,
                              size_t nChannel, size_t nChannelSubtr) {
  storage.unknowns.resize(nDirection * nStation * 8);
  storage.uvw.resize({nStation, 3});
  storage.model.resize(nDirection);
  for (unsigned int dir = 0; dir < nDirection; ++dir) {
    storage.model[dir].resize(4, nChannel, nBaseline);
  }
  storage.model_subtr.resize(4, nChannelSubtr, nBaseline);
  storage.count_converged = 0;
  storage.variance_before = 0.0f;
  storage.variance_after = 0.0f;
}
}  // end unnamed namespace

void Demixer::demix() {
  const size_t nThread = aocommon::ThreadPool::GetInstance().NThreads();
  const size_t nTime = itsAvgResults[0]->size();
  const size_t nTimeSubtr = itsAvgResultSubtr->size();
  const size_t multiplier = itsNTimeAvg / itsNTimeAvgSubtr;
  const size_t nDr = itsNModel;
  const size_t nDrSubtr = itsSubtrSources.size();
  const size_t nSt = itsNStation;
  const size_t nBl = itsBaselines.size();
  const size_t nCh = itsFreqDemix.size();
  const size_t nChSubtr = itsFreqSubtr.size();
  const size_t nCr = 4;

  std::vector<ThreadPrivateStorage> threadStorage(nThread);
  for (std::vector<ThreadPrivateStorage>::iterator it = threadStorage.begin(),
                                                   end = threadStorage.end();
       it != end; ++it) {
    initThreadPrivateStorage(*it, nDr, nSt, nBl, nCh, nChSubtr);

    // Copy the previous solution to the thread private vectors of unknowns.
    // When solution propagation is disabled, itsPrevSolution is never
    // updated. It then contains 1.0+0.0i for the diagonal terms and
    // 0.0+0.0i for the off-diagonal terms. Thus, when solution propagation
    // is disabled this statement effectively re-initializes the thread
    // private vectors of unknowns.
    copy(itsPrevSolution.begin(), itsPrevSolution.end(), it->unknowns.begin());
  }

  base::const_cursor<base::Baseline> cr_baseline(&(itsBaselines[0]));

  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, nTime, [&](size_t ts, size_t thread) {
    ThreadPrivateStorage& storage = threadStorage[thread];

    // If solution propagation is disabled, re-initialize the thread-private
    // vector of unknowns.
    if (!itsPropagateSolutions) {
      std::copy(itsPrevSolution.begin(), itsPrevSolution.end(),
                storage.unknowns.begin());
    }

    // Simulate.
    //
    // Model visibilities for each direction of interest will be computed
    // and stored.
    size_t stride_model[3] = {1, nCr, nCr * nCh};
    std::fill(storage.model.begin(), storage.model.end(), 0.0);

    for (size_t dr = 0; dr < nDr; ++dr) {
      base::nsplitUVW(itsUVWSplitIndex, itsBaselines,
                      itsAvgResults[dr]->get()[ts]->GetUvw(), storage.uvw);

      base::Simulator simulator(itsPatchList[dr]->Direction(), nSt,
                                itsBaselines, itsFreqDemix, {}, storage.uvw,
                                storage.model[dr], false, false);
      for (size_t i = 0; i < itsPatchList[dr]->NComponents(); ++i) {
        simulator.simulate(itsPatchList[dr]->component(i));
      }
    }

    // Estimate Jones matrices.
    //
    // A Jones matrix will be estimated for each pair of station and
    // direction.
    //
    // A single (overdetermined) non-linear set of equations for all
    // stations and directions is solved iteratively. The influence of
    // each direction on each other direction is given by the mixing
    // matrix.
    // Legacy port the XTensors to Casacore structures until usage of
    // (const) cursors gets deprecated in the future.
    std::unique_ptr<DPBuffer>& source_buffer = itsAvgResults[0]->get()[ts];
    const IPosition cube_shape(3, source_buffer->GetData().shape(2),
                               source_buffer->GetData().shape(1),
                               source_buffer->GetData().shape(0));
    base::const_cursor<bool> cr_flag =
        base::casa_const_cursor(casacore::Cube<bool>(
            cube_shape, source_buffer->GetFlags().data(), casacore::SHARE));
    base::const_cursor<float> cr_weight =
        base::casa_const_cursor(casacore::Cube<float>(
            cube_shape, source_buffer->GetWeights().data(), casacore::SHARE));

    const IPosition shape(5, itsFactors[ts].shape(4), itsFactors[ts].shape(3),
                          itsFactors[ts].shape(2), itsFactors[ts].shape(1),
                          itsFactors[ts].shape(0));
    const casacore::Array<casacore::DComplex> demix_factor(
        shape, itsFactors[ts].data(), casacore::SHARE);
    base::const_cursor<std::complex<double>> cr_mix =
        base::casa_const_cursor(demix_factor);

    std::vector<base::const_cursor<std::complex<float>>> cr_data(nDr);
    std::vector<base::const_cursor<std::complex<double>>> cr_model(nDr);
    for (size_t dr = 0; dr < nDr; ++dr) {
      cr_data[dr] = base::casa_const_cursor(casacore::Cube<casacore::Complex>(
          cube_shape, itsAvgResults[dr]->get()[ts]->GetData().data(),
          casacore::SHARE));
      cr_model[dr] = base::const_cursor<std::complex<double>>(
          storage.model[dr].data(), 3, stride_model);
    }

    const bool converged =
#ifdef HAVE_LIBDIRAC
        (itsUseLBFGS
             ? estimate(nDr, nSt, nBl, nCh, cr_baseline, cr_data, cr_model,
                        cr_flag, cr_weight, cr_mix, &(storage.unknowns[0]),
                        itsLBFGShistory, itsLBFGSrobustdof, itsRangeLBFGSsol[0],
                        itsRangeLBFGSsol[1], itsMaxIter)
             : estimate(nDr, nSt, nBl, nCh, cr_baseline, cr_data, cr_model,
                        cr_flag, cr_weight, cr_mix, &(storage.unknowns[0]),
                        itsMaxIter));
#else
        estimate(nDr, nSt, nBl, nCh, cr_baseline, cr_data, cr_model, cr_flag,
                 cr_weight, cr_mix, &(storage.unknowns[0]), itsMaxIter);
#endif /* HAVE_LIBDIRAC */

    if (converged) {
      ++storage.count_converged;
    }

    // Compute the residual.
    //
    // All the so-called "subtract sources" are subtracted from the
    // observed data. The previously estimated Jones matrices, as well as
    // the appropriate mixing weight are applied before subtraction.
    //
    // Note that the resolution of the residual can differ from the
    // resolution at which the Jones matrices were estimated.
    for (size_t ts_subtr = multiplier * ts,
                ts_subtr_end = std::min(ts_subtr + multiplier, nTimeSubtr);
         ts_subtr != ts_subtr_end; ++ts_subtr) {
      for (size_t dr = 0; dr < nDrSubtr; ++dr) {
        // Re-use simulation used for estimating Jones matrices if possible.
        base::cursor<std::complex<double>> cr_model_subtr(
            storage.model[dr].data(), 3, stride_model);

        // Re-simulate if required.
        if (multiplier != 1 || nCh != nChSubtr) {
          base::nsplitUVW(itsUVWSplitIndex, itsBaselines,
                          itsAvgResultSubtr->get()[ts_subtr]->GetUvw(),
                          storage.uvw);

          if (itsMovingPhaseRef) {
            // Convert phase reference to J2000
            itsMeasFrame.set(
                MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
            MDirection dirJ2000(MDirection::Convert(
                info().phaseCenter(),
                MDirection::Ref(MDirection::J2000, itsMeasFrame))());
            Quantum<casacore::Vector<double>> angles = dirJ2000.getAngle();
            itsPhaseRef = base::Direction(angles.getBaseValue()[0],
                                          angles.getBaseValue()[1]);
          }
          // Rotate the UVW coordinates for the target direction to the
          // direction of source to subtract. This is required because at
          // the resolution of the residual the UVW coordinates for
          // directions other than the target are unavailable (unless the
          // resolution of the residual is equal to the resolution at which
          // the Jones matrices were estimated, of course).
          rotateUVW(itsPhaseRef, itsPatchList[dr]->Direction(), nSt,
                    storage.uvw.data());

          // Zero the visibility buffer.
          storage.model_subtr = std::complex<double>();

          // Simulate visibilities at the resolution of the residual.
          size_t stride_model_subtr[3] = {1, nCr, nCr * nChSubtr};
          cr_model_subtr = base::cursor<std::complex<double>>(
              storage.model_subtr.data(), 3, stride_model_subtr);

          base::Simulator simulator(itsPatchList[dr]->Direction(), nSt,
                                    itsBaselines, itsFreqSubtr, {}, storage.uvw,
                                    storage.model_subtr, false, false);
          for (size_t i = 0; i < itsPatchList[dr]->NComponents(); ++i) {
            simulator.simulate(itsPatchList[dr]->component(i));
          }
        }

        // Apply Jones matrices.
        size_t stride_unknowns[2] = {1, 8};
        base::const_cursor<double> cr_unknowns(
            &(storage.unknowns[dr * nSt * 8]), 2, stride_unknowns);

        apply(nBl, nChSubtr, cr_baseline, cr_unknowns, cr_model_subtr);

        // Subtract the source contribution from the data.
        const IPosition cube_shape(
            3, itsAvgResultSubtr->get()[ts_subtr]->GetData().shape(2),
            itsAvgResultSubtr->get()[ts_subtr]->GetData().shape(1),
            itsAvgResultSubtr->get()[ts_subtr]->GetData().shape(0));
        casacore::Cube<casacore::Complex> cube_residual(
            cube_shape, itsAvgResultSubtr->get()[ts_subtr]->GetData().data(),
            casacore::SHARE);
        base::cursor<std::complex<float>> cr_residual =
            base::casa_cursor(cube_residual);

        // Construct a cursor to iterate over a slice of the mixing matrix
        // at the resolution of the residual. The "to" and "from" direction
        // are fixed. Since the full mixing matrix is 5-D, the slice is
        // therefore 3-D. Each individual value in the slice quantifies the
        // influence of the source to subtract on the target direction for
        // a particular correlation, channel, and baseline.
        //
        // The target direction is the direction with the highest index by
        // convention, i.e. index itsNDir - 1. The directions to subtract
        // have the lowest indices by convention, i.e. indices
        // [0, nDrSubtr).
        const std::array<size_t, 5> shape = itsFactorsSubtr[ts_subtr].shape();
        size_t stride_mix_subtr_slice[3] = {
            static_cast<size_t>(itsNDir * itsNDir),
            static_cast<size_t>(itsNDir * itsNDir * nCr),
            static_cast<size_t>(itsNDir * itsNDir * nCr * nChSubtr)};
        assert(shape[4] == itsNDir && shape[3] == itsNDir && shape[2] == nCr &&
               shape[1] == nChSubtr);

        const casacore::Array<casacore::DComplex> demix_factors_subtr(
            IPosition(5, shape[4], shape[3], shape[2], shape[1], shape[0]),
            itsFactorsSubtr[ts_subtr].data(), casacore::SHARE);

        IPosition offset(5, itsNDir - 1, dr, 0, 0, 0);
        base::const_cursor<std::complex<double>> cr_mix_subtr(
            &(demix_factors_subtr(offset)), 3, stride_mix_subtr_slice);

        // Subtract the source.
        float variance_before, variance_after;
        subtract(nBl, nChSubtr, cr_baseline, cr_residual, cr_model_subtr,
                 cr_mix_subtr, variance_before, variance_after);
        storage.variance_before += variance_before;
        storage.variance_after += variance_after;
      }
    }

    // Copy solutions to global solution array.
    std::copy(storage.unknowns.begin(), storage.unknowns.end(),
              &(itsUnknowns[(itsTimeIndex + ts) * nDr * nSt * 8]));
  });

  // Store last known solutions.
  if (itsPropagateSolutions && nTime > 0) {
    std::copy(&(itsUnknowns[(itsTimeIndex + nTime - 1) * nDr * nSt * 8]),
              &(itsUnknowns[(itsTimeIndex + nTime) * nDr * nSt * 8]),
              itsPrevSolution.begin());
  }

  // Update convergence count.
  for (size_t i = 0; i < nThread; ++i) {
    itsNConverged += threadStorage[i].count_converged;
  }

  // Calculate data before/after variance ratio
  float variance_before = 0.0f;
  float variance_after = 0.0f;
  for (std::vector<ThreadPrivateStorage>::iterator it = threadStorage.begin(),
                                                   end = threadStorage.end();
       it != end; ++it) {
    variance_before += it->variance_before;
    variance_after += it->variance_after;
  }
  // Update self
  itsVarianceRatio += (variance_before / (variance_after + 1e-6f));
}

void Demixer::dumpSolutions() {
  // Construct solution grid.
  const std::vector<double>& freq = getInfo().chanFreqs();
  const std::vector<double>& freqWidth = getInfo().chanWidths();
  parmdb::Axis::ShPtr freqAxis(
      new parmdb::RegularAxis(freq[0] - freqWidth[0] * 0.5, freqWidth[0], 1));
  parmdb::Axis::ShPtr timeAxis(new parmdb::RegularAxis(
      getInfo().startTime() - getInfo().timeInterval() * 0.5,
      itsTimeIntervalAvg, itsNTimeDemix));
  parmdb::Grid solGrid(freqAxis, timeAxis);

  // Create and initialize ParmDB.
  parmdb::ParmDB parmDB(parmdb::ParmDBMeta("casa", itsInstrumentName), true);
  parmdb::ParmSet parmSet;
  parmdb::ParmCache parmCache(parmSet, solGrid.getBoundingBox());

  // Store the (freq, time) resolution of the solutions.
  std::vector<double> resolution(2);
  resolution[0] = freqWidth[0];
  resolution[1] = itsTimeIntervalAvg;
  parmDB.setDefaultSteps(resolution);

  // Map station indices in the solution array to the corresponding antenna
  // names. This is required because solutions are only produced for
  // stations that participate in one or more baselines. Due to the baseline
  // selection or missing baselines, solutions may be available for less
  // than the total number of station available in the observation.
  const DPInfo& info = itsFilter->getInfo();
  const std::vector<int>& antennaUsed = info.antennaUsed();
  const std::vector<std::string>& antennaNames = info.antennaNames();

  std::vector<parmdb::Parm> parms;
  for (size_t dr = 0; dr < itsNModel; ++dr) {
    for (size_t st = 0; st < itsNStation; ++st) {
      std::string name(antennaNames[antennaUsed[st]]);
      std::string suffix(name + ":" + itsAllSources[dr]);

      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:0:0:Real:" + suffix)));
      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:0:0:Imag:" + suffix)));

      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:0:1:Real:" + suffix)));
      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:0:1:Imag:" + suffix)));

      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:1:0:Real:" + suffix)));
      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:1:0:Imag:" + suffix)));

      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:1:1:Real:" + suffix)));
      parms.push_back(parmdb::Parm(
          parmCache,
          parmSet.addParm(parmDB, "DirectionalGain:1:1:Imag:" + suffix)));
    }
  }

  // Cache parameter values.
  parmCache.cacheValues();

  // Assign solution grid to parameters.
  for (size_t i = 0; i < parms.size(); ++i) {
    parms[i].setSolveGrid(solGrid);
  }

  // Write solutions.
  for (size_t ts = 0; ts < itsNTimeDemix; ++ts) {
    double* unknowns = &(itsUnknowns[ts * itsNModel * itsNStation * 8]);
    for (size_t i = 0; i < parms.size(); ++i) {
      parms[i].setCoeff(parmdb::Grid::Location(0, ts), unknowns + i, 1);
    }
  }

  // Flush solutions to disk.
  parmCache.flush();
}

void Demixer::addToMS(const std::string& msName) {
  Step::addToMS(msName);
  casacore::Table histtab(msName + "/HISTORY", casacore::Table::Update);
  casacore::ScalarColumn<casacore::String> message(histtab, "MESSAGE");
  casacore::ScalarColumn<casacore::String> application(histtab, "APPLICATION");
  casacore::ArrayColumn<casacore::String> parms(histtab, "APP_PARAMS");

  unsigned int n_row = histtab.nrow();
  histtab.addRow();
  application.put(n_row, "DP3");
  message.put(n_row, "Noise ratio before/after demixing");
  casacore::Vector<casacore::String> appvec(1);
  appvec[0] = std::to_string(itsVarianceRatio);
  parms.put(n_row, appvec);
}

namespace {
std::string toString(double value) {
  std::ostringstream os;
  os << std::setprecision(16) << value;
  return os.str();
}
}  // end unnamed namespace

}  // namespace steps
}  // namespace dp3
