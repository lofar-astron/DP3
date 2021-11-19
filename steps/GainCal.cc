// GainCal.cc: DPPP step class to do a gain calibration
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "GainCal.h"
#include "ApplyCal.h"
#include "ColumnReader.h"

#include <Version.h>

#include "../base/Simulate.h"
#include "../base/PhaseFitter.h"
#include "../base/CalType.h"
#include "../base/CursorUtilCasa.h"
#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/DPLogger.h"
#include "../base/Exceptions.h"
#include "../base/SourceDBUtil.h"

#include "../parmdb/ParmDB.h"
#include "../parmdb/ParmValue.h"
#include "../parmdb/SourceDB.h"

#include "../common/ParameterSet.h"
#include "../common/StringTools.h"

#include <aocommon/threadpool.h>

#include <fstream>
#include <ctime>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/OS/File.h>

#include <boost/make_unique.hpp>

#include <vector>
#include <algorithm>

#include <limits>
#include <iostream>
#include <iomanip>

using casacore::Cube;
using casacore::IPosition;
using casacore::Matrix;
using casacore::Table;

using dp3::base::CalType;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::base::GainCalAlgorithm;

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace dp3 {
namespace steps {

GainCal::GainCal(InputStep& input, const common::ParameterSet& parset,
                 const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsUseModelColumn(parset.getBool(prefix + "usemodelcolumn", false)),
      itsModelColumnName(),
      itsParmDBName(parset.getString(prefix + "parmdb", "")),
      itsUseH5Parm(itsParmDBName.find(".h5") != string::npos),
      itsDebugLevel(parset.getInt(prefix + "debuglevel", 0)),
      itsDetectStalling(parset.getBool(prefix + "detectstalling", true)),
      itsApplySolution(parset.getBool(prefix + "applysolution", false)),
      itsUVWFlagStep(&input, parset, prefix),
      itsParallelFor(1),
      itsThreadPool(1),
      itsFirstSubStep(),
      itsResultStep(std::make_shared<ResultStep>()),
      itsBaselineSelection(parset, prefix),
      itsMaxIter(parset.getInt(prefix + "maxiter", 50)),
      itsTolerance(parset.getDouble(prefix + "tolerance", 1.e-5)),
      itsPropagateSolutions(
          parset.getBool(prefix + "propagatesolutions", true)),
      itsSolInt(parset.getInt(prefix + "solint", 1)),
      itsNFreqCells(0),
      itsConverged(0),
      itsNonconverged(0),
      itsFailed(0),
      itsStalled(0),
      itsStepInParmUpdate(0),
      itsChunkStartTime(0),
      itsStepInSolInt(0),
      itsAllSolutions() {
  std::stringstream ss;
  ss << parset;
  itsParsetString = ss.str();

  if (itsParmDBName == "") {
    itsParmDBName = parset.getString("msin") + "/instrument";
  }

  if (!itsUseH5Parm) {
    itsTimeSlotsPerParmUpdate =
        parset.getInt(prefix + "timeslotsperparmupdate", 500);
  } else {
    itsTimeSlotsPerParmUpdate = 0;
  }

  itsDataResultStep = std::make_shared<ResultStep>();
  itsUVWFlagStep.setNextStep(itsDataResultStep);

  if (!itsUseModelColumn) {
    auto predict_step = boost::make_unique<Predict>(input, parset, prefix);
    predict_step->SetThreadData(itsThreadPool, nullptr);
    predict_step->setNextStep(itsResultStep);
    itsFirstSubStep = std::move(predict_step);
  } else {
    // Remain compatible with the old situation, where the input step read the
    // model data and the input step had the model column name.
    const std::string column_key = prefix + "modelcolumn";
    if (parset.isDefined("msin.modelcolumn")) {
      if (parset.isDefined(column_key)) {
        throw std::runtime_error(
            "The input contains both the deprecated msin.modelcolumn setting "
            "and the " +
            column_key +
            " setting. Please remove the deprecated setting from the input.");
      }
      DPLOG_WARN_STR(
          "Warning: The input contains the deprecated msin.modelcolumn "
          "setting. Please use " +
          column_key + " instead.");
      itsModelColumnName = parset.getString("msin.modelcolumn");
    } else {
      itsModelColumnName = parset.getString(column_key, "MODEL_DATA");
    }
    itsApplyBeamToModelColumn =
        parset.getBool(prefix + "applybeamtomodelcolumn", false);

    auto column_reader_step = boost::make_unique<ColumnReader>(
        input, parset, prefix, itsModelColumnName);
    if (itsApplyBeamToModelColumn) {
      auto apply_beam_step =
          std::make_shared<ApplyBeam>(&input, parset, prefix, true);
      column_reader_step->setNextStep(apply_beam_step);
      apply_beam_step->setNextStep(itsResultStep);
    } else {
      column_reader_step->setNextStep(itsResultStep);
    }
    itsFirstSubStep = std::move(column_reader_step);
  }

  itsNIter.resize(4, 0);

  if (itsApplySolution) {
    itsBuf.resize(itsSolInt);
  } else {
    itsBuf.resize(1);
  }

  string modestr = parset.getString(prefix + "caltype");
  itsMode = base::StringToCalType(modestr);
  unsigned int defaultNChan = 0;
  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    defaultNChan = 1;
  } else if (itsMode == CalType::kTecScreen)
    throw std::runtime_error("Can't solve with mode TECSCREEN");
  itsNChan = parset.getInt(prefix + "nchan", defaultNChan);
}

GainCal::~GainCal() {}

void GainCal::setAntennaUsed() {
  Matrix<bool> selbl(itsBaselineSelection.apply(info()));
  unsigned int nBl = info().getAnt1().size();
  itsAntennaUsed.resize(info().antennaNames().size());
  itsAntennaUsed = false;
  for (unsigned int bl = 0; bl < nBl; ++bl) {
    if (selbl(info().getAnt1()[bl], info().getAnt2()[bl])) {
      itsAntennaUsed[info().getAnt1()[bl]] = true;
      itsAntennaUsed[info().getAnt2()[bl]] = true;
    }
  }
}

void GainCal::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  info().setNeedVisData();

  // By giving a thread pool to the predicter, the threads are
  // sustained.
  itsThreadPool.SetNThreads(info().nThreads());
  itsParallelFor.SetNThreads(info().nThreads());
  itsUVWFlagStep.updateInfo(infoIn);

  itsFirstSubStep->setInfo(infoIn);

  if (itsApplySolution) {
    info().setWriteData();
    info().setWriteFlags();
  }

  if (itsSolInt == 0) {
    itsSolInt = info().ntime();
  }
  if (itsTimeSlotsPerParmUpdate == 0) {
    itsTimeSlotsPerParmUpdate = info().ntime();
  }

  if (itsNChan == 0) {
    itsNChan = info().nchan();
  }
  if (itsNChan > info().nchan()) {
    itsNChan = info().nchan();
  }
  itsNFreqCells = info().nchan() / itsNChan;
  if (itsNChan * itsNFreqCells <
      info().nchan()) {  // If last freq cell is smaller
    itsNFreqCells++;
  }

  itsSols.reserve(itsTimeSlotsPerParmUpdate);

  itsSelectedBL = itsBaselineSelection.applyVec(info());
  setAntennaUsed();

  // Compute average frequency for every freqcell
  itsFreqData.resize(itsNFreqCells);
  for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
    double meanfreq = 0;
    unsigned int chmin = itsNChan * freqCell;
    unsigned int chmax = std::min(info().nchan(), chmin + itsNChan);

    meanfreq = std::accumulate(info().chanFreqs().data() + chmin,
                               info().chanFreqs().data() + chmax, 0.0);

    itsFreqData[freqCell] = meanfreq / (chmax - chmin);
  }

  // Initialize phase fitters, set their frequency data
  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    itsTECSols.reserve(itsTimeSlotsPerParmUpdate);

    itsPhaseFitters.reserve(
        itsNFreqCells);  // TODO: could be numthreads instead

    unsigned int nSt = info().antennaUsed().size();
    for (unsigned int st = 0; st < nSt; ++st) {
      itsPhaseFitters.push_back(boost::make_unique<PhaseFitter>());
      itsPhaseFitters[st]->Initialize(itsFreqData);
    }
  }

  iS.reserve(itsNFreqCells);
  unsigned int chMax = itsNChan;
  for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
    if ((freqCell + 1) * itsNChan >
        info().nchan()) {  // Last cell can be smaller
      chMax -= ((freqCell + 1) * itsNChan) % info().nchan();
    }

    GainCalAlgorithm::Mode smode;
    switch (itsMode) {
      case CalType::kDiagonal:
        smode = GainCalAlgorithm::DEFAULT;
        break;
      case CalType::kFullJones:
        smode = GainCalAlgorithm::FULLJONES;
        break;
      case CalType::kScalarPhase:
      case CalType::kDiagonalPhase:
      case CalType::kTec:
      case CalType::kTecAndPhase:
        smode = GainCalAlgorithm::PHASEONLY;
        break;
      case CalType::kDiagonalAmplitude:
      case CalType::kScalarAmplitude:
        smode = GainCalAlgorithm::AMPLITUDEONLY;
        break;
      default:
        throw Exception("Unhandled mode");
    }

    iS.emplace_back(GainCalAlgorithm(
        itsSolInt, chMax, smode, scalarMode(itsMode), itsTolerance,
        info().antennaUsed().size(), itsDetectStalling, itsDebugLevel));
  }

  itsFlagCounter.init(getInfo());

  itsChunkStartTime = info().startTime();

  if (itsDebugLevel > 0) {
    if (getInfo().nThreads() != 1)
      throw std::runtime_error("nthreads should be 1 in debug mode");
    assert(itsTimeSlotsPerParmUpdate >= info().ntime());
    itsAllSolutions.resize(IPosition(
        6, iS[0].numCorrelations(), info().antennaUsed().size(),
        (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) ? 2 : 1,
        itsNFreqCells, itsMaxIter, info().ntime()));
  }
}

void GainCal::show(std::ostream& os) const {
  os << "GainCal " << itsName << '\n';
  if (itsUseH5Parm) {
    os << "  H5Parm:              " << itsParmDBName;
  } else {
    os << "  parmdb:              " << itsParmDBName;
    if (Table::isReadable(itsParmDBName)) {
      os << " (existing)";
    } else {
      os << " (will be created)";
    }
  }
  os << '\n';
  os << "  solint:              " << itsSolInt << '\n';
  os << "  nchan:               " << itsNChan << '\n';
  os << "  max iter:            " << itsMaxIter << '\n';
  os << "  tolerance:           " << itsTolerance << '\n';
  os << "  caltype:             " << ToString(itsMode) << '\n';
  os << "  apply solution:      " << std::boolalpha << itsApplySolution << '\n';
  os << "  propagate solutions: " << std::boolalpha << itsPropagateSolutions
     << '\n';
  if (!itsUseH5Parm) {
    os << "  timeslotsperparmupdate: " << itsTimeSlotsPerParmUpdate << '\n';
  }
  os << "  detect stalling:     " << std::boolalpha << itsDetectStalling
     << '\n';
  os << "  use model column:    " << std::boolalpha << itsUseModelColumn
     << '\n';
  os << "  model column name:   " << itsModelColumnName << '\n';
  itsBaselineSelection.show(os);
  for (Step* step = itsFirstSubStep.get(); step != nullptr;
       step = step->getNextStep().get()) {
    step->show(os);
  }
  itsUVWFlagStep.show(os);
}

void GainCal::showTimings(std::ostream& os, double duration) const {
  double totaltime = itsTimer.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " GainCal " << itsName << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerPredict.getElapsed(), totaltime);
  os << " of it spent in predict" << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerFill.getElapsed(), totaltime);
  os << " of it spent in reordering visibility data" << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerSolve.getElapsed(), totaltime);
  os << " of it spent in estimating gains and computing residuals" << '\n';

  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    os << "          ";
    FlagCounter::showPerc1(os, itsTimerPhaseFit.getElapsed(), totaltime);
    os << " of it spent in fitting phases" << '\n';
  }

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerWrite.getElapsed(), totaltime);
  os << " of it spent in writing gain solutions to disk" << '\n';

  os << "        ";
  os << "Converged: " << itsConverged << ", stalled: " << itsStalled
     << ", non converged: " << itsNonconverged << ", failed: " << itsFailed
     << '\n';
  os << "        ";
  os << "Iters converged: "
     << (itsConverged == 0 ? 0 : itsNIter[0] / itsConverged);
  os << ", stalled: " << (itsStalled == 0 ? 0 : itsNIter[1] / itsStalled);
  os << ", non converged: "
     << (itsNonconverged == 0 ? 0 : itsNIter[2] / itsNonconverged);
  os << ", failed: " << (itsFailed == 0 ? 0 : itsNIter[3] / itsFailed) << '\n';
}

bool GainCal::process(const DPBuffer& bufin) {
  itsTimer.start();

  unsigned int bufIndex = 0;

  if (itsApplySolution) {
    // Need to keep a copy of all solint buffers in this step
    bufIndex = itsStepInSolInt;
    itsBuf[bufIndex].copy(bufin);
  } else {
    // We'll read the necessary info from the buffer and pass it on
    itsBuf[bufIndex].referenceFilled(bufin);
  }
  itsInput.fetchUVW(bufin, itsBuf[bufIndex], itsTimer);
  itsInput.fetchWeights(bufin, itsBuf[bufIndex], itsTimer);
  itsInput.fetchFullResFlags(bufin, itsBuf[bufIndex], itsTimer);

  // UVW flagging happens on a copy of the buffer, so these flags are not
  // written
  itsUVWFlagStep.process(itsBuf[bufIndex]);

  casacore::Complex* data = itsBuf[bufIndex].getData().data();
  float* weight = itsBuf[bufIndex].getWeights().data();
  const bool* flag = itsBuf[bufIndex].getFlags().data();

  // Simulate.
  //
  // Model visibilities for each direction of interest will be computed
  // and stored.

  itsTimerPredict.start();
  itsFirstSubStep->process(itsBuf[bufIndex]);
  itsTimerPredict.stop();

  itsTimerFill.start();

  if (itsStepInSolInt == 0) {
    // Start new solution interval

    for (unsigned int freqCell = 0; freqCell < itsNFreqCells; freqCell++) {
      iS[freqCell].clearStationFlagged();
      iS[freqCell].resetVis();
    }
  }

  // Store data in the GainCalAlgorithm object
  fillMatrices(itsResultStep->get().getData().data(), data, weight, flag);
  itsTimerFill.stop();

  ++itsStepInSolInt;
  if (itsStepInSolInt == itsSolInt) {
    // Solve past solution interval
    calibrate();
    itsStepInParmUpdate++;

    if (itsApplySolution) {
      Cube<casacore::Complex> invsol = invertSol(itsSols.back());
      for (unsigned int stepInSolInt = 0; stepInSolInt < itsSolInt;
           stepInSolInt++) {
        applySolution(itsBuf[stepInSolInt], invsol);
        getNextStep()->process(itsBuf[stepInSolInt]);
      }
    }

    itsStepInSolInt = 0;
  }

  itsTimer.stop();

  if (!itsUseH5Parm && (itsStepInParmUpdate == itsTimeSlotsPerParmUpdate)) {
    writeSolutionsParmDB(itsChunkStartTime);
    itsChunkStartTime +=
        itsSolInt * itsTimeSlotsPerParmUpdate * info().timeInterval();
    itsSols.clear();
    itsTECSols.clear();
    itsStepInParmUpdate = 0;
  }

  if (!itsApplySolution) {
    getNextStep()->process(itsBuf[bufIndex]);
  }
  return false;
}

Cube<casacore::Complex> GainCal::invertSol(const Cube<casacore::Complex>& sol) {
  Cube<casacore::Complex> invsol = sol.copy();
  unsigned int nCr = invsol.shape()[0];

  // Invert copy of solutions
  unsigned int nSt = invsol.shape()[1];
  for (unsigned int st = 0; st < nSt; ++st) {
    for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
      if (nCr == 4) {
        ApplyCal::invert(&invsol(0, st, freqCell));
      } else {
        for (unsigned int cr = 0; cr < nCr; ++cr) {
          invsol(cr, st, freqCell) = 1.0f / invsol(cr, st, freqCell);
        }
      }
    }
  }

  return invsol;
}

void GainCal::applySolution(DPBuffer& buf,
                            const Cube<casacore::Complex>& invsol) {
  unsigned int nbl = buf.getData().shape()[2];
  casacore::Complex* data = buf.getData().data();
  float* weight = buf.getWeights().data();  // Not initialized yet
  bool* flag = buf.getFlags().data();
  unsigned int nchan = buf.getData().shape()[1];

  unsigned int nCr = invsol.shape()[0];

  for (size_t bl = 0; bl < nbl; ++bl) {
    for (size_t chan = 0; chan < nchan; chan++) {
      unsigned int antA = info().antennaMap()[info().getAnt1()[bl]];
      unsigned int antB = info().antennaMap()[info().getAnt2()[bl]];
      unsigned int freqCell = chan / itsNChan;
      if (nCr > 2) {
        ApplyCal::applyFull(
            &invsol(0, antA, freqCell), &invsol(0, antB, freqCell),
            &data[bl * 4 * nchan + chan * 4],
            &weight[bl * 4 * nchan +
                    chan * 4],  // Not passing weights, any pointer should do
            &flag[bl * 4 * nchan + chan * 4], bl, chan, false,
            itsFlagCounter);  // Update weights is disabled here
      } else if (scalarMode(itsMode)) {
        ApplyCal::applyScalar(
            &invsol(0, antA, freqCell), &invsol(0, antB, freqCell),
            &data[bl * 4 * nchan + chan * 4],
            &weight[bl * 4 * nchan +
                    chan * 4],  // Not passing weights, any pointer should do
            &flag[bl * 4 * nchan + chan * 4], bl, chan, false,
            itsFlagCounter);  // Update weights is disabled here
      } else {
        ApplyCal::applyDiag(
            &invsol(0, antA, freqCell), &invsol(0, antB, freqCell),
            &data[bl * 4 * nchan + chan * 4],
            &weight[bl * 4 * nchan +
                    chan * 4],  // Not passing weights, any pointer should do
            &flag[bl * 4 * nchan + chan * 4], bl, chan, false,
            itsFlagCounter);  // Update weights is disabled here
      }
    }
  }
}

// Fills itsVis and itsMVis as matrices with all 00 polarizations in the
// top left, all 11 polarizations in the bottom right, etc.
// For TEC fitting, it also sets weights for the frequency cells
void GainCal::fillMatrices(casacore::Complex* model, casacore::Complex* data,
                           float* weight, const casacore::Bool* flag) {
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = info().ncorr();
  assert(nCr == 4 || nCr == 2 || nCr == 1);

  for (unsigned int ch = 0; ch < nCh; ++ch) {
    for (unsigned int bl = 0; bl < nBl; ++bl) {
      if (itsSelectedBL[bl]) {
        int ant1 = info().antennaMap()[info().getAnt1()[bl]];
        int ant2 = info().antennaMap()[info().getAnt2()[bl]];
        assert(ant1 >= 0 && ant2 >= 0);
        if (ant1 == ant2 || iS[ch / itsNChan].getStationFlagged()[ant1] ||
            iS[ch / itsNChan].getStationFlagged()[ant2] ||
            flag[bl * nCr * nCh + ch * nCr]) {  // Only check flag of cr==0
          continue;
        }

        if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
          iS[ch / itsNChan].incrementWeight(weight[bl * nCr * nCh + ch * nCr]);
        }

        for (unsigned int cr = 0; cr < nCr; ++cr) {
          // The nCrDiv is there such that for nCr==2 the visibilities end up at
          // (0,0) for cr==0, (1,1) for cr==1
          unsigned int nCrDiv = (nCr == 4 ? 2 : 1);
          iS[ch / itsNChan].getVis()(IPosition(6, ant1, cr / nCrDiv,
                                               itsStepInSolInt, ch % itsNChan,
                                               cr % 2, ant2)) =
              casacore::DComplex(data[bl * nCr * nCh + ch * nCr + cr]) *
              casacore::DComplex(sqrt(weight[bl * nCr * nCh + ch * nCr + cr]));
          iS[ch / itsNChan].getMVis()(IPosition(6, ant1, cr / nCrDiv,
                                                itsStepInSolInt, ch % itsNChan,
                                                cr % 2, ant2)) =
              casacore::DComplex(model[bl * nCr * nCh + ch * nCr + cr]) *
              casacore::DComplex(sqrt(weight[bl * nCr * nCh + ch * nCr + cr]));

          // conjugate transpose
          iS[ch / itsNChan].getVis()(IPosition(6, ant2, cr % 2, itsStepInSolInt,
                                               ch % itsNChan, cr / nCrDiv,
                                               ant1)) =
              casacore::DComplex(conj(data[bl * nCr * nCh + ch * nCr + cr])) *
              casacore::DComplex(sqrt(weight[bl * nCr * nCh + ch * nCr + cr]));
          iS[ch / itsNChan].getMVis()(IPosition(6, ant2, cr % 2,
                                                itsStepInSolInt, ch % itsNChan,
                                                cr / nCrDiv, ant1)) =
              casacore::DComplex(conj(model[bl * nCr * nCh + ch * nCr + cr])) *
              casacore::DComplex(sqrt(weight[bl * nCr * nCh + ch * nCr + cr]));
        }
      }
    }
  }
}

bool GainCal::scalarMode(CalType caltype) {
  return (caltype == CalType::kTecAndPhase || caltype == CalType::kTec ||
          caltype == CalType::kScalarPhase ||
          caltype == CalType::kScalarAmplitude);
}

bool GainCal::diagonalMode(CalType caltype) {
  return (caltype == CalType::kDiagonal || caltype == CalType::kDiagonalPhase ||
          caltype == CalType::kDiagonalAmplitude);
}

void GainCal::calibrate() {
  itsTimerSolve.start();

  for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
    if (itsPropagateSolutions) {
      iS[freqCell].init(false);
    } else {
      iS[freqCell].init(true);
    }
  }

  unsigned int iter = 0;

  casacore::Matrix<double> tecsol(itsMode == CalType::kTecAndPhase ? 2 : 1,
                                  info().antennaUsed().size(), 0);

  std::vector<GainCalAlgorithm::Status> converged(
      itsNFreqCells, GainCalAlgorithm::NOTCONVERGED);

  for (; iter < itsMaxIter; ++iter) {
    bool allConverged = true;
    itsParallelFor.Run(
        0, itsNFreqCells, [&](size_t freqCell, size_t /*thread*/) {
          // Do another step when stalled and not all converged
          if (converged[freqCell] != GainCalAlgorithm::CONVERGED) {
            converged[freqCell] = iS[freqCell].doStep(iter);
            // Only continue if there are steps worth continuing
            // (so not converged, failed or stalled)
            if (converged[freqCell] == GainCalAlgorithm::NOTCONVERGED) {
              allConverged = false;
            }
          }
        });

    if (itsDebugLevel > 0) {
      for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
        Matrix<casacore::DComplex> fullSolution =
            iS[freqCell].getSolution(false);
        std::copy(fullSolution.begin(), fullSolution.end(),
                  &(itsAllSolutions(IPosition(6, 0, 0, 0, freqCell, iter,
                                              itsStepInParmUpdate))));
      }
    }

    if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
      itsTimerSolve.stop();
      itsTimerPhaseFit.start();
      casacore::Matrix<casacore::DComplex> sols_f(itsNFreqCells,
                                                  info().antennaUsed().size());

      unsigned int nSt = info().antennaUsed().size();

      // TODO: set phase reference to something smarter than station 0
      for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
        casacore::Matrix<casacore::DComplex> sol =
            iS[freqCell].getSolution(false);
        if (iS[freqCell].getStationFlagged()[0]) {
          // If reference station flagged, flag whole channel
          for (unsigned int st = 0; st < info().antennaUsed().size(); ++st) {
            iS[freqCell].getStationFlagged()[st] = true;
          }
        } else {
          for (unsigned int st = 0; st < info().antennaUsed().size(); ++st) {
            sols_f(freqCell, st) = sol(st, 0) / sol(0, 0);
            assert(casacore::isFinite(sols_f(freqCell, st)));
          }
        }
      }

      // Fit the data for each station
      itsParallelFor.Run(0, nSt, [&](size_t st, size_t /*thread*/) {
        unsigned int numpoints = 0;
        double* phases = itsPhaseFitters[st]->PhaseData();
        double* weights = itsPhaseFitters[st]->WeightData();
        for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
          if (iS[freqCell].getStationFlagged()[st % nSt] ||
              converged[freqCell] == GainCalAlgorithm::FAILED) {
            phases[freqCell] = 0;
            weights[freqCell] = 0;
          } else {
            phases[freqCell] = arg(sols_f(freqCell, st));
            if (!std::isfinite(phases[freqCell])) {
              std::cout << "Yuk, phases[freqCell]=" << phases[freqCell]
                        << ", sols_f(freqCell, st)=" << sols_f(freqCell, st)
                        << '\n';
              assert(std::isfinite(phases[freqCell]));
            }
            assert(iS[freqCell].getWeight() > 0);
            weights[freqCell] = iS[freqCell].getWeight();
            numpoints++;
          }
        }

        for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
          assert(std::isfinite(phases[freqCell]));
        }

        if (numpoints > 1) {  // TODO: limit should be higher
          // cout<<"tecsol(0,"<<st<<")="<<tecsol(0,st)<<",
          // tecsol(1,"<<st<<")="<<tecsol(1,st)<<'\n';
          if (itsMode == CalType::kTecAndPhase) {
            itsPhaseFitters[st]->FitDataToTEC2Model(tecsol(0, st),
                                                    tecsol(1, st));
          } else {  // itsMode==kTec
            itsPhaseFitters[st]->FitDataToTEC1Model(tecsol(0, st));
          }
          // Update solution in GainCalAlgorithm object
          for (unsigned int freqCell = 0; freqCell < itsNFreqCells;
               ++freqCell) {
            assert(std::isfinite(phases[freqCell]));
            iS[freqCell].getSolution(false)(st, 0) =
                std::polar(1., phases[freqCell]);
          }
        } else {
          tecsol(0, st) = 0;  // std::numeric_limits<double>::quiet_NaN();
          if (itsMode == CalType::kTecAndPhase) {
            tecsol(1, st) = 0;  // std::numeric_limits<double>::quiet_NaN();
          }
        }

        if (itsDebugLevel > 0) {
          for (unsigned int freqCell = 0; freqCell < itsNFreqCells;
               ++freqCell) {
            Matrix<casacore::DComplex> fullSolution =
                iS[freqCell].getSolution(false);
            std::copy(fullSolution.begin(), fullSolution.end(),
                      &(itsAllSolutions(IPosition(6, 0, 0, 1, freqCell, iter,
                                                  itsStepInParmUpdate))));
          }
        }
      });
      itsTimerPhaseFit.stop();
      itsTimerSolve.start();
    }

    if (allConverged) {
      break;
    }

  }  // End niter

  for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
    switch (converged[freqCell]) {
      case GainCalAlgorithm::CONVERGED: {
        itsConverged++;
        itsNIter[0] += iter;
        break;
      }
      case GainCalAlgorithm::STALLED: {
        itsStalled++;
        itsNIter[1] += iter;
        break;
      }
      case GainCalAlgorithm::NOTCONVERGED: {
        itsNonconverged++;
        itsNIter[2] += iter;
        break;
      }
      case GainCalAlgorithm::FAILED: {
        itsFailed++;
        itsNIter[3] += iter;
        break;
      }
      default:
        throw Exception("Unknown converged status");
    }
  }

  // Calibrate terminated (either by maxiter or by converging)

  Cube<casacore::Complex> sol(iS[0].numCorrelations(),
                              info().antennaUsed().size(), itsNFreqCells);

  unsigned int transpose[2][4] = {{0, 1, 0, 0}, {0, 2, 1, 3}};

  unsigned int nSt = info().antennaUsed().size();

  for (unsigned int freqCell = 0; freqCell < itsNFreqCells; ++freqCell) {
    casacore::Matrix<casacore::DComplex> tmpsol =
        iS[freqCell].getSolution(true);

    for (unsigned int st = 0; st < nSt; st++) {
      for (unsigned int cr = 0; cr < iS[0].nCr(); ++cr) {
        unsigned int crt = transpose[iS[0].numCorrelations() / 4]
                                    [cr];  // Conjugate transpose ! (only for
                                           // numCorrelations = 4)
        sol(crt, st, freqCell) = conj(tmpsol(st, cr));  // Conjugate transpose
        if (itsMode == CalType::kDiagonal ||
            itsMode == CalType::kDiagonalPhase ||
            itsMode == CalType::kDiagonalAmplitude) {
          sol(crt + 1, st, freqCell) =
              conj(tmpsol(st + nSt, cr));  // Conjugate transpose
        }
      }
    }
  }
  itsSols.push_back(sol);
  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    itsTECSols.push_back(tecsol);
  }

  itsTimerSolve.stop();
}  // End calibrate()

void GainCal::initParmDB() {
  itsParmDB = std::shared_ptr<parmdb::ParmDB>(
      new parmdb::ParmDB(parmdb::ParmDBMeta("casa", itsParmDBName), false));
  itsParmDB->lock();
  // Store the (freq, time) resolution of the solutions.

  double freqWidth = getInfo().chanWidths()[0];
  if (getInfo().chanFreqs().size() >
      1) {  // Handle data with evenly spaced gaps between channels
    freqWidth = info().chanFreqs()[1] - info().chanFreqs()[0];
  }

  std::vector<double> resolution(2);
  resolution[0] = freqWidth * itsNChan;
  resolution[1] = info().timeInterval() * itsSolInt;
  itsParmDB->setDefaultSteps(resolution);

  string parmname = parmName() + "*";

  if (!itsParmDB->getNames(parmname).empty()) {
    DPLOG_WARN_STR("Solutions for " << parmname << " already in "
                                    << itsParmDBName << ", these are removed");
    // Specify entire domain of this MS; only to specify that the existing
    // values should be deleted for this domain
    parmdb::Axis::ShPtr tdomAxis(new parmdb::RegularAxis(
        info().startTime(), info().ntime() * info().timeInterval(), 1));
    parmdb::Axis::ShPtr fdomAxis(
        new parmdb::RegularAxis(info().chanFreqs()[0] - freqWidth * 0.5,
                                freqWidth * getInfo().chanFreqs().size(), 1));

    itsParmDB->deleteValues(
        parmname, parmdb::Box(fdomAxis->start(), tdomAxis->start(),
                              fdomAxis->end(), tdomAxis->end(), true));
  }

  // Write out default values, if they don't exist yet
  parmdb::ParmMap parmset;

  // Write out default amplitudes
  if (itsMode == CalType::kDiagonalPhase || itsMode == CalType::kScalarPhase) {
    itsParmDB->getDefValues(parmset, "Gain:0:0:Ampl");
    if (parmset.empty()) {
      parmdb::ParmValueSet pvset(parmdb::ParmValue(1.0));
      itsParmDB->putDefValue("Gain:0:0:Ampl", pvset);
      itsParmDB->putDefValue("Gain:1:1:Ampl", pvset);
    }
  }

  // Write out default phases
  if (itsMode == CalType::kDiagonalAmplitude ||
      itsMode == CalType::kScalarAmplitude) {
    itsParmDB->getDefValues(parmset, "Gain:0:0:Phase");
    if (parmset.empty()) {
      parmdb::ParmValueSet pvset(parmdb::ParmValue(0.0));
      itsParmDB->putDefValue("Gain:0:0:Phase", pvset);
      itsParmDB->putDefValue("Gain:1:1:Phase", pvset);
    }
  }

  // Write out default gains
  if (itsMode == CalType::kDiagonal || itsMode == CalType::kFullJones) {
    itsParmDB->getDefValues(parmset, "Gain:0:0:Real");
    if (parmset.empty()) {
      parmdb::ParmValueSet pvset(parmdb::ParmValue(1.0));
      itsParmDB->putDefValue("Gain:0:0:Real", pvset);
      itsParmDB->putDefValue("Gain:1:1:Real", pvset);
    }
  }
}

std::string GainCal::parmName() {
  std::string name;
  if (itsMode == CalType::kScalarPhase) {
    name = "CommonScalarPhase:";
  } else if (itsMode == CalType::kScalarAmplitude) {
    name = "CommonScalarAmplitude:";
  } else if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    name = "TEC:";
  } else {
    name = "Gain:";
  }
  return name;
}

void GainCal::writeSolutionsH5Parm(double) {
  itsTimer.start();
  itsTimerWrite.start();

  H5Parm h5parm(itsParmDBName, true);

  // Fill antenna info in H5Parm, need to convert from casa types to std types
  std::vector<std::string> allAntennaNames(info().antennaNames().size());
  std::vector<std::array<double, 3>> antennaPos(info().antennaPos().size());
  for (unsigned int i = 0; i < info().antennaNames().size(); ++i) {
    allAntennaNames[i] = info().antennaNames()[i];
    casacore::Quantum<casacore::Vector<double>> pos =
        info().antennaPos()[i].get("m");
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  h5parm.AddAntennas(allAntennaNames, antennaPos);

  std::vector<std::pair<double, double>> pointingPosition(1);
  casacore::MDirection phasecenter = info().phaseCenter();
  pointingPosition[0].first = phasecenter.getValue().get()[0];
  pointingPosition[0].second = phasecenter.getValue().get()[1];
  std::vector<string> pointingName(1, "POINTING");

  h5parm.AddSources(pointingName, pointingPosition);

  unsigned int nPol;
  std::vector<string> polarizations;
  if (scalarMode(itsMode)) {
    nPol = 1;
  } else if (diagonalMode(itsMode)) {
    nPol = 2;
    polarizations.push_back("XX");
    polarizations.push_back("YY");
  } else {
    assert(itsMode == CalType::kFullJones);
    polarizations.push_back("XX");
    polarizations.push_back("XY");
    polarizations.push_back("YX");
    polarizations.push_back("YY");
    nPol = 4;
  }

  // Construct time axis
  unsigned int nSolTimes = (info().ntime() + itsSolInt - 1) / itsSolInt;
  std::vector<double> solTimes(nSolTimes);
  assert(nSolTimes == itsSols.size());
  double starttime = info().startTime();
  for (unsigned int t = 0; t < nSolTimes; ++t) {
    solTimes[t] = starttime + (t + 0.5) * info().timeInterval() * itsSolInt;
  }

  // Construct frequency axis
  unsigned int nSolFreqs;
  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    nSolFreqs = 1;
  } else {
    nSolFreqs = itsNFreqCells;
  }

  std::vector<AxisInfo> axes;
  axes.push_back(AxisInfo("time", itsSols.size()));
  axes.push_back(AxisInfo("freq", nSolFreqs));
  axes.push_back(AxisInfo("ant", info().antennaUsed().size()));
  if (nPol > 1) {
    axes.push_back(AxisInfo("pol", nPol));
  }

  std::vector<SolTab> soltabs = makeSolTab(h5parm, itsMode, axes);

  std::vector<std::string> antennaUsedNames;
  for (unsigned int st = 0; st < info().antennaUsed().size(); ++st) {
    antennaUsedNames.push_back(info().antennaNames()[info().antennaUsed()[st]]);
  }

  std::vector<SolTab>::iterator soltabiter = soltabs.begin();
  for (; soltabiter != soltabs.end(); ++soltabiter) {
    (*soltabiter).SetAntennas(antennaUsedNames);
    if (nPol > 1) {
      (*soltabiter).SetPolarizations(polarizations);
    }
    if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
      // Set channel to frequency of middle channel
      // TODO: fix this for nchan
      std::vector<double> oneFreq(1);
      oneFreq[0] = info().chanFreqs()[info().nchan() / 2];
      (*soltabiter).SetFreqs(oneFreq);
    } else {
      (*soltabiter).SetFreqs(itsFreqData);
    }
    (*soltabiter).SetTimes(solTimes);
  }

  // Put solutions in a contiguous piece of memory
  string historyString = "CREATE by " + DP3Version::AsString() + "\n" +
                         "step " + itsName + " in parset: \n" + itsParsetString;

  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    std::vector<double> tecsols(nSolFreqs * antennaUsedNames.size() *
                                nSolTimes * nPol);
    std::vector<double> weights(
        nSolFreqs * antennaUsedNames.size() * nSolTimes * nPol, 1.);
    std::vector<double> phasesols;
    if (itsMode == CalType::kTecAndPhase) {
      phasesols.resize(nSolFreqs * antennaUsedNames.size() * nSolTimes * nPol);
    }
    size_t i = 0;
    for (unsigned int time = 0; time < nSolTimes; ++time) {
      for (unsigned int freqCell = 0; freqCell < nSolFreqs; ++freqCell) {
        for (unsigned int ant = 0; ant < info().antennaUsed().size(); ++ant) {
          for (unsigned int pol = 0; pol < nPol; ++pol) {
            assert(!itsTECSols[time].empty());
            tecsols[i] = itsTECSols[time](0, ant) / 8.44797245e9;
            if (!std::isfinite(tecsols[i])) {
              weights[i] = 0.;
            }
            if (itsMode == CalType::kTecAndPhase) {
              phasesols[i] = -itsTECSols[time](0, ant);
            }
            ++i;
          }
        }
      }
    }
    soltabs[0].SetValues(tecsols, weights, historyString);
    if (itsMode == CalType::kTecAndPhase) {
      soltabs[1].SetValues(phasesols, weights, historyString);
    }
  } else {
    std::vector<casacore::DComplex> sols(nSolFreqs * antennaUsedNames.size() *
                                         nSolTimes * nPol);
    std::vector<double> weights(
        nSolFreqs * antennaUsedNames.size() * nSolTimes * nPol, 1.);
    size_t i = 0;
    for (unsigned int time = 0; time < nSolTimes; ++time) {
      for (unsigned int freqCell = 0; freqCell < nSolFreqs; ++freqCell) {
        for (unsigned int ant = 0; ant < info().antennaUsed().size(); ++ant) {
          for (unsigned int pol = 0; pol < nPol; ++pol) {
            assert(!itsSols[time].empty());
            sols[i] = itsSols[time](pol, ant, freqCell);
            if (!std::isfinite(sols[i].real())) {
              weights[i] = 0.;
            }
            ++i;
          }
        }
      }
    }

    if (itsMode != CalType::kDiagonalAmplitude) {
      soltabs[0].SetComplexValues(sols, weights, false, historyString);
    } else {
      soltabs[0].SetComplexValues(sols, weights, true, historyString);
    }
    if (soltabs.size() > 1) {
      // Also write amplitudes
      soltabs[1].SetComplexValues(sols, weights, true, historyString);
    }
  }

  itsTimerWrite.stop();
  itsTimer.stop();
}

std::vector<SolTab> GainCal::makeSolTab(H5Parm& h5parm, CalType caltype,
                                        std::vector<AxisInfo>& axes) {
  unsigned int numsols = 1;
  // For [scalar]complexgain, store two soltabs: phase and amplitude
  if (caltype == CalType::kDiagonal || caltype == CalType::kScalar ||
      caltype == CalType::kTecAndPhase || caltype == CalType::kFullJones) {
    numsols = 2;
  }
  std::vector<SolTab> soltabs;
  for (unsigned int solnum = 0; solnum < numsols; ++solnum) {
    string solTabName;
    SolTab soltab;
    switch (caltype) {
      case CalType::kScalarPhase:
      case CalType::kDiagonalPhase:
        solTabName = "phase000";
        soltab = h5parm.CreateSolTab(solTabName, "phase", axes);
        break;
      case CalType::kScalar:
      case CalType::kDiagonal:
      case CalType::kFullJones:
        if (solnum == 0) {
          solTabName = "phase000";
          soltab = h5parm.CreateSolTab(solTabName, "phase", axes);
        } else {
          solTabName = "amplitude000";
          soltab = h5parm.CreateSolTab(solTabName, "amplitude", axes);
        }
        break;
      case CalType::kScalarAmplitude:
      case CalType::kDiagonalAmplitude:
        solTabName = "amplitude000";
        soltab = h5parm.CreateSolTab(solTabName, "amplitude", axes);
        break;
      case CalType::kTec:
      case CalType::kTecAndPhase:
        if (solnum == 0) {
          solTabName = "tec000";
          soltab = h5parm.CreateSolTab(solTabName, "tec", axes);
        } else {
          solTabName = "phase000";
          soltab = h5parm.CreateSolTab(solTabName, "phase", axes);
        }
        break;
      default:
        throw Exception("Unhandled mode in writing H5Parm output: " +
                        ToString(caltype));
    }
    soltabs.push_back(soltab);
  }
  return soltabs;
}

void GainCal::writeSolutionsParmDB(double startTime) {
  itsTimer.start();
  itsTimerWrite.start();

  // Open the ParmDB at the first write.
  // In that way the instrumentmodel ParmDB can be in the MS directory.
  if (!itsParmDB) {
    initParmDB();
  }  // End initialization of parmdb

  unsigned int ntime = itsSols.size();
  unsigned int nchan, nfreqs;
  if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
    nfreqs = 1;
    nchan = info().nchan();
  } else {
    nfreqs = itsNFreqCells;
    nchan = itsNChan;
  }

  // Construct solution grid for the current chunk
  double freqWidth = getInfo().chanWidths()[0];
  if (getInfo().chanFreqs().size() >
      1) {  // Handle data with evenly spaced gaps between channels
    freqWidth = info().chanFreqs()[1] - info().chanFreqs()[0];
  }

  // Get end time of the current chunk. For the last chunk, this
  // is chopped off at the end of the MS (only if solint > 1)
  double endTime =
      std::min(startTime + ntime * info().timeInterval() * itsSolInt,
               info().startTime() + info().ntime() * info().timeInterval());

  // Make time axis (can be non regular for last chunk if solint > 1)
  std::vector<double> lowtimes(ntime), hightimes(ntime);
  for (unsigned int t = 0; t < ntime; ++t) {
    lowtimes[t] = startTime + info().timeInterval() * itsSolInt * t;
    hightimes[t] = std::min(
        startTime + info().timeInterval() * itsSolInt * (t + 1), endTime);
  }
  parmdb::Axis::ShPtr timeAxis = parmdb::Axis::makeAxis(lowtimes, hightimes);

  parmdb::Axis::ShPtr freqAxis(new parmdb::RegularAxis(
      getInfo().chanFreqs()[0] - freqWidth * 0.5, freqWidth * nchan, nfreqs));
  parmdb::Grid solGrid(freqAxis, timeAxis);

  // Construct domain grid for the current chunk
  parmdb::Axis::ShPtr tdomAxis(
      new parmdb::RegularAxis(startTime, endTime - startTime, 1));
  parmdb::Axis::ShPtr fdomAxis(
      new parmdb::RegularAxis(info().chanFreqs()[0] - freqWidth * 0.5,
                              freqWidth * getInfo().chanFreqs().size(), 1));
  parmdb::Grid domainGrid(fdomAxis, tdomAxis);

  // Write the solutions per parameter.
  const char* str0101[] = {"0:0:", "0:1:", "1:0:", "1:1:"};
  const char* strri[] = {"Real:", "Imag:"};
  Matrix<double> values(nfreqs, ntime);

  casacore::DComplex sol;

  unsigned int nSt = info().antennaUsed().size();

  for (size_t st = 0; st < nSt; ++st) {
    // Do not write NaN solutions for stations that were not used
    if (!itsAntennaUsed[info().antennaUsed()[st]]) {
      // itsAntennaUsed is indexed with real antenna numbers, so antennaUsed()
      // is needed
      continue;
    }
    for (int pol = 0; pol < 4; ++pol) {  // For 0101
      if (scalarMode(itsMode) && pol > 0) {
        continue;
      } else if (diagonalMode(itsMode) && (pol == 1 || pol == 2)) {
        continue;
      }
      int realimmax = 2;  // For tecandphase, this functions as dummy between
                          // tec and commonscalarphase
      if (itsMode == CalType::kDiagonalPhase ||
          itsMode == CalType::kScalarPhase ||
          itsMode == CalType::kDiagonalAmplitude ||
          itsMode == CalType::kScalarAmplitude || itsMode == CalType::kTec) {
        realimmax = 1;
      }
      for (int realim = 0; realim < realimmax;
           ++realim) {  // For real and imaginary
        string name = parmName();

        if (itsMode != CalType::kScalarPhase &&
            itsMode != CalType::kScalarAmplitude) {
          name += str0101[pol];
          if (itsMode == CalType::kDiagonalPhase) {
            name = name + "Phase:";
          } else if (itsMode == CalType::kDiagonalAmplitude) {
            name = name + "Ampl:";
          } else {
            name = name + strri[realim];
          }
        }
        if (itsMode == CalType::kTec || itsMode == CalType::kTecAndPhase) {
          if (realim == 0) {
            name = "TEC:";
          } else {
            name = "CommonScalarPhase:";
          }
        }

        name += info().antennaNames()[info().antennaUsed()[st]];

        // Collect its solutions for all times and frequency cells in a single
        // array.
        for (unsigned int ts = 0; ts < ntime; ++ts) {
          for (unsigned int freqCell = 0; freqCell < nfreqs; ++freqCell) {
            if (itsMode == CalType::kFullJones) {
              if (realim == 0) {
                values(freqCell, ts) = real(itsSols[ts](pol, st, freqCell));
              } else {
                values(freqCell, ts) = imag(itsSols[ts](pol, st, freqCell));
              }
            } else if (itsMode == CalType::kDiagonal) {
              if (realim == 0) {
                values(freqCell, ts) = real(itsSols[ts](pol / 3, st, freqCell));
              } else {
                values(freqCell, ts) = imag(itsSols[ts](pol / 3, st, freqCell));
              }
            } else if (itsMode == CalType::kScalarPhase ||
                       itsMode == CalType::kDiagonalPhase) {
              values(freqCell, ts) = arg(itsSols[ts](pol / 3, st, freqCell));
            } else if (itsMode == CalType::kScalarAmplitude ||
                       itsMode == CalType::kDiagonalAmplitude) {
              values(freqCell, ts) = abs(itsSols[ts](pol / 3, st, freqCell));
            } else if (itsMode == CalType::kTec ||
                       itsMode == CalType::kTecAndPhase) {
              if (realim == 0) {
                values(freqCell, ts) =
                    itsTECSols[ts](realim, st) / 8.44797245e9;
              } else {
                values(freqCell, ts) = -itsTECSols[ts](
                    realim, st);  // TODO: why is there a minus here?
              }
            } else {
              throw Exception("Unhandled mode");
            }
          }
        }
        parmdb::ParmValue::ShPtr pv(new parmdb::ParmValue());
        pv->setScalars(solGrid, values);

        parmdb::ParmValueSet pvs(domainGrid,
                                 std::vector<parmdb::ParmValue::ShPtr>(1, pv));
        std::map<std::string, int>::const_iterator pit =
            itsParmIdMap.find(name);

        if (pit == itsParmIdMap.end()) {
          // First time, so a new nameId will be set.
          // Check if the name was defined in the parmdb previously
          int nameId = itsParmDB->getNameId(name);
          itsParmDB->putValues(name, nameId, pvs);
          itsParmIdMap[name] = nameId;
        } else {
          // Parm has been put before.
          int nameId = pit->second;
          itsParmDB->putValues(name, nameId, pvs);
        }
      }
    }
  }

  itsTimerWrite.stop();
  itsTimer.stop();
}

void GainCal::finish() {
  itsTimer.start();

  // Solve remaining time slots if any
  if (itsStepInSolInt != 0) {
    calibrate();

    if (itsApplySolution) {
      Cube<casacore::Complex> invsol = invertSol(itsSols.back());
      for (unsigned int stepInSolInt = 0; stepInSolInt < itsStepInSolInt;
           stepInSolInt++) {
        applySolution(itsBuf[stepInSolInt], invsol);
        getNextStep()->process(itsBuf[stepInSolInt]);
      }
    }
  }

  itsTimer.stop();

  if (!itsSols.empty()) {
    if (itsUseH5Parm) {
      writeSolutionsH5Parm(itsChunkStartTime);
    } else {
      writeSolutionsParmDB(itsChunkStartTime);
    }
    if (itsDebugLevel > 0) {
      H5::H5File hdf5file = H5::H5File("debug.h5", H5F_ACC_TRUNC);
      std::vector<hsize_t> dims(6);
      for (unsigned int i = 0; i < 6; ++i) {
        dims[i] = itsAllSolutions.shape()[5 - i];
      }
      H5::DataSpace dataspace(dims.size(), &(dims[0]), NULL);
      H5::CompType complex_data_type(sizeof(casacore::DComplex));
      complex_data_type.insertMember("r", 0, H5::PredType::IEEE_F64LE);
      complex_data_type.insertMember("i", sizeof(double),
                                     H5::PredType::IEEE_F64LE);
      H5::DataSet dataset =
          hdf5file.createDataSet("val", complex_data_type, dataspace);
      dataset.write(itsAllSolutions.data(), complex_data_type);
      hdf5file.close();
    }
  }

  // Let the next steps finish.
  getNextStep()->finish();
}

}  // namespace steps
}  // namespace dp3
