// DDECal.cc: DPPP step class to do a direction dependent gain calibration
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema & André Offringa

#include "DDECal.h"

#include <Version.h>

#include "../base/CalType.h"
#include "../base/DP3.h"
#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/DPLogger.h"
#include "../base/Simulate.h"
#include "../base/SourceDBUtil.h"

#include "../steps/IDGPredict.h"
#include "../steps/MSReader.h"
#include "../steps/ColumnReader.h"

#include "../ddecal/SolverFactory.h"
#ifdef HAVE_ARMADILLO
#include "../ddecal/constraints/ScreenConstraint.h"
#endif
#include "../ddecal/constraints/SmoothnessConstraint.h"
#include "../ddecal/gain_solvers/SolverBuffer.h"
#include "../ddecal/linear_solvers/LLSSolver.h"

#include <schaapcommon/facets/facet.h>

#include <aocommon/matrix2x2.h>

#include "../parmdb/ParmDB.h"
#include "../parmdb/ParmValue.h"
#include "../parmdb/SourceDB.h"

#include "../common/Memory.h"
#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"
#include "../common/StringTools.h"

#include <aocommon/threadpool.h>

#include <boost/make_unique.hpp>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/OS/File.h>

#include <algorithm>
#include <cassert>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

using aocommon::FitsReader;

using schaapcommon::facets::Facet;
using schaapcommon::h5parm::SolTab;

using dp3::base::CalType;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::ddecal::LLSSolver;
using dp3::ddecal::LLSSolverType;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

DDECal::DDECal(InputStep* input, const common::ParameterSet& parset,
               const string& prefix)
    : itsInput(*input),
      itsSettings(parset, prefix),
      itsAvgTime(0),
      itsSols(),
      itsSolutionWriter(itsSettings.h5parm_name),
      itsTimeStep(0),
      itsRequestedSolInt(itsSettings.solution_interval),
      itsSolIntCount(1),
      itsNSolInts(0),
      itsBufferedSolInts(0),
      itsNChan(itsSettings.n_channels),
      itsUVWFlagStep(input, parset, prefix),
      itsSolver(ddecal::CreateRegularSolver(itsSettings, parset, prefix)),
      itsStatStream() {
  if (!itsSettings.stat_filename.empty()) {
    itsStatStream =
        boost::make_unique<std::ofstream>(itsSettings.stat_filename);
  }

  // Initialize steps
  initializeColumnReaders(parset, prefix);
  initializeIDG(parset, prefix);
  initializePredictSteps(parset, prefix);

  if (itsDirections.size() == 0)
    throw std::runtime_error(
        "DDECal initialized with 0 directions: something is wrong with your "
        "parset or your sourcedb");
}

DDECal::~DDECal() {}

void DDECal::initializeColumnReaders(const common::ParameterSet& parset,
                                     const string& prefix) {
  for (const std::string& col : itsSettings.model_data_columns) {
    if (itsSettings.model_data_columns.size() == 1) {
      itsDirections.emplace_back(1, "pointing");
    } else {
      itsDirections.emplace_back(1, col);
    }
    itsSteps.push_back(
        std::make_shared<ColumnReader>(itsInput, parset, prefix, col));
    setModelNextSteps(*itsSteps.back(), col, parset, prefix);
  }
}

void DDECal::initializeIDG(const common::ParameterSet& parset,
                           const string& prefix) {
  // TODO it would be nicer to get a new method in the DS9 reader to first get
  // names of directions and pass that to idgpredict. It will then read it
  // itself instead of DDECal having to do everything. It is better to do it all
  // in IDGPredict, so we can also make it
  if (itsSettings.idg_region_filename.empty() &&
      itsSettings.idg_image_filenames.empty()) {
    return;
  }

  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<float>>>
      readers = IDGPredict::GetReaders(itsSettings.idg_image_filenames);
  std::vector<Facet> facets = IDGPredict::GetFacets(
      itsSettings.idg_region_filename, readers.first.front());

  for (size_t i = 0; i < facets.size(); ++i) {
    std::string dir_name = "dir" + std::to_string(i);
    if (!facets[i].Direction().empty()) {
      dir_name = facets[i].Direction();
    }
    itsDirections.emplace_back(1, dir_name);

    itsSteps.push_back(std::make_shared<IDGPredict>(
        itsInput, parset, prefix, readers, std::vector<Facet>{facets[i]}));
    setModelNextSteps(*itsSteps.back(), facets[i].Direction(), parset, prefix);
  }
}

void DDECal::initializePredictSteps(const common::ParameterSet& parset,
                                    const string& prefix) {
  size_t start = itsDirections.size();

  // Default directions are all patches
  if (itsSettings.directions.empty() && !itsSettings.source_db.empty()) {
    parmdb::SourceDB sourceDB(parmdb::ParmDBMeta("", itsSettings.source_db),
                              true, false);
    std::vector<string> patchNames =
        base::makePatchList(sourceDB, std::vector<string>());
    for (const string& patch : patchNames) {
      itsDirections.emplace_back(1, patch);
    }
  } else {
    for (const string& direction : itsSettings.directions) {
      common::ParameterValue dirStr(direction);
      itsDirections.emplace_back(dirStr.getStringVector());
    }
  }

  for (size_t dir = start; dir < itsDirections.size(); ++dir) {
    itsSteps.push_back(std::make_shared<Predict>(itsInput, parset, prefix,
                                                 itsDirections[dir]));
    setModelNextSteps(*itsSteps.back(), itsDirections[dir][0], parset, prefix);
  }
}

void DDECal::setModelNextSteps(Step& step, const std::string& direction,
                               const common::ParameterSet& parset,
                               const string& prefix) const {
  std::string step_names_key = prefix + "modelnextsteps." + direction;
  if (!parset.isDefined(step_names_key)) {
    step_names_key = prefix + "modelnextsteps";  // Fall back setting.
  }

  if (parset.isDefined(step_names_key)) {
    Step::ShPtr first_step =
        base::DP3::makeStepsFromParset(parset, "", step_names_key, itsInput,
                                       false, steps::Step::MsType::kRegular);

    if (first_step) {
      step.setNextStep(first_step);
    }
  }
}

void DDECal::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);
  info().setNeedVisData();
  if (itsSettings.subtract) info().setWriteData();

  itsUVWFlagStep.updateInfo(infoIn);
  itsThreadPool =
      boost::make_unique<aocommon::ThreadPool>(getInfo().nThreads());

  // Update info for substeps and set other required parameters
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsSteps[dir]->setInfo(infoIn);

    if (auto s = std::dynamic_pointer_cast<Predict>(itsSteps[dir])) {
      s->SetThreadData(*itsThreadPool, &itsMeasuresMutex);
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(itsSteps[dir])) {
      itsSolIntCount =
          std::max(itsSolIntCount,
                   s->GetBufferSize() / itsSteps.size() / itsRequestedSolInt);
      // We increment by one so the IDGPredict will not flush in its process
      s->SetBufferSize(itsRequestedSolInt * itsSolIntCount + 1);
    } else if (!std::dynamic_pointer_cast<ColumnReader>(itsSteps[dir])) {
      throw std::runtime_error("DDECal received an invalid first model step");
    }
  }

  itsSolver->SetNThreads(getInfo().nThreads());

  if (itsRequestedSolInt == 0) {
    itsRequestedSolInt = info().ntime();
  }

  itsDataResultStep = std::make_shared<ResultStep>();
  itsUVWFlagStep.setNextStep(itsDataResultStep);

  // Each step should be coupled to a resultstep
  itsResultSteps.resize(itsSteps.size());
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsResultSteps[dir] =
        std::make_shared<MultiResultStep>(itsRequestedSolInt * itsSolIntCount);

    // Add the resultstep to the end of the model next steps
    std::shared_ptr<Step> step = itsSteps[dir];
    while (step->getNextStep()) {
      step = step->getNextStep();
    }
    step->setNextStep(itsResultSteps[dir]);
  }

  if (itsNChan == 0 || itsNChan > info().nchan()) {
    itsNChan = info().nchan();
  }

  // Convert from casacore::Vector to std::vector, pass only used antennas to
  // multidirsolver
  std::vector<int> ant1(info().getAnt1().size());
  std::vector<int> ant2(info().getAnt2().size());
  for (unsigned int i = 0; i < ant1.size(); ++i) {
    ant1[i] = info().antennaMap()[info().getAnt1()[i]];
    ant2[i] = info().antennaMap()[info().getAnt2()[i]];
  }

  // Fill antenna info in H5Parm, need to convert from casa types to std types
  // Fill in metadata for all antennas, also those that may be filtered out.
  std::vector<std::string> antennaNames(info().antennaNames().size());
  std::vector<std::array<double, 3>> antennaPos(info().antennaPos().size());
  for (unsigned int i = 0; i < info().antennaNames().size(); ++i) {
    antennaNames[i] = info().antennaNames()[i];
    casacore::Quantum<casacore::Vector<double>> pos =
        info().antennaPos()[i].get("m");
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  itsSolutionWriter.AddAntennas(antennaNames, antennaPos);

  // Prepare positions for only used antennas, to be used for constraints
  std::vector<std::array<double, 3>> usedAntennaPositions(
      info().antennaUsed().size());
  for (size_t i = 0; i < info().antennaUsed().size(); ++i) {
    usedAntennaPositions[i] = antennaPos[info().antennaUsed()[i]];
  }

  size_t nSolTimes =
      (info().ntime() + itsRequestedSolInt - 1) / itsRequestedSolInt;
  size_t nChannelBlocks = info().nchan() / itsNChan;
  itsSols.resize(nSolTimes);
  itsNIter.resize(nSolTimes);
  itsNApproxIter.resize(nSolTimes);
  itsConstraintSols.resize(nSolTimes);
  itsVisInInterval.assign(nChannelBlocks, std::pair<size_t, size_t>(0, 0));

  itsChanBlockStart.resize(nChannelBlocks + 1);
  itsChanBlockFreqs.resize(nChannelBlocks);
  itsChanBlockStart.front() = 0;
  for (size_t chBlock = 0; chBlock != nChannelBlocks; ++chBlock) {
    itsChanBlockStart[chBlock + 1] =
        (chBlock + 1) * info().nchan() / nChannelBlocks;
    const size_t blockSize =
        itsChanBlockStart[chBlock + 1] - itsChanBlockStart[chBlock];
    const double* freqStart =
        info().chanFreqs().data() + itsChanBlockStart[chBlock];
    const double meanFreq =
        std::accumulate(freqStart, freqStart + blockSize, 0.0) / blockSize;
    itsChanBlockFreqs[chBlock] = meanFreq;
  }

  itsWeightsPerAntenna.assign(
      itsChanBlockFreqs.size() * info().antennaUsed().size(), 0.0);

  std::vector<ddecal::SolverBase*> solvers = itsSolver->ConstraintSolvers();
  for (ddecal::SolverBase* solver : solvers) {
    for (const std::unique_ptr<ddecal::Constraint>& constraint :
         solver->GetConstraints()) {
      // Initialize the constraint with some common metadata
      constraint->Initialize(info().antennaUsed().size(), itsDirections.size(),
                             itsChanBlockFreqs);

      // Different constraints need different information. Determine if the
      // constraint is of a type that needs more information, and if so
      // initialize the constraint.
      ddecal::AntennaConstraint* antConstraint =
          dynamic_cast<ddecal::AntennaConstraint*>(constraint.get());
      if (antConstraint != nullptr) {
        if (itsSettings.antenna_constraint.empty()) {
          // Set the antenna constraint to all stations within certain distance
          // specified by 'coreconstraint' parameter.
          // Take the first used station as reference station
          const double refX = usedAntennaPositions[0][0],
                       refY = usedAntennaPositions[0][1],
                       refZ = usedAntennaPositions[0][2];
          std::vector<std::set<size_t>> antConstraintList(1);
          std::set<size_t>& coreAntennaIndices = antConstraintList.front();
          const double coreDistSq =
              itsSettings.core_constraint * itsSettings.core_constraint;
          for (size_t ant = 0; ant != usedAntennaPositions.size(); ++ant) {
            const double dx = refX - usedAntennaPositions[ant][0],
                         dy = refY - usedAntennaPositions[ant][1],
                         dz = refZ - usedAntennaPositions[ant][2],
                         distSq = dx * dx + dy * dy + dz * dz;
            if (distSq <= coreDistSq) coreAntennaIndices.insert(ant);
          }
          antConstraint->SetAntennaSets(std::move(antConstraintList));
        } else {
          // Set the antenna constraint to a list of stations indices that
          // are to be kept the same during the solve.
          const casacore::Vector<casacore::String>& antNames =
              info().antennaNames();
          std::vector<std::string> antNamesStl(
              antNames.begin(),
              antNames.end());  // casacore vector doesn't support find properly
          std::vector<std::set<size_t>> constraintList;
          for (const std::set<std::string>& constraintNameSet :
               itsSettings.antenna_constraint) {
            constraintList.emplace_back();
            for (const std::string& constraintName : constraintNameSet) {
              auto iter = std::find(antNamesStl.begin(), antNamesStl.end(),
                                    constraintName);
              if (iter != antNamesStl.end())
                constraintList.back().insert(
                    info().antennaMap()[iter - antNamesStl.begin()]);
            }
            if (constraintList.back().size() <= 1)
              throw std::runtime_error(
                  "Error in antenna constraint: at least two antennas "
                  "expected");
          }
          antConstraint->SetAntennaSets(std::move(constraintList));
        }
      }

#ifdef HAVE_ARMADILLO
      ddecal::ScreenConstraint* screenConstraint =
          dynamic_cast<ddecal::ScreenConstraint*>(constraint.get());
      if (screenConstraint != 0) {
        screenConstraint->setAntennaPositions(usedAntennaPositions);
        screenConstraint->setDirections(GetSourcePositions());
        screenConstraint->initPiercePoints();
        const size_t ref_antenna = 0;
        const double refX = usedAntennaPositions[ref_antenna][0],
                     refY = usedAntennaPositions[ref_antenna][1],
                     refZ = usedAntennaPositions[ref_antenna][2];
        std::vector<size_t> coreAntennaIndices;
        std::vector<size_t> otherAntennaIndices;
        const double coreDistSq = itsSettings.screen_core_constraint *
                                  itsSettings.screen_core_constraint;
        for (size_t ant = 0; ant != usedAntennaPositions.size(); ++ant) {
          const double dx = refX - usedAntennaPositions[ant][0],
                       dy = refY - usedAntennaPositions[ant][1],
                       dz = refZ - usedAntennaPositions[ant][2],
                       distSq = dx * dx + dy * dy + dz * dz;
          if (distSq <= coreDistSq)
            coreAntennaIndices.emplace_back(info().antennaMap()[ant]);
          else
            otherAntennaIndices.emplace_back(info().antennaMap()[ant]);
        }
        screenConstraint->setCoreAntennas(coreAntennaIndices);
        screenConstraint->setOtherAntennas(otherAntennaIndices);
      }
#endif

      ddecal::SmoothnessConstraint* sConstraint =
          dynamic_cast<ddecal::SmoothnessConstraint*>(constraint.get());
      if (sConstraint) {
        std::vector<double> distanceFactors;
        // If no smoothness reference distance is specified, the smoothing is
        // made independent of the distance
        if (itsSettings.smoothness_ref_distance == 0.0) {
          distanceFactors.assign(usedAntennaPositions.size(), 1.0);
        } else {
          // Make a list of factors such that more distant antennas apply a
          // smaller smoothing kernel.
          distanceFactors.reserve(usedAntennaPositions.size());
          for (size_t i = 1; i != usedAntennaPositions.size(); ++i) {
            const double dx =
                usedAntennaPositions[0][0] - usedAntennaPositions[i][0];
            const double dy =
                usedAntennaPositions[0][1] - usedAntennaPositions[i][1];
            const double dz =
                usedAntennaPositions[0][2] - usedAntennaPositions[i][2];
            const double factor = itsSettings.smoothness_ref_distance /
                                  std::sqrt(dx * dx + dy * dy + dz * dz);
            distanceFactors.push_back(factor);
            // For antenna 0, the distance of antenna 1 is used:
            if (i == 1) distanceFactors.push_back(factor);
          }
        }
        sConstraint->SetDistanceFactors(std::move(distanceFactors));
      }
    }
  }

  size_t nSt = info().antennaUsed().size();
  // Give renumbered antennas to solver
  itsSolver->Initialize(nSt, itsDirections.size(), info().nchan(),
                        nChannelBlocks, ant1, ant2);

  for (size_t i = 0; i < nSolTimes; ++i) {
    itsSols[i].resize(nChannelBlocks);
  }
}

void DDECal::show(std::ostream& os) const {
  os << "DDECal " << itsSettings.name << '\n'
     << "  mode (constraints):  " << ToString(itsSettings.mode) << '\n'
     << "  algorithm:           "
     << ddecal::ToString(itsSettings.solver_algorithm) << '\n'
     << "  H5Parm:              " << itsSettings.h5parm_name << '\n'
     << "  solint:              " << itsRequestedSolInt << '\n'
     << "  nchan:               " << itsNChan << '\n'
     << "  directions:          " << itsDirections << '\n';
  if (itsSettings.min_vis_ratio != 0.0) {
    os << "  min visib. ratio:    " << itsSettings.min_vis_ratio << '\n';
  }
  os << "  tolerance:           " << itsSolver->GetAccuracy() << '\n'
     << "  max iter:            " << itsSolver->GetMaxIterations() << '\n'
     << "  flag unconverged:    " << std::boolalpha
     << itsSettings.flag_unconverged << '\n'
     << "     diverged only:    " << std::boolalpha
     << itsSettings.flag_diverged_only << '\n'
     << "  propagate solutions: " << std::boolalpha
     << itsSettings.propagate_solutions << '\n'
     << "       converged only: " << std::boolalpha
     << itsSettings.propagate_converged_only << '\n'
     << "  detect stalling:     " << std::boolalpha
     << itsSolver->GetDetectStalling() << '\n'
     << "  step size:           " << itsSolver->GetStepSize() << '\n';
  if (!itsSettings.antenna_constraint.empty())
    os << "  antennaconstraint:   " << itsSettings.antenna_constraint << '\n';
  if (itsSettings.core_constraint != 0.0)
    os << "  coreconstraint:      " << itsSettings.core_constraint << '\n';
  if (itsSettings.smoothness_constraint != 0.0)
    os << "  smoothnessconstraint:" << itsSettings.smoothness_constraint
       << '\n';
  os << "  approximate fitter:  " << itsSettings.approximate_tec << '\n'
     << "  only predict:        " << itsSettings.only_predict << '\n'
     << "  subtract model:      " << itsSettings.subtract << '\n';
  for (unsigned int i = 0; i < itsSteps.size(); ++i) {
    std::shared_ptr<Step> step = itsSteps[i];
    os << "Model steps for direction " << itsDirections[i][0] << '\n';
    do {
      step->show(os);
    } while (step = step->getNextStep());
    os << '\n';
  }
  itsUVWFlagStep.show(os);
}

void DDECal::showTimings(std::ostream& os, double duration) const {
  double totaltime = itsTimer.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " DDECal " << itsSettings.name << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerPredict.getElapsed(), totaltime);
  os << " of it spent in predict" << '\n';

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerSolve.getElapsed(), totaltime);
  os << " of it spent in estimating gains and computing residuals" << '\n';

  itsSolver->GetTimings(os, itsTimerSolve.getElapsed());

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerWrite.getElapsed(), totaltime);
  os << " of it spent in writing gain solutions to disk" << '\n';

  os << "          ";
  os << "Substeps taken:" << '\n';
  for (auto& step : itsSteps) {
    os << "          ";
    step->showTimings(os, duration);
  }

  os << "Iterations taken: [";
  for (size_t i = 0; i < itsNIter.size() - 1; ++i) {
    os << itsNIter[i];
    if (itsNApproxIter[i] != 0) os << '|' << itsNApproxIter[i];
    os << ",";
  }
  os << itsNIter[itsNIter.size() - 1];
  if (itsNApproxIter[itsNIter.size() - 1] != 0)
    os << '|' << itsNApproxIter[itsNIter.size() - 1];
  os << "]" << '\n';
}

void DDECal::InitializeScalarOrDiagonalSolutions(size_t bufferIndex) {
  if (itsSolIntBuffers[bufferIndex].NSolution() > 0 &&
      itsSettings.propagate_solutions) {
    if (itsNIter[itsSolIntBuffers[bufferIndex].NSolution() - 1] >
            itsSolver->GetMaxIterations() &&
        itsSettings.propagate_converged_only) {
      // initialize solutions with 1.
      const size_t n = itsDirections.size() * info().antennaUsed().size() *
                       itsSolver->NSolutionPolarizations();
      for (std::vector<casacore::DComplex>& solvec :
           itsSols[itsSolIntBuffers[bufferIndex].NSolution()]) {
        solvec.assign(n, 1.0);
      }
    } else {
      // initialize solutions with those of the previous step
      itsSols[itsSolIntBuffers[bufferIndex].NSolution()] =
          itsSols[itsSolIntBuffers[bufferIndex].NSolution() - 1];
    }
  } else {
    // initialize solutions with 1.
    const size_t n = itsDirections.size() * info().antennaUsed().size() *
                     itsSolver->NSolutionPolarizations();
    for (std::vector<casacore::DComplex>& solvec :
         itsSols[itsSolIntBuffers[bufferIndex].NSolution()]) {
      solvec.assign(n, 1.0);
    }
  }
}

void DDECal::initializeFullMatrixSolutions(size_t bufferIndex) {
  if (itsSolIntBuffers[bufferIndex].NSolution() > 0 &&
      itsSettings.propagate_solutions) {
    if (itsNIter[itsSolIntBuffers[bufferIndex].NSolution() - 1] >
            itsSolver->GetMaxIterations() &&
        itsSettings.propagate_converged_only) {
      // initialize solutions with unity matrix [1 0 ; 0 1].
      size_t n = itsDirections.size() * info().antennaUsed().size();
      for (std::vector<casacore::DComplex>& solvec :
           itsSols[itsSolIntBuffers[bufferIndex].NSolution()]) {
        solvec.resize(n * 4);
        for (size_t i = 0; i != n; ++i) {
          solvec[i * 4 + 0] = 1.0;
          solvec[i * 4 + 1] = 0.0;
          solvec[i * 4 + 2] = 0.0;
          solvec[i * 4 + 3] = 1.0;
        }
      }
    } else {
      // initialize solutions with those of the previous step
      itsSols[itsSolIntBuffers[bufferIndex].NSolution()] =
          itsSols[itsSolIntBuffers[bufferIndex].NSolution() - 1];
    }
  } else {
    // initialize solutions with unity matrix [1 0 ; 0 1].
    size_t n = itsDirections.size() * info().antennaUsed().size();
    for (std::vector<casacore::DComplex>& solvec :
         itsSols[itsSolIntBuffers[bufferIndex].NSolution()]) {
      solvec.resize(n * 4);
      for (size_t i = 0; i != n; ++i) {
        solvec[i * 4 + 0] = 1.0;
        solvec[i * 4 + 1] = 0.0;
        solvec[i * 4 + 2] = 0.0;
        solvec[i * 4 + 3] = 1.0;
      }
    }
  }
}

void DDECal::flagChannelBlock(size_t cbIndex, size_t bufferIndex) {
  const size_t nBl = info().nbaselines();
  const size_t nChanBlocks = itsChanBlockFreqs.size();
  // Set the antenna-based weights to zero
  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t ant1 = info().antennaMap()[info().getAnt1()[bl]];
    size_t ant2 = info().antennaMap()[info().getAnt2()[bl]];
    for (size_t ch = itsChanBlockStart[cbIndex];
         ch != itsChanBlockStart[cbIndex + 1]; ++ch) {
      itsWeightsPerAntenna[ant1 * nChanBlocks + cbIndex] = 0.0;
      itsWeightsPerAntenna[ant2 * nChanBlocks + cbIndex] = 0.0;
    }
  }
  // Set the visibility weights to zero
  for (DPBuffer& buffer : itsSolIntBuffers[bufferIndex].DataBuffers()) {
    for (size_t bl = 0; bl < nBl; ++bl) {
      float* begin = &buffer.getWeights()(0, itsChanBlockStart[cbIndex], bl);
      float* end = &buffer.getWeights()(0, itsChanBlockStart[cbIndex + 1], bl);
      std::fill(begin, end, 0.0f);
    }
  }
}

void DDECal::checkMinimumVisibilities(size_t bufferIndex) {
  for (size_t cb = 0; cb != itsChanBlockFreqs.size(); ++cb) {
    double fraction =
        double(itsVisInInterval[cb].first) / itsVisInInterval[cb].second;
    if (fraction < itsSettings.min_vis_ratio) flagChannelBlock(cb, bufferIndex);
  }
}

void DDECal::doSolve() {
  std::vector<std::vector<std::vector<DPBuffer>>> model_buffers(
      itsSolIntBuffers.size());
  for (size_t i = 0; i < itsSolIntBuffers.size(); ++i) {
    model_buffers[i].resize(itsSolIntBuffers[i].Size());
    for (std::vector<DPBuffer>& dir_buffers : model_buffers[i]) {
      dir_buffers.reserve(itsDirections.size());
    }
  }

  for (size_t dir = 0; dir < itsDirections.size(); ++dir) {
    if (auto s = dynamic_cast<IDGPredict*>(itsSteps[dir].get())) {
      itsTimerPredict.start();
      s->flush();
      itsTimerPredict.stop();
    }
    for (size_t i = 0; i < itsResultSteps[dir]->size(); ++i) {
      const size_t sol_int = i / itsRequestedSolInt;
      const size_t timestep = i % itsRequestedSolInt;
      model_buffers[sol_int][timestep].emplace_back(
          std::move(itsResultSteps[dir]->get()[i]));
    }
  }

  std::vector<ddecal::SolverBase*> solvers = itsSolver->ConstraintSolvers();
  // Declare solver_buffer outside the loop, so it can reuse its memory.
  ddecal::SolverBuffer solver_buffer;

  for (size_t i = 0; i < itsSolIntBuffers.size(); ++i) {
    // When the model data is subtracted after calibration, the model data
    // needs to be stored before solving, because the solver modifies it.
    // This is done conditionally to prevent using memory when it is
    // not required (the model data can be large).
    if (itsSettings.subtract || itsSettings.only_predict) {
      storeModelData(model_buffers[i]);
    }

    solver_buffer.AssignAndWeight(itsSolIntBuffers[i].DataBuffers(),
                                  std::move(model_buffers[i]));

    ddecal::RegularSolverBase::SolveResult solveResult;
    if (!itsSettings.only_predict) {
      checkMinimumVisibilities(i);

      for (ddecal::SolverBase* solver : solvers) {
        for (const std::unique_ptr<ddecal::Constraint>& constraint :
             solver->GetConstraints()) {
          constraint->SetWeights(itsWeightsPerAntenna);
        }
      }

      if (itsSolver->NSolutionPolarizations() == 4)
        initializeFullMatrixSolutions(i);
      else
        InitializeScalarOrDiagonalSolutions(i);

      itsTimerSolve.start();

      solveResult = itsSolver->Solve(
          solver_buffer, itsSols[itsSolIntBuffers[i].NSolution()],
          itsAvgTime / itsRequestedSolInt, itsStatStream.get());

      itsTimerSolve.stop();

      itsNIter[itsSolIntBuffers[i].NSolution()] = solveResult.iterations;
      itsNApproxIter[itsSolIntBuffers[i].NSolution()] =
          solveResult.constraint_iterations;
    }

    if (itsSettings.subtract || itsSettings.only_predict) {
      subtractCorrectedModel(i);
    }

    // Check for nonconvergence and flag if desired. Unconverged solutions are
    // identified by the number of iterations being one more than the max
    // allowed number
    if (solveResult.iterations > itsSolver->GetMaxIterations() &&
        itsSettings.flag_unconverged) {
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          if (itsSettings.flag_diverged_only) {
            // Set weights with negative values (indicating unconverged
            // solutions that diverged) to zero (all other unconverged
            // solutions remain unflagged)
            for (double& weight : result.weights) {
              if (weight < 0.) weight = 0.;
            }
          } else {
            // Set all weights to zero
            result.weights.assign(result.weights.size(), 0.);
          }
        }
      }
    } else {
      // Set any negative weights (indicating unconverged solutions that
      // diverged) to one (all other unconverged solutions are unflagged
      // already)
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          for (double& weight : result.weights) {
            if (weight < 0.0) weight = 1.0;
          }
        }
      }
    }

    // Store constraint solutions if any constaint has a non-empty result
    bool someConstraintHasResult = false;
    for (const auto& constraint_results : solveResult.results) {
      if (!constraint_results.empty()) {
        someConstraintHasResult = true;
        break;
      }
    }
    if (someConstraintHasResult) {
      itsConstraintSols[itsSolIntBuffers[i].NSolution()] = solveResult.results;
    }
  }

  itsTimer.stop();

  for (size_t i = 0; i < itsSolIntBuffers.size(); ++i) {
    itsSolIntBuffers[i].RestoreFlagsAndWeights();
    for (size_t step = 0; step < itsSolIntBuffers[i].Size(); ++step) {
      // Push data (possibly changed) to next step
      getNextStep()->process(itsSolIntBuffers[i][step]);
    }
  }

  itsTimer.start();
}

bool DDECal::process(const DPBuffer& bufin) {
  itsTimer.start();

  // Create a new solution interval if needed
  if (itsSolIntBuffers.empty() ||
      itsSolIntBuffers.back().Size() == itsRequestedSolInt) {
    itsSolIntBuffers.emplace_back(itsInput, itsNSolInts, itsRequestedSolInt,
                                  itsDirections.size(), itsTimer);
  }

  const size_t currentIntervalIndex = itsSolIntBuffers.back().Size();
  itsSolIntBuffers[itsBufferedSolInts].PushBack(bufin);
  doPrepare(itsSolIntBuffers.back()[currentIntervalIndex], itsBufferedSolInts,
            currentIntervalIndex);

  if (currentIntervalIndex + 1 == itsRequestedSolInt) {
    ++itsBufferedSolInts;
    ++itsNSolInts;
  }

  if (itsBufferedSolInts == itsSolIntCount) {
    doSolve();

    // Clean up, prepare for next iteration
    itsAvgTime = 0;
    itsBufferedSolInts = 0;
    itsVisInInterval.assign(itsVisInInterval.size(),
                            std::pair<size_t, size_t>(0, 0));
    itsWeightsPerAntenna.assign(itsWeightsPerAntenna.size(), 0.0);

    for (size_t dir = 0; dir < itsResultSteps.size(); ++dir) {
      itsResultSteps[dir]->clear();
    }
    itsSolIntBuffers.clear();
  }

  ++itsTimeStep;
  itsTimer.stop();

  return false;
}

void DDECal::doPrepare(const DPBuffer& bufin, size_t sol_int, size_t step) {
  // UVW flagging happens on the copy of the buffer
  // These flags are later restored and therefore not written
  itsUVWFlagStep.process(bufin);

  itsTimerPredict.start();

  itsThreadPool->For(0, itsSteps.size(), [&](size_t dir, size_t) {
    itsSteps[dir]->process(bufin);
  });

  // Handle weights and flags
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = 4;

  size_t nchanblocks = itsChanBlockFreqs.size();

  double weightFactor =
      1. / (nCh * (info().antennaUsed().size() - 1) * nCr * itsRequestedSolInt);

  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t chanblock = 0, ant1 = info().antennaMap()[info().getAnt1()[bl]],
           ant2 = info().antennaMap()[info().getAnt2()[bl]];
    for (size_t ch = 0; ch < nCh; ++ch) {
      if (ch == itsChanBlockStart[chanblock + 1]) {
        chanblock++;
      }
      for (size_t cr = 0; cr < nCr; ++cr) {
        itsVisInInterval[chanblock].second++;  // total nr of vis
        if (bufin.getFlags()(cr, ch, bl)) {
          // Flagged points: set weight to 0
          DPBuffer& buf_sol = itsSolIntBuffers[sol_int].DataBuffers()[step];
          buf_sol.getWeights()(cr, ch, bl) = 0;
        } else {
          // Add this weight to both involved antennas
          double weight = bufin.getWeights()(cr, ch, bl);
          itsWeightsPerAntenna[ant1 * nchanblocks + chanblock] += weight;
          itsWeightsPerAntenna[ant2 * nchanblocks + chanblock] += weight;
          if (weight != 0.0) {
            itsVisInInterval[chanblock].first++;  // unflagged nr of vis
          }
        }
      }
    }
  }

  for (auto& weight : itsWeightsPerAntenna) {
    weight *= weightFactor;
  }

  itsTimerPredict.stop();

  itsAvgTime += itsAvgTime + bufin.getTime();
}

std::vector<std::pair<double, double>> DDECal::GetSourcePositions() const {
  std::vector<std::pair<double, double>> source_positions;
  source_positions.reserve(itsSteps.size());
  for (const std::shared_ptr<ModelDataStep>& s : itsSteps) {
    source_positions.push_back(s->GetFirstDirection());
  }
  return source_positions;
}

void DDECal::WriteSolutions() {
  itsTimer.start();
  itsTimerWrite.start();

  // Create antenna info for H5Parm, used antennas only.
  std::vector<std::string> used_antenna_names;
  used_antenna_names.reserve(info().antennaUsed().size());
  for (size_t used_antenna : info().antennaUsed()) {
    used_antenna_names.emplace_back(info().antennaNames()[used_antenna]);
  }

  const std::string history = "CREATE by " + DP3Version::AsString() + "\n" +
                              "step " + itsSettings.name + " in parset: \n" +
                              itsSettings.parset_string;

  itsSolutionWriter.Write(itsSols, itsConstraintSols, info().startTime(),
                          info().timeInterval() * itsSettings.solution_interval,
                          itsSettings.mode, used_antenna_names,
                          GetSourcePositions(), itsDirections,
                          info().chanFreqs(), itsChanBlockFreqs, history);

  itsTimerWrite.stop();
  itsTimer.stop();
}

void DDECal::finish() {
  itsTimer.start();

  if (itsSolIntBuffers.size() > 0) {
    doSolve();
  }

  if (!itsSettings.only_predict) WriteSolutions();

  itsSolIntBuffers.clear();
  itsTimer.stop();

  // Let the next steps finish.
  getNextStep()->finish();
}

void DDECal::storeModelData(
    const std::vector<std::vector<DPBuffer>>& input_model_buffers) {
  const size_t n_times = input_model_buffers.size();
  const size_t n_directions = itsSteps.size();
  // The model data might already have data in it. By using resize(),
  // reallocation between timesteps is avoided.
  itsModelData.resize(n_times);
  for (size_t timestep = 0; timestep != n_times; ++timestep) {
    itsModelData[timestep].resize(n_directions);
    for (size_t dir = 0; dir < n_directions; ++dir) {
      const casacore::Cube<std::complex<float>>& input_model_data =
          input_model_buffers[timestep][dir].getData();
      itsModelData[timestep][dir].assign(input_model_data.begin(),
                                         input_model_data.end());
    }
  }
}

void DDECal::subtractCorrectedModel(size_t bufferIndex) {
  // The original data is still in the data buffers (the solver
  // doesn't change those). Here we apply the solutions to all the model data
  // directions and subtract them from the data.
  std::vector<std::vector<casacore::DComplex>>& solutions =
      itsSols[itsSolIntBuffers[bufferIndex].NSolution()];
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nDir = itsDirections.size();
  for (size_t time = 0; time != itsSolIntBuffers[bufferIndex].Size(); ++time) {
    DPBuffer& data_buffer = itsSolIntBuffers[bufferIndex].DataBuffers()[time];
    const std::vector<std::vector<std::complex<float>>>& modelData =
        itsModelData[time];
    for (size_t bl = 0; bl < nBl; ++bl) {
      const size_t ant1 = info().getAnt1()[bl];
      const size_t ant2 = info().getAnt2()[bl];
      size_t chanblock = 0;

      for (size_t ch = 0; ch < nCh; ++ch) {
        if (ch == itsChanBlockStart[chanblock + 1]) {
          chanblock++;
        }

        std::complex<float>* data = &data_buffer.getData()(0, ch, bl);
        const size_t index = data - data_buffer.getData().data();
        if (itsSettings.only_predict) {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());

          for (size_t dir = 0; dir != nDir; ++dir)
            value += aocommon::MC2x2(&modelData[dir][index]);

          for (size_t cr = 0; cr < 4; ++cr) data[cr] = value[cr];
        } else {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());
          for (size_t dir = 0; dir != nDir; ++dir) {
            if (itsSolver->NSolutionPolarizations() == 4) {
              const aocommon::MC2x2 sol1(
                  &solutions[chanblock][(ant1 * nDir + dir) * 4]);
              const aocommon::MC2x2 sol2(
                  &solutions[chanblock][(ant2 * nDir + dir) * 4]);
              value += sol1.Multiply(aocommon::MC2x2(&modelData[dir][index]))
                           .MultiplyHerm(sol2);
            } else if (itsSolver->NSolutionPolarizations() == 2) {
              const aocommon::MC2x2 sol1(
                  solutions[chanblock][(ant1 * nDir + dir) * 2 + 0], 0.0, 0.0,
                  solutions[chanblock][(ant1 * nDir + dir) * 2 + 1]);
              const aocommon::MC2x2 sol2(
                  solutions[chanblock][(ant2 * nDir + dir) * 2 + 0], 0.0, 0.0,
                  solutions[chanblock][(ant2 * nDir + dir) * 2 + 1]);
              value += sol1.Multiply(aocommon::MC2x2(&modelData[dir][index]))
                           .MultiplyHerm(sol2);
            } else {
              std::complex<double> solfactor(
                  solutions[chanblock][ant1 * nDir + dir] *
                  std::conj(solutions[chanblock][ant2 * nDir + dir]));
              for (size_t cr = 0; cr < 4; ++cr)
                value[cr] += solfactor *
                             std::complex<double>(modelData[dir][index + cr]);
            }
          }
          for (size_t cr = 0; cr < 4; ++cr) data[cr] -= value[cr];
        }
      }  // channel loop
    }    // bl loop
  }      // time loop
}

}  // namespace steps
}  // namespace dp3
