//# DDECal.cc: DPPP step class to do a direction dependent gain calibration
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
//# $Id: DDECal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/Package__Version.h>
#include <DPPP_DDECal/DDECal.h>

#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/DPLogger.h>
#include <DPPP/MSReader.h>
#include <DPPP/Simulate.h>
#include <DPPP/SourceDBUtil.h>

#include <DPPP_DDECal/ScreenConstraint.h>
#include <DPPP_DDECal/TECConstraint.h>
#include <DPPP_DDECal/RotationConstraint.h>
#include <DPPP_DDECal/RotationAndDiagonalConstraint.h>

#include <ParmDB/ParmDB.h>
#include <ParmDB/ParmValue.h>
#include <ParmDB/SourceDB.h>

#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>
#include <Common/ParameterSet.h>
#include <Common/StreamUtil.h>
#include <Common/StringUtil.h>

#include <fstream>
#include <ctime>
#include <utility>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/OS/File.h>

#include <vector>
#include <algorithm>

#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace casacore;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    DDECal::DDECal (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsAvgTime       (0),
        itsSols          (),
        itsH5ParmName    (parset.getString (prefix + "h5parm",
                                            parset.getString("msin")+
                                              "/instrument.h5")),
        itsH5Parm        (itsH5ParmName, true),
        itsPropagateSolutions (parset.getBool (prefix + "propagatesolutions",
                                               false)),
        itsTimeStep      (0),
        itsSolInt        (parset.getInt (prefix + "solint", 1)),
        itsStepInSolInt  (0),
        itsNChan         (parset.getInt (prefix + "nchan", 1)),
        itsUVWFlagStep   (input, parset, prefix),
        itsCoreConstraint(parset.getDouble (prefix + "coreconstraint", 0.0)),
        itsScreenCoreConstraint(parset.getDouble (prefix + "tecscreen.coreconstraint", 0.0)),
        itsFullMatrixMinimalization(false),
        itsApproximateTEC(false),
        itsStatFilename(parset.getString(prefix + "statfilename", ""))
    {
      stringstream ss;
      ss << parset;
      itsParsetString = ss.str();

      vector<string> strDirections = 
         parset.getStringVector (prefix + "directions",
                                 vector<string> ());
         
      itsMultiDirSolver.set_max_iterations(parset.getInt(prefix + "maxiter", 50));
      double tolerance = parset.getDouble(prefix + "tolerance", 1.e-5);
      itsMultiDirSolver.set_accuracy(tolerance);
      itsMultiDirSolver.set_constraint_accuracy(parset.getDouble(prefix + "approxtolerance", tolerance*10.0));
      itsMultiDirSolver.set_step_size(parset.getDouble(prefix + "stepsize", 0.2));
      itsMultiDirSolver.set_detect_stalling(parset.getBool(prefix + "detectstalling", true));

      if(!itsStatFilename.empty())
	itsStatStream.reset(new std::ofstream(itsStatFilename));
      
      // Default directions are all patches
      if (strDirections.empty()) {
        string sourceDBName = parset.getString(prefix+"sourcedb");
        BBS::SourceDB sourceDB(BBS::ParmDBMeta("", sourceDBName), false);
        vector<string> patchNames = makePatchList(sourceDB, vector<string>());
        itsDirections.resize(patchNames.size());
        for (uint i=0; i<patchNames.size(); ++i) {
          itsDirections[i] = vector<string>(1, patchNames[i]);
        }
      } else {
        itsDirections.resize(strDirections.size());
        for (uint i=0; i<strDirections.size(); ++i) {
          ParameterValue dirStr(strDirections[i]);
          itsDirections[i] = dirStr.getStringVector();
        }
      }

      itsMode = GainCal::stringToCalType(
                   toLower(parset.getString(prefix + "mode",
                                            "complexgain")));
      if(itsCoreConstraint != 0.0) {
        itsConstraints.push_back(casacore::CountedPtr<Constraint>(
          new CoreConstraint()));
      }
      switch(itsMode) {
        case GainCal::COMPLEXGAIN:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new DiagonalConstraint(4)));
          itsMultiDirSolver.set_phase_only(false);
          itsFullMatrixMinimalization = true;
          break;
        case GainCal::SCALARCOMPLEXGAIN:
          // no constraints
          itsMultiDirSolver.set_phase_only(false);
          itsFullMatrixMinimalization = false;
          break;
        case GainCal::FULLJONES:
          // no constraints
          itsMultiDirSolver.set_phase_only(false);
          itsFullMatrixMinimalization = true;
          break;
        case GainCal::PHASEONLY:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new PhaseOnlyConstraint()));
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new DiagonalConstraint(4)));
          itsMultiDirSolver.set_phase_only(true);
          itsFullMatrixMinimalization = true;
          break;
        case GainCal::SCALARPHASE:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new PhaseOnlyConstraint()));
          itsMultiDirSolver.set_phase_only(true);
          break;
        case GainCal::AMPLITUDEONLY:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new DiagonalConstraint(4)));
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new AmplitudeOnlyConstraint()));
          itsMultiDirSolver.set_phase_only(false);
          itsFullMatrixMinimalization = true;
          break;
        case GainCal::SCALARAMPLITUDE:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new AmplitudeOnlyConstraint()));
          itsMultiDirSolver.set_phase_only(false);
          itsFullMatrixMinimalization = false;
          break;
        case GainCal::TEC:
        case GainCal::TECANDPHASE:
          itsApproximateTEC = parset.getBool(prefix + "approximatetec", false);
          if(itsApproximateTEC)
          {
            int iters = parset.getInt(prefix + "maxapproxiter", itsMultiDirSolver.max_iterations()/2);
            int chunksize = parset.getInt(prefix + "approxchunksize", 0);
            casacore::CountedPtr<ApproximateTECConstraint> ptr;
            if(itsMode == GainCal::TEC)
              ptr = casacore::CountedPtr<ApproximateTECConstraint>(
                new ApproximateTECConstraint(TECConstraint::TECOnlyMode));
            else
              ptr = casacore::CountedPtr<ApproximateTECConstraint>(
                new ApproximateTECConstraint(TECConstraint::TECAndCommonScalarMode));
            ptr->SetMaxApproximatingIterations(iters);
            ptr->SetFittingChunkSize(chunksize);
            itsConstraints.push_back(ptr);
          }
          else {
            if(itsMode == GainCal::TEC)
              itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                new TECConstraint(TECConstraint::TECOnlyMode)));
              else
              itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                new TECConstraint(TECConstraint::TECAndCommonScalarMode)));
          }
          itsMultiDirSolver.set_phase_only(true);
          itsFullMatrixMinimalization = false;
          break;
        case GainCal::TECSCREEN:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                    new ScreenConstraint(parset, prefix+"tecscreen.")));
          itsMultiDirSolver.set_phase_only(true);
          itsFullMatrixMinimalization = false;
          break;
        case GainCal::ROTATIONANDDIAGONAL:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
              new RotationAndDiagonalConstraint()));
          itsFullMatrixMinimalization = true;
          break;
        case GainCal::ROTATION:
          itsConstraints.push_back(casacore::CountedPtr<Constraint>(
              new RotationConstraint()));
          itsFullMatrixMinimalization = true;
          break;
        default:
          THROW (Exception, "Unexpected mode: " << 
                          GainCal::calTypeToString(itsMode));
      }

      const size_t nDir = itsDirections.size();
      itsPredictSteps.resize(nDir);
      for (size_t dir=0; dir<nDir; ++dir) {
        itsPredictSteps[dir]=Predict(input, parset, prefix, itsDirections[dir]);
      }
    }

    DDECal::~DDECal()
    {}

    DPStep::ShPtr DDECal::makeStep (DPInput* input,
                                    const ParameterSet& parset,
                                    const std::string& prefix)
    {
      return DPStep::ShPtr(new DDECal(input, parset, prefix));
    }

    void DDECal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();

      const size_t nDir=itsDirections.size();

      itsUVWFlagStep.updateInfo(infoIn);
      for (size_t dir=0; dir<nDir; ++dir) {
        itsPredictSteps[dir].updateInfo(infoIn);
      }

      if (itsSolInt==0) {
        itsSolInt=info().ntime();
      }

      itsDataPtrs.resize(itsSolInt);
      itsModelDataPtrs.resize(itsSolInt);
      for (uint t=0; t<itsSolInt; ++t) {
        itsModelDataPtrs[t].resize(nDir);
      }
      for  (uint i=0;i<itsConstraints.size();i++) {
        itsMultiDirSolver.add_constraint(itsConstraints[i].get());
      }


      itsBufs.resize(itsSolInt);

      itsDataResultStep = ResultStep::ShPtr(new ResultStep());
      itsUVWFlagStep.setNextStep(itsDataResultStep);

      itsResultSteps.resize(nDir);
      for (size_t dir=0; dir<nDir; ++dir) {
        itsResultSteps[dir] = MultiResultStep::ShPtr(new MultiResultStep(itsSolInt));
        itsPredictSteps[dir].setNextStep(itsResultSteps[dir]);
      }

      if (itsNChan==0 || itsNChan>info().nchan()) {
        itsNChan = info().nchan();
      }

      // Convert from casacore::Vector to std::vector
      vector<int> ant1(info().getAnt1().size());
      vector<int> ant2(info().getAnt2().size());
      for (uint i=0; i<ant1.size(); ++i) {
        ant1[i]=info().getAnt1()[i];
        ant2[i]=info().getAnt2()[i];
      }

      // Fill antenna info in H5Parm, need to convert from casa types to std types
      std::vector<std::string> antennaNames(info().antennaNames().size());
      std::vector<std::vector<double> > antennaPos(info().antennaPos().size());
      for (uint i=0; i<info().antennaNames().size(); ++i) {
        antennaNames[i]=info().antennaNames()[i];
        casacore::Quantum<casacore::Vector<double> > pos = info().antennaPos()[i].get("m");
        antennaPos[i].resize(3);
        antennaPos[i][0] = pos.getValue()[0];
        antennaPos[i][1] = pos.getValue()[1];
        antennaPos[i][2] = pos.getValue()[2];
      }

      itsH5Parm.addAntennas(antennaNames, antennaPos);

      std::vector<std::pair<double, double> > sourcePositions(itsDirections.size());
      for (uint i=0; i<itsDirections.size(); ++i) {
        sourcePositions[i] = itsPredictSteps[i].getFirstDirection();
      }
      itsH5Parm.addSources(getDirectionNames(), sourcePositions);

      uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
      itsSols.resize(nSolTimes);
      itsNIter.resize(nSolTimes);
      itsNApproxIter.resize(nSolTimes);
      itsConstraintSols.resize(nSolTimes);

      size_t nChannelBlocks = info().nchan()/itsNChan;
      itsChanBlockFreqs.resize(nChannelBlocks);
      for(size_t chBlock=0; chBlock!=nChannelBlocks; ++chBlock) {
        const size_t
          channelIndexStart = chBlock * info().nchan() / nChannelBlocks,
          channelIndexEnd = (chBlock+1) * info().nchan() / nChannelBlocks,
          curChannelBlockSize = channelIndexEnd - channelIndexStart;
        double  meanfreq = std::accumulate(
            info().chanFreqs().data()+channelIndexStart,
            info().chanFreqs().data()+channelIndexEnd,
            0.0) / curChannelBlockSize;
        itsChanBlockFreqs[chBlock] = meanfreq;
      }

      for (uint i=0; i<itsConstraints.size();++i) {
        // Different constraints need different information. Determine if the constraint is
        // of a type that needs more information, and if so initialize the constraint.
        CoreConstraint* coreConstraint = dynamic_cast<CoreConstraint*>(itsConstraints[i].get());
        if(coreConstraint != 0)
        {
          double
            refX = antennaPos[i][0],
            refY = antennaPos[i][1],
            refZ = antennaPos[i][2];
          std::set<size_t> coreAntennaIndices;
          const double coreDistSq = itsCoreConstraint*itsCoreConstraint;
          for(size_t ant=0; ant!=antennaPos.size(); ++ant)
          {
            double
              dx = refX - antennaPos[ant][0],
              dy = refY - antennaPos[ant][1],
              dz = refZ - antennaPos[ant][2],
              distSq = dx*dx + dy*dy + dz*dz;
            if(distSq <= coreDistSq)
              coreAntennaIndices.insert(ant);
          }
          coreConstraint->initialize(info().antennaNames().size(),
            itsDirections.size(),
            info().nchan(),
            coreAntennaIndices);
        }
        
        ScreenConstraint* screenConstraint = dynamic_cast<ScreenConstraint*>(itsConstraints[i].get());
        if(screenConstraint != 0)
        {
          screenConstraint->initialize(
              info().antennaNames().size(),
              itsDirections.size(),
              nChannelBlocks,
              &(itsChanBlockFreqs[0])
          );
          screenConstraint->setAntennaPositions(antennaPos);
          screenConstraint->setDirections(sourcePositions);
          screenConstraint->initPiercePoints();
          double
            refX = antennaPos[i][0],
            refY = antennaPos[i][1],
            refZ = antennaPos[i][2];
          std::vector<size_t> coreAntennaIndices;
          std::vector<size_t> otherAntennaIndices;
          const double coreDistSq = itsScreenCoreConstraint*itsScreenCoreConstraint;
          for(size_t ant=0; ant!=antennaPos.size(); ++ant)
          {
            double
              dx = refX - antennaPos[ant][0],
              dy = refY - antennaPos[ant][1],
              dz = refZ - antennaPos[ant][2],
              distSq = dx*dx + dy*dy + dz*dz;
            if(distSq <= coreDistSq)
              coreAntennaIndices.push_back(ant);
            else
              otherAntennaIndices.push_back(ant);
          }
          screenConstraint->setCoreAntennas(coreAntennaIndices);
          screenConstraint->setOtherAntennas(otherAntennaIndices);
        }
        
        TECConstraintBase* tecConstraint = dynamic_cast<TECConstraintBase*>(itsConstraints[i].get());
        if(tecConstraint != 0)
        {
          tecConstraint->initialize(info().antennaNames().size(),
              itsDirections.size(),
              nChannelBlocks,
              &(itsChanBlockFreqs[0]));
        }

        RotationAndDiagonalConstraint* rotationAndDiagonalConstraint = dynamic_cast<RotationAndDiagonalConstraint*>(itsConstraints[i].get());
        if(rotationAndDiagonalConstraint != 0)
        {
          rotationAndDiagonalConstraint->initialize(info().antennaNames().size(), itsDirections.size(), nChannelBlocks);
        }
        RotationConstraint* rotationConstraint = dynamic_cast<RotationConstraint*>(itsConstraints[i].get());
        if(rotationConstraint != 0)
        {
          rotationConstraint->initialize(info().antennaNames().size(), itsDirections.size(), nChannelBlocks);
        }
      }

      uint nSt = info().antennaNames().size();
      itsMultiDirSolver.init(nSt, nDir, info().nchan(), ant1, ant2);
      itsMultiDirSolver.set_channel_blocks(nChannelBlocks);

      for (uint i=0; i<nSolTimes; ++i) {
        itsSols[i].resize(nChannelBlocks);
      }
    }

    void DDECal::show (std::ostream& os) const
    {
      os << "DDECal " << itsName << endl;
      os << "  H5Parm:              " << itsH5ParmName <<endl;
      os << "  solint:              " << itsSolInt <<endl;
      os << "  nchan:               " << itsNChan <<endl;
      os << "  directions:          " << itsDirections << endl;
      os << "  mode (constraints):  " << GainCal::calTypeToString(itsMode) 
         << endl;
      os << "  coreconstraint:      " << itsCoreConstraint << endl;
      os << "  approximate fitter:  " << itsApproximateTEC << endl;
      for (uint i=0; i<itsPredictSteps.size(); ++i) {
        itsPredictSteps[i].show(os);
      }
      itsUVWFlagStep.show(os);
    }

    void DDECal::showTimings (std::ostream& os, double duration) const
    {
      double totaltime=itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " DDECal " << itsName << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerPredict.getElapsed(), totaltime);
      os << " of it spent in predict" << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerSolve.getElapsed(), totaltime);
      os << " of it spent in estimating gains and computing residuals" << endl;

      itsMultiDirSolver.showTimings(os, itsTimerSolve.getElapsed());

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerWrite.getElapsed(), totaltime);
      os << " of it spent in writing gain solutions to disk" << endl;

      os << "Iterations taken: [";
      for (uint i=0; i<itsNIter.size()-1; ++i) {
        os<<itsNIter[i];
        if(itsNApproxIter[i]!=0)
          os << '|' << itsNApproxIter[i];
        os<<",";
      }
      os<<itsNIter[itsNIter.size()-1];
        if(itsNApproxIter[itsNIter.size()-1]!=0)
          os << '|' << itsNApproxIter[itsNIter.size()-1];
      os<<"]"<<endl;
    }

    void DDECal::initializeScalarSolutions() {
      if (itsTimeStep/itsSolInt>0 && itsPropagateSolutions) {
        // initialize solutions with those of the previous step
        itsSols[itsTimeStep/itsSolInt] = itsSols[itsTimeStep/itsSolInt-1];
      } else {
        // initialize solutions with 1.
        size_t n = itsDirections.size()*info().antennaNames().size();
        for (vector<DComplex>& solvec : itsSols[itsTimeStep/itsSolInt]) {
          solvec.assign(n, 1.0);
        }
      }
    }

    void DDECal::initializeFullMatrixSolutions() {
      if (itsTimeStep/itsSolInt>0 && itsPropagateSolutions) {
        // initialize solutions with those of the previous step
        itsSols[itsTimeStep/itsSolInt] = itsSols[itsTimeStep/itsSolInt-1];
      } else {
        // initialize solutions with unity matrix [1 0 ; 0 1].
        size_t n = itsDirections.size()*info().antennaNames().size();
        for (vector<DComplex>& solvec : itsSols[itsTimeStep/itsSolInt]) {
          solvec.resize(n*4);
          for(size_t i=0; i!=n; ++i)
          {
            solvec[i*4 + 0] = 1.0;
            solvec[i*4 + 1] = 0.0;
            solvec[i*4 + 2] = 0.0;
            solvec[i*4 + 3] = 1.0;
          }
        }
      }
    }
    
    vector<string> DDECal::getDirectionNames() {
      vector<string> res;
      
      for (vector<vector<string> >::iterator dirIter = itsDirections.begin();
           dirIter != itsDirections.end();
           dirIter++) {
        stringstream ss;
        ss << (*dirIter);
        res.push_back(ss.str());
      }
      return res;
    }

    void DDECal::doSolve ()
    {
      if(itsFullMatrixMinimalization)
        initializeFullMatrixSolutions();
      else
        initializeScalarSolutions();

      itsTimerSolve.start();
      MultiDirSolver::SolveResult solveResult;
      if(itsFullMatrixMinimalization)
      {
        solveResult = itsMultiDirSolver.processFullMatrix(itsDataPtrs, itsModelDataPtrs,
          itsSols[itsTimeStep/itsSolInt],
	  itsAvgTime / itsSolInt, itsStatStream.get());
      }
      else {
        solveResult = itsMultiDirSolver.processScalar(itsDataPtrs, itsModelDataPtrs,
          itsSols[itsTimeStep/itsSolInt],
	  itsAvgTime / itsSolInt, itsStatStream.get());
      }
      itsTimerSolve.stop();

      itsNIter[itsTimeStep/itsSolInt] = solveResult.iterations;
      itsNApproxIter[itsTimeStep/itsSolInt] = solveResult.constraintIterations;

      // Store constraint solutions if any constaint has a non-empty result
      bool someConstraintHasResult = false;
      for (uint constraintnum=0; constraintnum<solveResult._results.size(); ++constraintnum) {
        if (!solveResult._results[constraintnum].empty()) {
          someConstraintHasResult = true;
          break;
        }
      }
      if (someConstraintHasResult) {
        itsConstraintSols[itsTimeStep/itsSolInt]=solveResult._results;
      }
    }

    bool DDECal::process (const DPBuffer& bufin)
    {
      itsTimer.start();

      itsBufs[itsStepInSolInt].copy(bufin);
      itsDataPtrs[itsStepInSolInt] = itsBufs[itsStepInSolInt].getData().data();

      // Fetch inputs because parallel PredictSteps should not read it from disk
      itsInput->fetchUVW(bufin, itsBufs[itsStepInSolInt], itsTimer);
      itsInput->fetchWeights(bufin, itsBufs[itsStepInSolInt], itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBufs[itsStepInSolInt], itsTimer);

      // UVW flagging happens on a copy of the buffer, so these flags are not written
      itsUVWFlagStep.process(itsBufs[itsStepInSolInt]);

      itsTimerPredict.start();

//      if(itsPredictSteps.size() < LOFAR::OpenMP::maxThreads())
//        LOFAR::OpenMP::setNested(true);
#pragma omp parallel for schedule(dynamic) if(itsPredictSteps.size()>1)
      for (size_t dir=0; dir<itsPredictSteps.size(); ++dir) {
        itsPredictSteps[dir].process(itsBufs[itsStepInSolInt]);
        itsModelDataPtrs[itsStepInSolInt][dir] =
                 itsResultSteps[dir]->get()[itsStepInSolInt].getData().data();
      }

      // Handle weights and flags
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;
      
      for (size_t ch=0; ch<nCh; ++ch) {
        for (size_t bl=0; bl<nBl; ++bl) {
          for (size_t cr=0; cr<nCr; ++cr) {
            if (itsBufs[itsStepInSolInt].getFlags().data()[bl*nCr*nCh+ch*nCr+cr]) {
              // Flagged points: set data and model to 0
              itsDataPtrs[itsStepInSolInt][bl*nCr*nCh+ch*nCr+cr] = 0;
              for (size_t dir=0; dir<itsPredictSteps.size(); ++dir) {
                itsModelDataPtrs[itsStepInSolInt][dir][bl*nCr*nCh+ch*nCr+cr] = 0;
              }
            } else {
              // Premultiply non-flagged data with sqrt(weight)
              double weight = sqrt(itsBufs[itsStepInSolInt].getWeights().data()[bl*nCr*nCh+ch*nCr+cr]);
              itsDataPtrs[itsStepInSolInt][bl*nCr*nCh+ch*nCr+cr] *= weight;
              for (size_t dir=0; dir<itsPredictSteps.size(); ++dir) {
                itsModelDataPtrs[itsStepInSolInt][dir][bl*nCr*nCh+ch*nCr+cr] *= weight;
              }
            }
          }
        }
      }

      itsTimerPredict.stop();

      itsAvgTime += itsAvgTime + bufin.getTime();

      if (itsStepInSolInt==itsSolInt-1) {
        doSolve();
            itsStepInSolInt=0;
        itsAvgTime=0;
        for (size_t dir=0; dir<itsResultSteps.size(); ++dir) {
          itsResultSteps[dir]->clear();
        }
      } else {
        itsStepInSolInt++;
      }

      itsTimeStep++;
      itsTimer.stop();

      getNextStep()->process(bufin);

      return false;
    }

    void DDECal::writeSolutions() {
      itsTimer.start();
      itsTimerWrite.start();

      uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
      uint nDir = itsDirections.size();
      ASSERT(nSolTimes==itsSols.size());
      vector<double> solTimes(nSolTimes);
      double starttime=info().startTime();
      for (uint t=0; t<nSolTimes; ++t) {
        solTimes[t] = starttime+(t+0.5)*info().timeInterval()*itsSolInt;
      }

      if (itsConstraintSols[0].empty()) {
        // Record the actual iterands of the solver, not constraint results

        uint nPol;

        vector<string> polarizations;
        if(itsMode == GainCal::COMPLEXGAIN ||
           itsMode == GainCal::PHASEONLY ||
           itsMode == GainCal::AMPLITUDEONLY) {
          nPol = 2;
          polarizations.push_back("XX");
          polarizations.push_back("YY");
        } else if (itsMode == GainCal::FULLJONES) {
          polarizations.push_back("XX");
          polarizations.push_back("XY");
          polarizations.push_back("YX");
          polarizations.push_back("YY");
          nPol = 4;
        } else {
          nPol = 1;
        }

        uint nSolChan = itsSols[0].size();
        ASSERT(nSolChan == itsChanBlockFreqs.size());

        vector<DComplex> sols(nSolChan*info().nantenna()*nSolTimes*nDir*nPol);
        size_t i=0;

        // For nPol=1, loop over pol runs just once
        // For nPol=2, it runs over values 0 and 2 (picking diagonal elements from 4 pols)
        // For nPol=4, it runs over 0, 1, 2, 3
        uint polIncr= (itsMode==GainCal::FULLJONES?1:3);
        uint maxPol = (nPol>1?4:1);
        // Put solutions in a contiguous piece of memory
        for (uint time=0; time<nSolTimes; ++time) {
          for (uint chan=0; chan<nSolChan; ++chan) {
            for (uint ant=0; ant<info().nantenna(); ++ant) {
              for (uint dir=0; dir<nDir; ++dir) {
                for (uint pol=0; pol<maxPol; pol+=polIncr) {
                  ASSERT(!itsSols[time].empty());
                  ASSERT(!itsSols[time][chan].empty());
                  ASSERT(time<itsSols.size());
                  ASSERT(chan<itsSols[time].size());
                  ASSERT(ant*nDir*maxPol+dir*maxPol+pol<itsSols[time][chan].size());
                  ASSERT(i<sols.size());
                  sols[i] = itsSols[time][chan][ant*nDir*maxPol+dir*maxPol+pol];
                  ++i;
                }
              }
            }
          }
        }
        vector<H5Parm::AxisInfo> axes;
        axes.push_back(H5Parm::AxisInfo("time", itsSols.size()));
        axes.push_back(H5Parm::AxisInfo("freq", nSolChan));
        axes.push_back(H5Parm::AxisInfo("ant", info().nantenna()));
        axes.push_back(H5Parm::AxisInfo("dir", nDir));
        if (nPol>1) {
          axes.push_back(H5Parm::AxisInfo("pol", nPol));
        }

        string historyString = "CREATE by DPPP\n" +
            Version::getInfo<DPPPVersion>("DPPP", "top") + "\n" +
            "step " + itsName + " in parset: \n" + itsParsetString;
        uint numsols = 1;
        // For [scalar]complexgain, store two soltabs: phase and amplitude
        if (itsMode == GainCal::COMPLEXGAIN ||
            itsMode == GainCal::SCALARCOMPLEXGAIN || itsMode == GainCal::FULLJONES) {
          numsols = 2;
        }
        for (uint solnum=0; solnum<numsols; ++solnum) {
          string solTabName;
          H5Parm::SolTab soltab;
          switch (itsMode) {
            case GainCal::SCALARPHASE:
            case GainCal::PHASEONLY:
            case GainCal::FULLJONES:
              if (solnum==0) {
                solTabName = "phase000";
                soltab = itsH5Parm.createSolTab(solTabName, "phase", axes);
                soltab.setComplexValues(sols, vector<double>(), false, historyString);
              } else {
                solTabName = "amplitude000";
                soltab = itsH5Parm.createSolTab(solTabName, "amplitude", axes);
                soltab.setComplexValues(sols, vector<double>(), true, historyString);
              }
              break;
            case GainCal::SCALARCOMPLEXGAIN:
            case GainCal::COMPLEXGAIN:
              if (solnum==0) {
                solTabName = "phase000";
                soltab = itsH5Parm.createSolTab(solTabName, "phase", axes);
                soltab.setComplexValues(sols, vector<double>(), false, historyString);
              } else {
                solTabName = "amplitude000";
                soltab = itsH5Parm.createSolTab(solTabName, "amplitude", axes);
                soltab.setComplexValues(sols, vector<double>(), true, historyString);
              }
              break;
            case GainCal::SCALARAMPLITUDE:
            case GainCal::AMPLITUDEONLY:
              solTabName = "amplitude000";
              soltab = itsH5Parm.createSolTab(solTabName, "amplitude", axes);
              soltab.setComplexValues(sols, vector<double>(), true, historyString);
              break;
            default: 
              THROW(Exception, "Constraint should have produced output");
          }

          // Tell H5Parm that all antennas and directions were used 
          std::vector<std::string> antennaNames(info().antennaNames().size());
          for (uint i=0; i<info().antennaNames().size(); ++i) {
            antennaNames[i]=info().antennaNames()[i];
          }
          soltab.setAntennas(antennaNames);
    
          soltab.setSources(getDirectionNames());

          if (nPol>1) {
            soltab.setPolarizations(polarizations);
          }
   
          soltab.setFreqs(itsChanBlockFreqs);

          soltab.setTimes(solTimes);
        } // solnums loop
      } else {
        // Record the Constraint::Result in the H5Parm

        uint nConstraints = itsConstraintSols[0].size();

        for (uint constraintNum=0; constraintNum<nConstraints; ++constraintNum) {
          // Number of solution names, e.g. 2 for "TEC" and "ScalarPhase"
          uint nSolNames = itsConstraintSols[0][constraintNum].size();
          for (uint solNameNum=0; solNameNum<nSolNames; ++solNameNum) {
            // Get the result of the constraint solution at first time to get metadata
            Constraint::Result firstResult = itsConstraintSols[0][constraintNum][solNameNum];

            vector<hsize_t> dims(firstResult.dims.size()+1); // Add time dimension at beginning
            dims[0]=itsConstraintSols.size(); // Number of times
            size_t numSols=dims[0];
            for (uint i=1; i<dims.size(); ++i) {
              dims[i] = firstResult.dims[i-1];
              numSols *= dims[i];
            }

            vector<string> firstaxesnames = StringUtil::tokenize(firstResult.axes,",");

            vector<H5Parm::AxisInfo> axes;
            axes.push_back(H5Parm::AxisInfo("time", itsConstraintSols.size()));
            for (size_t axisNum=0; axisNum<firstaxesnames.size(); ++axisNum) {
              axes.push_back(H5Parm::AxisInfo(firstaxesnames[axisNum], firstResult.dims[axisNum]));
            }

            // Put solutions in a contiguous piece of memory
            vector<double> sols(numSols);
            vector<double>::iterator nextpos = sols.begin();
            for (uint time=0; time<itsSols.size(); ++time) {
              nextpos = std::copy(
                itsConstraintSols[time][constraintNum][solNameNum].vals.begin(),
                itsConstraintSols[time][constraintNum][solNameNum].vals.end(),
                nextpos);
            }

            string solTabName = firstResult.name+"000";
            H5Parm::SolTab soltab = itsH5Parm.createSolTab(solTabName, firstResult.name, axes);
            soltab.setValues(sols, vector<double>(),
                             "CREATE by DPPP\n" +
                             Version::getInfo<DPPPVersion>("DPPP", "top") + "\n" +
                             "step " + itsName + " in parset: \n" +
                             itsParsetString);

            // Tell H5Parm that all antennas and directions were used 
            std::vector<std::string> antennaNames(info().antennaNames().size());
            for (uint i=0; i<info().antennaNames().size(); ++i) {
              antennaNames[i]=info().antennaNames()[i];
            }
            soltab.setAntennas(antennaNames);
      
            soltab.setSources(getDirectionNames());
 
            if (soltab.hasAxis("pol")) {
              vector<string> polarizations;
              switch (soltab.getAxis("pol").size) {
                case 2: polarizations.push_back("XX");
                        polarizations.push_back("YY");
                        break;
                case 4: polarizations.push_back("XX");
                        polarizations.push_back("XY");
                        polarizations.push_back("YX");
                        polarizations.push_back("YY");
                        break;
                default:
                        THROW (Exception, "No metadata for numpolarizations = " << soltab.getAxis("pol").size);
              }
              
              soltab.setPolarizations(polarizations);
            }

            // Set channel to frequencies
            // Do not use itsChanBlockFreqs, because constraint may have changed size
            uint nChannelBlocks = 1;
            if (soltab.hasAxis("freq")) {
              nChannelBlocks = soltab.getAxis("freq").size;
            }
            vector<double> chanBlockFreqs;

            chanBlockFreqs.resize(nChannelBlocks);
            for(size_t chBlock=0; chBlock!=nChannelBlocks; ++chBlock) {
              const size_t
                channelIndexStart = chBlock * info().nchan() / nChannelBlocks,
                channelIndexEnd = (chBlock+1) * info().nchan() / nChannelBlocks,
                curChannelBlockSize = channelIndexEnd - channelIndexStart;
              double  meanfreq = std::accumulate(
                  info().chanFreqs().data()+channelIndexStart,
                  info().chanFreqs().data()+channelIndexEnd,
                  0.0) / curChannelBlockSize;
              chanBlockFreqs[chBlock] = meanfreq;
            }

            soltab.setFreqs(chanBlockFreqs);
 
            soltab.setTimes(solTimes);
            // End TODO 
          }
        }
      }

      itsTimerWrite.stop();
      itsTimer.stop();
    }

    void DDECal::finish()
    {
      itsTimer.start();

      if (itsStepInSolInt!=0) {
        //shrink itsDataPtrs, itsModelDataPtrs
        std::vector<casacore::Complex*>(itsDataPtrs.begin(),
            itsDataPtrs.begin()+itsStepInSolInt).swap(itsDataPtrs);
        std::vector<std::vector<casacore::Complex*> >(itsModelDataPtrs.begin(),
                    itsModelDataPtrs.begin()+itsStepInSolInt).swap(itsModelDataPtrs);

        doSolve();
      }

      writeSolutions();

      itsTimer.stop();

      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
