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
#include <DPPP_DDECal/DDECal.h>
#include <DPPP/Simulate.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/MSReader.h>
#include <DPPP/DPLogger.h>
#include <DPPP_DDECal/ScreenConstraint.h>
#include <ParmDB/ParmDB.h>
#include <ParmDB/ParmValue.h>
#include <ParmDB/SourceDB.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>

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
        itsH5Parm        (itsH5ParmName),
        itsPropagateSolutions (parset.getBool (prefix + "propagatesolutions",
                                               false)),
        itsTimeStep      (0),
        itsSolInt        (parset.getInt (prefix + "solint", 1)),
        itsStepInSolInt  (0),
        itsNChan         (parset.getInt (prefix + "nchan", 0)),
        itsNFreqCells    (0),
        itsCoreConstraint(parset.getDouble (prefix + "coreconstraint", 0.0)),
        itsScreenCoreConstraint(parset.getDouble (prefix + "tecscreen.coreconstraint", 0.0)),
        itsMultiDirSolver(parset.getInt (prefix + "maxiter", 50),
                          parset.getDouble (prefix + "tolerance", 1.e-5),
                          parset.getDouble (prefix + "stepsize", 0.2))
    {
      vector<string> strDirections = 
         parset.getStringVector (prefix + "directions",
                                 vector<string> ());

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
      if (itsMode == GainCal::PHASEONLY) {
        itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                  new PhaseConstraint()));
      } else if (itsMode == GainCal::TEC) {
        itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                  new TECConstraint(TECConstraint::TECOnlyMode)));
      } else if (itsMode == GainCal::TECANDPHASE){
        itsConstraints.push_back(casacore::CountedPtr<Constraint>(
                  new TECConstraint(TECConstraint::TECAndCommonScalarMode)));
      } else if (itsMode == GainCal::COMPLEXGAIN) {
        // no constraints
      } else if (itsMode == GainCal::TECSCREEN) {
        itsConstraints.push_back(casacore::CountedPtr<Constraint>(
		  new ScreenConstraint(parset, prefix+"tecscreen.")));
      } else {
        THROW (Exception, "Unexpected mode: " << 
                          GainCal::calTypeToString(itsMode));
      }

      itsDataPtrs.resize(itsSolInt);
      itsModelDataPtrs.resize(itsSolInt);
      for (uint t=0; t<itsSolInt; ++t) {
        itsModelDataPtrs[t].resize(itsDirections.size());
      }
      for  (uint i=0;i<itsConstraints.size();i++) {
        itsMultiDirSolver.add_constraint(itsConstraints[i].get());
      }

      const size_t nDir = itsDirections.size();

      itsBufs.resize(itsSolInt);

      itsPredictSteps.resize(nDir);

      itsResultSteps.resize(nDir);
      for (size_t dir=0; dir<nDir; ++dir) {
        itsResultSteps[dir] = MultiResultStep::ShPtr(new MultiResultStep(itsSolInt));
        itsPredictSteps[dir]=Predict(input, parset, prefix, itsDirections[dir]);
        itsPredictSteps[dir].setNextStep(itsResultSteps[dir]);
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

      for (size_t dir=0; dir<nDir; ++dir) {
        itsPredictSteps[dir].updateInfo(infoIn);
      }

      if (itsSolInt==0) {
        itsSolInt=info().ntime();
      }

      if (itsNChan==0) {
        itsNChan = info().nchan();
      }
      if (itsNChan>info().nchan()) {
        itsNChan=info().nchan();
      }
      itsNFreqCells = info().nchan() / itsNChan;
      if (itsNChan*itsNFreqCells<info().nchan()) { // If last freq cell is smaller
        itsNFreqCells++;
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

      std::vector<std::string> sourceNames(itsDirections.size());
      std::vector<std::pair<double, double> > sourcePositions(itsDirections.size());
      for (uint i=0; i<itsDirections.size(); ++i) {
        sourceNames[i]=itsDirections[i][0]; // This only gives the name of the first patch
        sourcePositions[i]=itsPredictSteps[i].getFirstDirection();
      }
      itsH5Parm.addSources(sourceNames, sourcePositions);

      uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
      itsSols.resize(nSolTimes);
      itsNIter.resize(nSolTimes);
      itsConstraintSols.resize(nSolTimes);

      vector<double> chanFreqs(info().nchan());  //nChannelBlocks
      for (uint i=0;i<info().chanFreqs().size();++i) {
        chanFreqs[i]=info().chanFreqs()[i];
      }

      for (uint i=0; i<itsConstraints.size();++i) {
        itsConstraints[i]->init(
            info().antennaNames().size(),
            itsDirections.size(),
            info().nchan(), //nChannelBlocks
            &(chanFreqs[0])
        );
        
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
          coreConstraint->setCoreAntennas(coreAntennaIndices);
        }
        
        ScreenConstraint* screenConstraint = dynamic_cast<ScreenConstraint*>(itsConstraints[i].get());
        if(screenConstraint != 0)
        {
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
      }

      uint nSt = info().antennaNames().size();
      itsMultiDirSolver.init(nSt, nDir, info().nchan(), ant1, ant2);
      for (uint i=0; i<nSolTimes; ++i) {
        itsSols[i].resize(info().nchan());
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
      for (uint i=0; i<itsPredictSteps.size(); ++i) {
        itsPredictSteps[i].show(os);
      }
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

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerWrite.getElapsed(), totaltime);
      os << " of it spent in writing gain solutions to disk" << endl;

      os << "Iterations taken: [";
      for (uint i=0; i<itsNIter.size()-1; ++i) {
        os<<itsNIter[i]<<",";
      }
      os<<itsNIter[itsNIter.size()-1]<<"]"<<endl;
    }

    void DDECal::initializeSolutions() {
      if (itsTimeStep/itsSolInt>0 && itsPropagateSolutions) {
        // initialize solutions with those of the previous step
        itsSols[itsTimeStep/itsSolInt] = itsSols[itsTimeStep/itsSolInt-1];
      } else {
        // initialize solutions with 1.
        for (vector<vector<DComplex> >::iterator solveciter = itsSols[itsTimeStep/itsSolInt].begin();
             solveciter != itsSols[itsTimeStep/itsSolInt].end(); ++solveciter) {
           (*solveciter).assign(itsDirections.size()*info().antennaNames().size(), 1.0);
        }
      }
    }

    bool DDECal::process (const DPBuffer& bufin)
    {
      itsTimer.start();

      itsBufs[itsStepInSolInt].copy(bufin);
      itsDataPtrs[itsStepInSolInt] = itsBufs[itsStepInSolInt].getData().data();

      itsInput->fetchUVW(bufin, itsBufs[itsStepInSolInt], itsTimer);
      itsInput->fetchWeights(bufin, itsBufs[itsStepInSolInt], itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBufs[itsStepInSolInt], itsTimer);

      itsTimerPredict.start();

#pragma omp parallel for
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
        initializeSolutions();

        itsTimerSolve.start();
        MultiDirSolver::SolveResult solveResult = 
                  itsMultiDirSolver.process(itsDataPtrs, itsModelDataPtrs,
                  itsSols[itsTimeStep/itsSolInt],
                  itsAvgTime / itsSolInt);
        itsTimerSolve.stop();

        itsNIter[itsTimeStep/itsSolInt] = solveResult.iterations;
        // Store constraint solutions
        if (itsMode!=GainCal::COMPLEXGAIN && itsMode!=GainCal::PHASEONLY) {
          itsConstraintSols[itsTimeStep/itsSolInt]=solveResult._results;
        }

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

      uint nDir = itsDirections.size();

      vector<double> weights(info().nchan()*info().nantenna()*info().ntime()*itsDirections.size(),1.);

      if (itsConstraintSols[0].empty()) {
        // Record the actual iterands of the solver
        vector<DComplex> sols(info().nchan()*info().nantenna()*info().ntime()*nDir);
        size_t i=0;

        // Put solutions in a contiguous piece of memory
        for (uint dir=0; dir<nDir; ++dir) {
          for (uint ant=0; ant<info().nantenna(); ++ant) {
            for (uint chan=0; chan<info().nchan(); ++chan) {
              for (uint time=0; time<itsSols.size(); ++time) {
                ASSERT(!itsSols[time].empty());
                sols[i] = itsSols[time][chan][ant*nDir+dir];
                ++i;
              }
            }
          }
        }
        string axesnames="dir,ant,freq,time";
        vector<hsize_t> dims(4);
        dims[0]=nDir;
        dims[1]=info().nantenna();
        dims[2]=info().nchan();
        dims[3]=itsSols.size();

        uint numsols=(itsMode==GainCal::COMPLEXGAIN?2:1);
        for (uint solnum=0; solnum<numsols; ++solnum) {
          string solTabName;
          if (itsMode==GainCal::PHASEONLY) {
            solTabName = "phaseonly000";
            itsH5Parm.addSolution(solTabName, "scalarphase", axesnames,
                                  dims, sols, weights, false);
          } else if (itsMode==GainCal::COMPLEXGAIN) {
            if (solnum==0) {
              solTabName = "scalarphase000";
              itsH5Parm.addSolution(solTabName, "scalarphase", axesnames,
                                    dims, sols, weights, false);
            } else {
              solTabName = "scalaramplitude000";
              itsH5Parm.addSolution(solTabName, "scalaramplitude", axesnames,
                                    dims, sols, weights, true);
            }
          } else {
            THROW(Exception, "Constraint should have produced output");
          }
          // Tell H5Parm that all antennas and directions were used 
          // TODO: do this more cleanly
          std::vector<std::string> antennaNames(info().antennaNames().size());
          for (uint i=0; i<info().antennaNames().size(); ++i) {
            antennaNames[i]=info().antennaNames()[i];
          }
          itsH5Parm.setSolAntennas(solTabName, antennaNames);
    
          std::vector<std::string> sourceNames(itsDirections.size());
          for (uint i=0; i<itsDirections.size(); ++i) {
            sourceNames[i]=itsDirections[i][0]; // This only gives the name of the first patch
          }
          itsH5Parm.setSolSources(solTabName, sourceNames);
   
          // Set channel to frequency of middle channel 
          vector<double> chanFreqs(info().nchan());
          for (uint chan=0; chan<info().nchan(); ++chan) {
            chanFreqs[chan] = info().chanFreqs()[chan];
          }
          itsH5Parm.setFreqs(solTabName, chanFreqs);
    
          uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
          vector<double> solTimes(nSolTimes);
          double starttime=info().startTime();
          for (uint t=0; t<nSolTimes; ++t) {
            solTimes[t] = starttime+t*info().timeInterval()*itsSolInt+0.5*info().timeInterval();
          }
          // End TODO
          itsH5Parm.setTimes(solTabName, solTimes);
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

            vector<hsize_t> dims(firstResult.dims.size()+1); // Add time dimension
            dims[dims.size()-1]=itsConstraintSols.size();
            size_t numSols=dims[dims.size()-1];
            for (uint i=0; i<dims.size()-1; ++i) {
              dims[i] = firstResult.dims[i]; 
              numSols *= dims[i];
            }

            string axesnames = firstResult.axes + ",time";

            // Put solutions in a contiguous piece of memory
            vector<double> sols(numSols);
            size_t posInFlatSol=0;
            for (uint i=0; i<firstResult.vals.size(); ++i) {
              for (uint time=0; time<itsSols.size(); ++time) {
                sols[posInFlatSol++] = 
                  itsConstraintSols[time][constraintNum][solNameNum].vals[i];
              }
            }

            string solTabName = firstResult.name+"000";
            itsH5Parm.addSolution(solTabName, firstResult.name, 
                                  axesnames, dims, sols, weights);

            // Tell H5Parm that all antennas and directions were used 
            // TODO: do this more cleanly
            std::vector<std::string> antennaNames(info().antennaNames().size());
            for (uint i=0; i<info().antennaNames().size(); ++i) {
              antennaNames[i]=info().antennaNames()[i];
            }
            itsH5Parm.setSolAntennas(solTabName, antennaNames);
      
            std::vector<std::string> sourceNames(itsDirections.size());
            for (uint i=0; i<itsDirections.size(); ++i) {
              sourceNames[i]=itsDirections[i][0]; // This only gives the name of the first patch
            }
            itsH5Parm.setSolSources(solTabName, sourceNames);
     
            // Set channel to frequency of middle channel 
            vector<double> oneFreq(1);
            oneFreq[0] = info().chanFreqs()[info().nchan()/2];
            itsH5Parm.setFreqs(solTabName, oneFreq);
      
            uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
            vector<double> solTimes(nSolTimes);
            double starttime=info().startTime();
            for (uint t=0; t<nSolTimes; ++t) {
              solTimes[t] = starttime+t*info().timeInterval()*itsSolInt+0.5*info().timeInterval();
            }
            itsH5Parm.setTimes(solTabName, solTimes);
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
         initializeSolutions();

        //shrink itsDataPtrs, itsModelDataPtrs
        std::vector<casacore::Complex*>(itsDataPtrs.begin(),
            itsDataPtrs.begin()+itsStepInSolInt).swap(itsDataPtrs);
        std::vector<std::vector<casacore::Complex*> >(itsModelDataPtrs.begin(),
                    itsModelDataPtrs.begin()+itsStepInSolInt).swap(itsModelDataPtrs);
        itsTimerSolve.start();
        itsMultiDirSolver.process(itsDataPtrs, itsModelDataPtrs,
                                  itsSols[itsTimeStep/itsSolInt],
                                  itsAvgTime/itsStepInSolInt);
        itsTimerSolve.stop();
      }

      writeSolutions();

      itsTimer.stop();

      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
