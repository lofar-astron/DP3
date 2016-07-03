//# GainCal.cc: DPPP step class to do a gain calibration
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/GainCal.h>
#include <DPPP/Simulate.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/phasefitter.h>
#include <DPPP/CursorUtilCasa.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/SourceDBUtil.h>
#include <DPPP/MSReader.h>
#include <DPPP/DPLogger.h>
#include <ParmDB/ParmDB.h>
#include <ParmDB/ParmValue.h>
#include <ParmDB/SourceDB.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>

#include <fstream>
#include <ctime>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>
#include <casa/OS/File.h>

#include <vector>
#include <algorithm>

#include <limits>
#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    GainCal::GainCal (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsUseModelColumn(parset.getBool (prefix + "usemodelcolumn", false)),
        itsParmDBName    (parset.getString (prefix + "parmdb", "")),
        itsMode          (parset.getString (prefix + "caltype")),
        itsDebugLevel    (parset.getInt (prefix + "debuglevel", 0)),
        itsDetectStalling (parset.getBool (prefix + "detectstalling", true)),
        itsApplySolution (parset.getBool (prefix + "applysolution", false)),
        itsBaselines     (),
        itsMaxIter       (parset.getInt (prefix + "maxiter", 50)),
        itsTolerance     (parset.getDouble (prefix + "tolerance", 1.e-5)),
        itsPropagateSolutions
                         (parset.getBool(prefix + "propagatesolutions", true)),
        itsSolInt        (parset.getInt(prefix + "solint", 1)),
        itsNChan         (parset.getInt(prefix + "nchan", 0)),
        itsNFreqCells    (0),
        itsMinBLperAnt   (parset.getInt(prefix + "minblperant", 4)),
        itsTimeSlotsPerParmUpdate
                         (parset.getInt(prefix + "timeslotsperparmupdate", 500)),
        itsConverged     (0),
        itsNonconverged  (0),
        itsFailed        (0),
        itsStalled       (0),
        itsStepInParmUpdate      (0),
        itsChunkStartTime(0),
        itsStepInSolInt        (0)
    {
      if (itsParmDBName=="") {
        itsParmDBName=parset.getString("msin")+"/instrument";
      }

      if (!itsUseModelColumn) {
        itsPredictStep=Predict(input, parset, prefix);
        itsResultStep=new ResultStep();
        itsPredictStep.setNextStep(DPStep::ShPtr(itsResultStep));
      } else {
        itsApplyBeamToModelColumn=parset.getBool(prefix +
                                              "applybeamtomodelcolumn", false);
        if (itsApplyBeamToModelColumn) {
          itsApplyBeamStep=ApplyBeam(input, parset, prefix, true);
          ASSERT(!itsApplyBeamStep.invert());
          itsResultStep=new ResultStep();
          itsApplyBeamStep.setNextStep(DPStep::ShPtr(itsResultStep));
        }
      }

      itsNIter.resize(4,0);

      if (itsApplySolution) {
        itsBuf.resize(itsSolInt);
      } else {
        itsBuf.resize(1);
      }

      ASSERT(itsMode=="diagonal" || itsMode=="phaseonly" ||
             itsMode=="fulljones" || itsMode=="scalarphase" ||
             itsMode=="amplitudeonly" || itsMode=="scalaramplitude" ||
             itsMode=="tec");
    }

    GainCal::~GainCal()
    {}

    void GainCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();

      const size_t nBl=info().nbaselines();

      if (itsUseModelColumn) {
        if (itsApplyBeamToModelColumn) {
          itsApplyBeamStep.updateInfo(infoIn);
        }
      } else {
        itsPredictStep.updateInfo(infoIn);
      }
      if (itsApplySolution) {
        info().setWriteData();
        info().setWriteFlags();
      }

      for (uint i=0; i<nBl; ++i) {
        itsBaselines.push_back (Baseline(info().getAnt1()[i],
                                         info().getAnt2()[i]));
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

      itsSols.reserve(itsTimeSlotsPerParmUpdate);

      // Initialize phase fitters, set their frequency data
      if (itsMode=="tec") {
        itsTECSols.reserve(itsTimeSlotsPerParmUpdate);

        itsPhaseFitters.reserve(itsNFreqCells); // TODO: could be numthreads instead
        vector<double> freqData(itsNFreqCells);
        for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
          double meanfreq=0;
          uint chmin=itsNChan*freqCell;
          uint chmax=min(info().nchan(), chmin+itsNChan);

          meanfreq = std::accumulate(info().chanFreqs().data()+chmin,
                                     info().chanFreqs().data()+chmax, 0.0);

          freqData[freqCell] = meanfreq / (chmax-chmin);
        }

        uint nSt=info().antennaNames().size();
        for (uint st=0; st<nSt; ++st) {
          itsPhaseFitters.push_back(CountedPtr<PhaseFitter>(new PhaseFitter(itsNFreqCells)));
          double* nu = itsPhaseFitters[st]->FrequencyData();
          for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
            nu[freqCell] = freqData[freqCell];
          }
        }
      }

      iS.reserve(itsNFreqCells);
      uint chMax = itsNChan;
      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        if ((freqCell+1)*itsNChan>info().nchan()) { // Last cell can be smaller
          chMax-=((freqCell+1)*itsNChan)%info().nchan();
        }
        iS.push_back(StefCal(itsSolInt, chMax, itsMode,
                             itsTolerance, info().antennaNames().size(),
                             itsDetectStalling, itsDebugLevel));
      }

      itsFlagCounter.init(getInfo());

      itsChunkStartTime = info().startTime();
    }

    void GainCal::show (std::ostream& os) const
    {
      os << "GainCal " << itsName << endl;
      os << "  parmdb:             " << itsParmDBName;
      if (Table::isReadable(itsParmDBName)) {
        os << " (existing)";
      } else {
        os << " (will be created)";
      }
      os << endl;
      os << "  solint:              " << itsSolInt <<endl;
      os << "  nchan:               " << itsNChan <<endl;
      os << "  max iter:            " << itsMaxIter << endl;
      os << "  tolerance:           " << itsTolerance << endl;
      os << "  caltype:             " << itsMode << endl;
      os << "  apply solution:      " << boolalpha << itsApplySolution << endl;
      os << "  propagate solutions: " << boolalpha << itsPropagateSolutions << endl;
      os << "  detect stalling:     " << boolalpha << itsDetectStalling << endl;
      os << "  use model column:    " << boolalpha << itsUseModelColumn << endl;
      if (!itsUseModelColumn) {
        itsPredictStep.show(os);
      } else if (itsApplyBeamToModelColumn) {
        itsApplyBeamStep.show(os);
      }
    }

    void GainCal::showTimings (std::ostream& os, double duration) const
    {
      double totaltime=itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " GainCal " << itsName << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerPredict.getElapsed(), totaltime);
      os << " of it spent in predict" << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerFill.getElapsed(), totaltime);
      os << " of it spent in reordering visibility data" << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerSolve.getElapsed(), totaltime);
      os << " of it spent in estimating gains and computing residuals" << endl;

      if (itsMode == "tec") {
        os << "          ";
        FlagCounter::showPerc1 (os, itsTimerPhaseFit.getElapsed(), totaltime);
        os << "of it spent in fitting phases" << endl;
      }

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerWrite.getElapsed(), totaltime);
      os << " of it spent in writing gain solutions to disk" << endl;

      os << "          ";
      os <<"Converged: "<<itsConverged<<", stalled: "<<itsStalled<<", non converged: "<<itsNonconverged<<", failed: "<<itsFailed<<endl;
      os <<"Iters converged: " << (itsConverged==0?0:itsNIter[0]/itsConverged);
      os << ", stalled: "<<      (itsStalled  ==0?0:itsNIter[1]/itsStalled);
      os << ", non converged: "<<(itsNonconverged==0?0:itsNIter[2]/itsNonconverged);
      os << ", failed: "<<(itsFailed==0?0:itsNIter[3]/itsFailed)<<endl;
    }

    bool GainCal::process (const DPBuffer& bufin)
    {
      itsTimer.start();

      uint bufIndex=0;

      if (itsApplySolution) {
        bufIndex=itsStepInSolInt;
        itsBuf[bufIndex].copy(bufin);
      } else {
        itsBuf[0].referenceFilled (bufin);
      }
      itsInput->fetchUVW(bufin, itsBuf[bufIndex], itsTimer);
      itsInput->fetchWeights(bufin, itsBuf[bufIndex], itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBuf[bufIndex], itsTimer);

      Cube<Complex> dataCube=itsBuf[bufIndex].getData();
      Complex* data=dataCube.data();
      float* weight = itsBuf[bufIndex].getWeights().data();
      const Bool* flag=itsBuf[bufIndex].getFlags().data();

      // Simulate.
      //
      // Model visibilities for each direction of interest will be computed
      // and stored.

      itsTimerPredict.start();

      if (itsUseModelColumn) {
        itsInput->getModelData (itsBuf[bufIndex].getRowNrs(), itsModelData);
        if (itsApplyBeamToModelColumn) { // TODO: double check this
          // Temporarily put model data in data column for applybeam step
          // ApplyBeam step will copy the buffer so no harm is done
          itsBuf[bufIndex].getData()=itsModelData;
          itsApplyBeamStep.process(itsBuf[bufIndex]);
          //Put original data back in data column
          itsBuf[bufIndex].getData()=dataCube;
        }
      } else { // Predict
        itsPredictStep.process(itsBuf[bufIndex]);
      }

      itsTimerPredict.stop();

      itsTimerFill.start();

      if (itsStepInSolInt==0) {
        // Start new solution interval

        for (uint freqCell=0; freqCell<itsNFreqCells; freqCell++) {
          setAntennaMaps(flag, freqCell);
          iS[freqCell].resetVis();
        }
      }

      // Store data in the stefcal object
      if (itsUseModelColumn && !itsApplyBeamToModelColumn) {
        fillMatrices(itsModelData.data(),data,weight,flag);
      } else {
        fillMatrices(itsResultStep->get().getData().data(),data,weight,flag);
      }
      itsTimerFill.stop();

      if (itsStepInSolInt==itsSolInt-1) {
        // Solve past solution interval
        stefcal();
        itsStepInParmUpdate++;

        if (itsApplySolution) {
          Cube<DComplex> invsol = invertSol(itsSols.back());
          for (uint stepInSolInt=0; stepInSolInt<itsSolInt; stepInSolInt++) {
            applySolution(itsBuf[stepInSolInt], invsol);
            getNextStep()->process(itsBuf[stepInSolInt]);
          }
        }

        itsStepInSolInt=0;
      } else {
        itsStepInSolInt++;
      }

      itsTimer.stop();

      if (itsStepInParmUpdate == itsTimeSlotsPerParmUpdate) {
        writeSolutions(itsChunkStartTime);
        itsChunkStartTime += itsSolInt * itsTimeSlotsPerParmUpdate * info().timeInterval();
        itsSols.clear();
        itsTECSols.clear();
        itsStepInParmUpdate = 0;
      }

      if (!itsApplySolution) {
        getNextStep()->process(itsBuf[bufIndex]);
      }
      return false;
    }

    Cube<DComplex> GainCal::invertSol(const Cube<DComplex>& sol) {
      Cube<DComplex> invsol = sol.copy();
      uint nCr = invsol.shape()[0];

      // Invert copy of solutions
      uint nSt = invsol.shape()[1];
      for (uint st=0; st<nSt; ++st) {
        for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
          if (nCr==4) {
            ApplyCal::invert(&invsol(0,st,freqCell));
          } else {
            for (uint cr=0; cr<nCr; ++cr) {
              invsol(cr, st, freqCell) = 1./invsol(cr, st, freqCell);
            }
          }
        }
      }

      return invsol;
    }

    void GainCal::applySolution(DPBuffer& buf, const Cube<DComplex>& invsol) {
      uint nbl = buf.getData().shape()[2];
      Complex* data = buf.getData().data();
      float* weight = buf.getWeights().data(); // Not initialized yet
      bool* flag = buf.getFlags().data();
      uint nchan = buf.getData().shape()[1];

      uint nCr = invsol.shape()[0];

      for (size_t bl=0; bl<nbl; ++bl) {
        for (size_t chan=0;chan<nchan;chan++) {
          uint antA = info().getAnt1()[bl];
          uint antB = info().getAnt2()[bl];
          uint freqCell = chan / itsNChan;
          if (nCr>2) {
            ApplyCal::applyFull( &invsol(0, antA, freqCell),
                       &invsol(0, antB, freqCell),
                       &data[bl * 4 * nchan + chan * 4 ],
                       &weight[bl * 4 * nchan + chan * 4 ], // Not passing weights, any pointer should do
                       &flag[  bl * 4 * nchan + chan * 4 ],
                       bl, chan, false, itsFlagCounter); // Update weights is disabled here
          }
          else {
            ApplyCal::applyDiag( &invsol(0, antA, freqCell),
                       &invsol(0, antB, freqCell),
                       &data[bl * 4 * nchan + chan * 4 ],
                       &weight[bl * 4 * nchan + chan * 4 ], // Not passing weights, any pointer should do
                       &flag[  bl * 4 * nchan + chan * 4 ],
                       bl, chan, false, itsFlagCounter); // Update weights is disabled here
          }
        }
      }
    }

    // Fills itsVis and itsMVis as matrices with all 00 polarizations in the
    // top left, all 11 polarizations in the bottom right, etc.
    // For TEC fitting, it also sets weights for the frequency cells
    void GainCal::fillMatrices (casa::Complex* model, casa::Complex* data, float* weight,
                                const casa::Bool* flag) {
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;

      for (uint ch=0;ch<nCh;++ch) {
        for (uint bl=0;bl<nBl;++bl) {
          int ant1=info().getAnt1()[bl];
          int ant2=info().getAnt2()[bl];
          if (ant1==ant2 ||
              iS[ch/itsNChan].getStationFlagged()[ant1] ||
              iS[ch/itsNChan].getStationFlagged()[ant2] ||
              flag[bl*nCr*nCh+ch*nCr]) { // Only check flag of cr==0
            continue;
          }

          if (itsMode=="tec") {
            iS[ch/itsNChan].incrementWeight(weight[bl*nCr*nCh+ch*nCr]);
          }

          for (uint cr=0;cr<nCr;++cr) {
            iS[ch/itsNChan].getVis() (IPosition(6,ant1,cr/2,itsStepInSolInt,ch%itsNChan,cr%2,ant2)) =
                DComplex(data [bl*nCr*nCh+ch*nCr+cr]) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
            iS[ch/itsNChan].getMVis()(IPosition(6,ant1,cr/2,itsStepInSolInt,ch%itsNChan,cr%2,ant2)) =
                DComplex(model[bl*nCr*nCh+ch*nCr+cr]) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));

            // conjugate transpose
            iS[ch/itsNChan].getVis() (IPosition(6,ant2,cr%2,itsStepInSolInt,ch%itsNChan,cr/2,ant1)) =
                DComplex(conj(data [bl*nCr*nCh+ch*nCr+cr])) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
            iS[ch/itsNChan].getMVis()(IPosition(6,ant2,cr%2,itsStepInSolInt,ch%itsNChan,cr/2,ant1)) =
                DComplex(conj(model[bl*nCr*nCh+ch*nCr+cr] )) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
          }
        }
      }
    }

    void GainCal::setAntennaMaps (const Bool* flag, uint freqCell) {
      uint nCr=info().ncorr();
      uint nCh=info().nchan();
      uint nBl=info().nbaselines();

      iS[freqCell].clearStationFlagged();

      casa::Vector<casa::uInt> dataPerAntenna; // nAnt
      dataPerAntenna.resize(info().antennaNames().size());

      // TODO: implement smarter graph algorithm here that requires
      // less than O(n_ant^3) worst case operations.
      bool flaggedAntsAdded=true;
      uint loopcount=0;
      while (flaggedAntsAdded) {
        dataPerAntenna=0;
        for (uint bl=0;bl<nBl;++bl) {
          uint ant1=info().getAnt1()[bl];
          uint ant2=info().getAnt2()[bl];
          if (ant1==ant2) {
            continue;
          }
          uint chmax=min((freqCell+1)*itsNChan, nCh);
          for (uint ch=freqCell*itsNChan;ch<chmax;++ch) {
            for (uint cr=0;cr<nCr;++cr) {
              if (!(flag[bl*nCr*nCh + ch*nCr + cr]
                    || iS[freqCell].getStationFlagged()[ant1]
                    || iS[freqCell].getStationFlagged()[ant2]
                   )) {
                dataPerAntenna(ant1)++;
                dataPerAntenna(ant2)++;
              }
            }
          }
        }

        flaggedAntsAdded=false;
        for (uint ant=0; ant<info().antennaNames().size(); ++ant) {
          if (dataPerAntenna(ant)>nCr*itsMinBLperAnt) {
            iS[freqCell].getStationFlagged()[ant]=false; // Index in stefcal numbering
          } else { // Not enough data
            //cout<<"flagging station "<<ant<<", "<<dataPerAntenna(ant)<<endl;
            if (!iS[freqCell].getStationFlagged()[ant]) {
              iS[freqCell].getStationFlagged()[ant]=true;
              flaggedAntsAdded=true;
            }
          }
        }
        loopcount++;
      }
    }

    void GainCal::stefcal () {
      itsTimerSolve.start();

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        if (itsPropagateSolutions) {
          iS[freqCell].init(false);
        } else {
          iS[freqCell].init(true);
        }
      }

      uint iter=0;

      casa::Matrix<double> tecsol(2, info().antennaNames().size(), 0);

      std::vector<StefCal::Status> converged(itsNFreqCells,StefCal::NOTCONVERGED);
      for (;iter<itsMaxIter;++iter) {
        bool allConverged=true;
#pragma omp parallel for
        for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
          if (converged[freqCell]==StefCal::CONVERGED) { // Do another step when stalled and not all converged
            continue;
          }
          converged[freqCell] = iS[freqCell].doStep(iter);
          if (converged[freqCell]==StefCal::NOTCONVERGED) { // Only continue if there are steps worth continuing (so not converged, failed or stalled)
            allConverged = false;
          } 
        }

        if (itsMode=="tec") {
          itsTimerSolve.stop();
          itsTimerPhaseFit.start();
          casa::Matrix<casa::DComplex> sols_f(itsNFreqCells, info().antennaNames().size());

          uint nSt = info().antennaNames().size();

          // TODO: set phase reference so something smarter that station 0
          for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
            casa::Matrix<casa::DComplex> sol = iS[freqCell].getSolution(false);
            if (iS[freqCell].getStationFlagged()[0]) {
              // If reference station flagged, flag whole channel
              for (uint st=0; st<info().antennaNames().size(); ++st) {
                iS[freqCell].getStationFlagged()[st] = true;
              }
            } else {
              for (uint st=0; st<info().antennaNames().size(); ++st) {
                sols_f(freqCell, st) = sol(st, 0)/sol(0, 0);
                ASSERT(isFinite(sols_f(freqCell, st)));
              }
            }
          }

#pragma omp parallel for
          for (uint st=0; st<nSt; ++st) {
            uint numpoints=0;
            double cost=0;
            double* phases = itsPhaseFitters[st]->PhaseData();
            double* weights = itsPhaseFitters[st]->WeightData();
            for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
              if (iS[freqCell].getStationFlagged()[st%nSt] || 
                  converged[freqCell]==StefCal::FAILED) {
                phases[freqCell] = 0;
                weights[freqCell] = 0;
              } else {
                phases[freqCell] = arg(sols_f(freqCell, st));
                if (!isFinite(phases[freqCell])) {
                  cout<<"Yuk, phases[freqCell]="<<phases[freqCell]<<", sols_f(freqCell, st)="<<sols_f(freqCell, st)<<endl;
                  ASSERT(isFinite(phases[freqCell]));
                }
                ASSERT(iS[freqCell].getWeight()>0);
                weights[freqCell] = iS[freqCell].getWeight();
                if (itsDebugLevel > 0 && st==34) {
                  cout<<"w["<<freqCell<<"]="<<weights[freqCell]<<endl;
                }
                numpoints++;
              }
            }

            if (itsDebugLevel>0) {
              cout<<"st="<<st<<", numpoints="<<numpoints<<endl;
            }

            if (itsDebugLevel>0 && st==34) {
              cout<<"t="<<itsStepInParmUpdate<<", st="<<st<<", unfitted["<<iter<<"]=[";
              uint freqCell2=0;
              for (; freqCell2<itsNFreqCells-1; ++freqCell2) {
                cout<<phases[freqCell2]<<",";
              }
              cout<<phases[freqCell2]<<"];"<<endl;
            }

            for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
              ASSERT(isFinite(phases[freqCell]));
            }

            if (numpoints>1) { // TODO: limit should be higher
              //cout<<"tecsol(0,"<<st<<")="<<tecsol(0,st)<<", tecsol(1,"<<st<<")="<<tecsol(1,st)<<endl;
              cost=itsPhaseFitters[st]->FitDataToTECModel(tecsol(0, st), tecsol(1,st));
              // Update solution in stefcal
              for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
                ASSERT(isFinite(phases[freqCell]));
                iS[freqCell].getSolution(false)(st, 0) = polar(1., phases[freqCell]);
              }
            } else {
              tecsol(0, st) = 0; //std::numeric_limits<double>::quiet_NaN();
              tecsol(1, st) = 0; //std::numeric_limits<double>::quiet_NaN();
            }

            if (st==34 && itsDebugLevel>0) {
              cout<<"fitted["<<st<<"]=[";
              uint freqCell=0;
              for (; freqCell<itsNFreqCells-1; ++freqCell) {
                cout<<phases[freqCell]<<",";
              }
              cout<<phases[freqCell]<<"]"<<endl;
              cout << "fitdata["<<st<<"]=[" << tecsol(0,1) << ", " << tecsol(1,1) << ", " << cost << "];" << endl;
            }
          }
          itsTimerPhaseFit.stop();
          itsTimerSolve.start();
        }

        if (allConverged) {
          break;
        }

#ifdef DEBUG
        if (itsDebugLevel>1) { // Only antenna 1
          cout<<"phases["<<iter<<"]=[";
          uint freqCell=0;
          for (; freqCell<itsNFreqCells-1; ++freqCell) {
            cout<<arg(iS[freqCell].getSolution(false)(1, 0)/iS[freqCell].getSolution(false)(0,0))<<",";
          }
          cout<<arg(iS[freqCell].getSolution(false)(1, 0)/iS[freqCell].getSolution(false)(0,0))<<"];"<<endl;
        }
        if (itsDebugLevel>2) { // Statistics about nStalled
          cout<<"convstatus["<<iter<<"]=[";
          uint freqCell=0;
          for (; freqCell<itsNFreqCells-1; ++freqCell) {
            cout<<converged[freqCell]<<",";
          }
          cout<<converged[freqCell]<<"]"<<endl;
        }
#endif
      } // End niter

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        switch (converged[freqCell]) {
        case StefCal::CONVERGED: {itsConverged++; break;}
        case StefCal::STALLED: {itsStalled++; itsNIter[1]+=iter; break;}
        case StefCal::NOTCONVERGED: {itsNonconverged++; itsNIter[2]+=iter; break;}
        case StefCal::FAILED: {itsFailed++; itsNIter[3]+=iter; break;}
        default:
          THROW(Exception, "Unknown converged status");
        }
      }

      // Stefcal terminated (either by maxiter or by converging)

      Cube<DComplex> sol(iS[0].numCorrelations(), info().antennaNames().size(), itsNFreqCells);

      uint transpose[2][4] = { { 0, 1, 0, 0 }, { 0, 2, 1, 3 } };

      uint nSt = info().antennaNames().size();

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        casa::Matrix<casa::DComplex> tmpsol = iS[freqCell].getSolution(true);

        for (uint st=0; st<nSt; st++) {
          for (uint cr=0; cr<iS[0].nCr(); ++cr) {
            uint crt=transpose[iS[0].numCorrelations()/4][cr];  // Conjugate transpose ! (only for numCorrelations = 4)
            sol(crt, st, freqCell) = conj(tmpsol(st, cr));        // Conjugate transpose
            if (itsMode=="diagonal" || itsMode=="phaseonly" || itsMode=="amplitudeonly") {
              sol(crt+1, st, freqCell) = conj(tmpsol(st+nSt,cr)); // Conjugate transpose
            }
          }
        }
      }
      itsSols.push_back(sol);
      if (itsMode=="tec") {
        itsTECSols.push_back(tecsol);
      }

      itsTimerSolve.stop();
    } // End stefcal()


    void GainCal::initParmDB() {
      itsParmDB = boost::shared_ptr<BBS::ParmDB>
        (new BBS::ParmDB(BBS::ParmDBMeta("casa", itsParmDBName),
                         false));
      itsParmDB->lock();
      // Store the (freq, time) resolution of the solutions.

      double freqWidth = getInfo().chanWidths()[0];
      if (getInfo().chanFreqs().size()>1) { // Handle data with evenly spaced gaps between channels
        freqWidth = info().chanFreqs()[1]-info().chanFreqs()[0];
      }

      vector<double> resolution(2);
      resolution[0] = freqWidth * itsNChan;
      resolution[1] = info().timeInterval() * itsSolInt;
      itsParmDB->setDefaultSteps(resolution);

      string parmname=parmName()+"*";

      if (!itsParmDB->getNames(parmname).empty()) {
        DPLOG_WARN_STR ("Solutions for "<<parmname<<" already in "<<itsParmDBName
                        <<", these are removed");
        // Specify entire domain of this MS; only to specify that the existing
        // values should be deleted for this domain
        BBS::Axis::ShPtr tdomAxis(
            new BBS::RegularAxis(
                info().startTime(),
                info().ntime() * info().timeInterval(),
                1));
        BBS::Axis::ShPtr fdomAxis(
            new BBS::RegularAxis(
                info().chanFreqs()[0] - freqWidth * 0.5,
                freqWidth * getInfo().chanFreqs().size(), 1));

        itsParmDB->deleteValues(parmname, BBS::Box(
                                    fdomAxis->start(), tdomAxis->start(),
                                    fdomAxis->end(), tdomAxis->end(), true));
      }

      // Write out default values, if they don't exist yet
      ParmMap parmset;

      // Write out default amplitudes
      if (itsMode=="phaseonly" || itsMode=="scalarphase") {
        itsParmDB->getDefValues(parmset, "Gain:0:0:Ampl");
        if (parmset.empty()) {
          ParmValueSet pvset(ParmValue(1.0));
          itsParmDB->putDefValue("Gain:0:0:Ampl",pvset);
          itsParmDB->putDefValue("Gain:1:1:Ampl",pvset);
        }
      }

      // Write out default phases
      if (itsMode=="amplitudeonly" || itsMode=="scalaramplitude") {
        itsParmDB->getDefValues(parmset, "Gain:0:0:Phase");
        if (parmset.empty()) {
          ParmValueSet pvset(ParmValue(0.0));
          itsParmDB->putDefValue("Gain:0:0:Phase",pvset);
          itsParmDB->putDefValue("Gain:1:1:Phase",pvset);
        }
      }

      // Write out default gains
      if (itsMode=="diagonal" || itsMode=="fulljones") {
        itsParmDB->getDefValues(parmset, "Gain:0:0:Real");
        if (parmset.empty()) {
          ParmValueSet pvset(ParmValue(1.0));
          itsParmDB->putDefValue("Gain:0:0:Real",pvset);
          itsParmDB->putDefValue("Gain:1:1:Real",pvset);
        }
      }
    }

    string GainCal::parmName() {
      string name;
      if (itsMode=="scalarphase") {
        name=string("CommonScalarPhase:");
      } else if (itsMode=="scalaramplitude") {
        name=string("CommonScalarAmplitude:");
      } else if (itsMode=="tec") {
        name=string("TEC:");
      }
      else {
        name=string("Gain:");
      }

      return name;
    }

    void GainCal::writeSolutions(double startTime) {
      itsTimer.start();
      itsTimerWrite.start();

      // Open the ParmDB at the first write.
      // In that way the instrumentmodel ParmDB can be in the MS directory.
      if (!itsParmDB) {
        initParmDB();
      } // End initialization of parmdb

      uint ntime=itsSols.size();
      uint nchan, nfreqs;
      if (itsMode=="tec") {
        nfreqs = 1;
        nchan = info().nchan();
      } else {
        nfreqs = itsNFreqCells;
        nchan = itsNChan;
      }

      // Construct solution grid for the current chunk
      double freqWidth = getInfo().chanWidths()[0];
      if (getInfo().chanFreqs().size()>1) { // Handle data with evenly spaced gaps between channels
        freqWidth = info().chanFreqs()[1]-info().chanFreqs()[0];
      }
      BBS::Axis::ShPtr freqAxis(
          new BBS::RegularAxis(
              getInfo().chanFreqs()[0] - freqWidth * 0.5,
              freqWidth*nchan,
              nfreqs));
      BBS::Axis::ShPtr timeAxis(
          new BBS::RegularAxis(
              startTime,
              info().timeInterval() * itsSolInt,
              ntime));
      BBS::Grid solGrid(freqAxis, timeAxis);

      // Construct domain grid for the current chunk
      BBS::Axis::ShPtr tdomAxis(
          new BBS::RegularAxis(
              startTime,
              ntime * info().timeInterval() * itsSolInt,
              1));
      BBS::Axis::ShPtr fdomAxis(
          new BBS::RegularAxis(
              info().chanFreqs()[0] - freqWidth * 0.5,
              freqWidth * getInfo().chanFreqs().size(), 1));
      BBS::Grid domainGrid(fdomAxis, tdomAxis);

      // Write the solutions per parameter.
      const char* str0101[] = {"0:0:","0:1:","1:0:","1:1:"};
      const char* strri[] = {"Real:","Imag:"};
      Matrix<double> values(nfreqs, ntime);

      DComplex sol;

      uint nSt=info().antennaNames().size();

      for (size_t st=0; st<nSt; ++st) {
        for (int pol=0; pol<4; ++pol) { // For 0101
          if ((itsMode=="diagonal" || itsMode=="phaseonly" ||
               itsMode=="amplitudeonly") && (pol==1||pol==2)) {
            continue;
          }
          if ((itsMode=="scalarphase" || itsMode=="scalaramplitude" ||
              itsMode=="tec") && pol>0) {
            continue;
          }
          int realimmax; // For tec, this functions as dummy between tec and commonscalarphase
          if (itsMode=="phaseonly" || itsMode=="scalarphase" ||
              itsMode=="amplitudeonly" || itsMode=="scalaramplitude") {
            realimmax=1;
          } else {
            realimmax=2;
          }
          for (int realim=0; realim<realimmax; ++realim) { // For real and imaginary
            string name = parmName();

            if (itsMode!="scalarphase" && itsMode!="scalaramplitude") {
              name+=str0101[pol];
              if (itsMode=="phaseonly") {
                name=name+"Phase:";
              } else if (itsMode=="amplitudeonly") {
                name=name+"Ampl:";
              } else {
                name=name+strri[realim];
              }
            }
            if (itsMode=="tec") {
              if (realim==0) {
                name="TEC:";
              } else {
                name="CommonScalarPhase:";
              }
            }

            name+=info().antennaNames()[st];

            // Collect its solutions for all times and frequency cells in a single array.
            for (uint ts=0; ts<ntime; ++ts) {
              for (uint freqCell=0; freqCell<nfreqs; ++freqCell) {
                if (itsMode=="fulljones") {
                  if (realim==0) {
                    values(freqCell, ts) = real(itsSols[ts](pol,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](pol,st,freqCell));
                  }
                } else if (itsMode=="diagonal") {
                  if (realim==0) {
                    values(freqCell, ts) = real(itsSols[ts](pol/3,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](pol/3,st,freqCell));
                  }
                } else if (itsMode=="scalarphase" || itsMode=="phaseonly") {
                  values(freqCell, ts) = arg(itsSols[ts](pol/3,st,freqCell));
                } else if (itsMode=="scalaramplitude" || itsMode=="amplitudeonly") {
                  values(freqCell, ts) = abs(itsSols[ts](pol/3,st,freqCell));
                } else if (itsMode=="tec") {
                  if (realim==0) {
                    values(freqCell, ts) = itsTECSols[ts](realim,st) / 8.44797245e9;
                  } else {
                    values(freqCell, ts) = -itsTECSols[ts](realim,st); // TODO: why is there a minus here?
                  }
                }
                else {
                  THROW (Exception, "Unhandled mode");
                }
              }
            }
            BBS::ParmValue::ShPtr pv(new BBS::ParmValue());
            pv->setScalars (solGrid, values);

            BBS::ParmValueSet pvs(domainGrid,
                                  vector<BBS::ParmValue::ShPtr>(1, pv));
            map<string,int>::const_iterator pit = itsParmIdMap.find(name);

            if (pit == itsParmIdMap.end()) {
              // First time, so a new nameId will be set.
              // Check if the name was defined in the parmdb previously
              int nameId = itsParmDB->getNameId(name);
              itsParmDB->putValues (name, nameId, pvs);
              itsParmIdMap[name] = nameId;
            } else {
              // Parm has been put before.
              int nameId = pit->second;
              itsParmDB->putValues (name, nameId, pvs);
            }
          }
        }
      }

      itsTimerWrite.stop();
      itsTimer.stop();
    }

    void GainCal::finish()
    {
      itsTimer.start();

      //Solve remaining time slots if any
      if (itsStepInSolInt!=0) {
        stefcal();

        if (itsApplySolution) {
          Cube<DComplex> invsol = invertSol(itsSols.back());
          for (uint stepInSolInt=0; stepInSolInt<itsStepInSolInt; stepInSolInt++) {
            applySolution(itsBuf[stepInSolInt], invsol);
            getNextStep()->process(itsBuf[stepInSolInt]);
          }
        }
      }

      itsTimer.stop();

      if (!itsSols.empty()) {
        writeSolutions(itsChunkStartTime);
      }

      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
