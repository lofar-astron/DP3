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

#include "GainCal.h"
#include "Simulate.h"
#include "ApplyCal.h"
#include "PhaseFitter.h"
#include "CursorUtilCasa.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "DPLogger.h"
#include "Exceptions.h"
#include "SourceDBUtil.h"
#include "MSReader.h"
#include "Version.h"

#include "../ParmDB/ParmDB.h"
#include "../ParmDB/ParmValue.h"
#include "../ParmDB/SourceDB.h"

#include "../Common/ParallelFor.h"
#include "../Common/ParameterSet.h"
#include "../Common/StringUtil.h"

#include <fstream>
#include <ctime>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/OS/File.h>

#include <vector>
#include <algorithm>

#include <limits>
#include <iostream>
#include <iomanip>

using namespace casacore;
using namespace DP3::BBS;

namespace DP3 {
  namespace DPPP {

    GainCal::GainCal (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
      : itsInput         (input),
        itsName          (prefix),
        itsUseModelColumn(parset.getBool (prefix + "usemodelcolumn", false)),
        itsParmDBName    (parset.getString (prefix + "parmdb", "")),
        itsUseH5Parm     (itsParmDBName.find(".h5") != string::npos),
        itsDebugLevel    (parset.getInt (prefix + "debuglevel", 0)),
        itsDetectStalling (parset.getBool (prefix + "detectstalling", true)),
        itsApplySolution (parset.getBool (prefix + "applysolution", false)),
        itsUVWFlagStep   (input, parset, prefix),
        itsBaselineSelection (parset, prefix),
        itsMaxIter       (parset.getInt (prefix + "maxiter", 50)),
        itsTolerance     (parset.getDouble (prefix + "tolerance", 1.e-5)),
        itsPropagateSolutions
                         (parset.getBool(prefix + "propagatesolutions", true)),
        itsSolInt        (parset.getInt(prefix + "solint", 1)),
        itsNFreqCells    (0),
        itsConverged     (0),
        itsNonconverged  (0),
        itsFailed        (0),
        itsStalled       (0),
        itsStepInParmUpdate      (0),
        itsChunkStartTime(0),
        itsStepInSolInt        (0),
        itsAllSolutions ()
    {
      stringstream ss;
      ss << parset;
      itsParsetString = ss.str();

      if (itsParmDBName=="") {
        itsParmDBName=parset.getString("msin")+"/instrument";
      }

      if (!itsUseH5Parm) {
        itsTimeSlotsPerParmUpdate = parset.getInt(prefix +
                                                  "timeslotsperparmupdate",
                                                  500);
      } else {
        itsTimeSlotsPerParmUpdate = 0;
      }

      itsDataResultStep = ResultStep::ShPtr(new ResultStep());
      itsUVWFlagStep.setNextStep(itsDataResultStep);

      if (!itsUseModelColumn) {
        itsPredictStep.reset(new Predict(input, parset, prefix));
        itsResultStep = ResultStep::ShPtr(new ResultStep());
        itsPredictStep->setNextStep(itsResultStep);
      } else {
#ifdef HAVE_LOFAR_BEAM
        itsApplyBeamToModelColumn=parset.getBool(prefix +
                                              "applybeamtomodelcolumn", false);
        if (itsApplyBeamToModelColumn) {
          itsApplyBeamStep=ApplyBeam(input, parset, prefix, true);
          assert(!itsApplyBeamStep.invert());
          itsResultStep = ResultStep::ShPtr(new ResultStep());
          itsApplyBeamStep.setNextStep(itsResultStep);
        }
#endif
      }

      itsNIter.resize(4,0);

      if (itsApplySolution) {
        itsBuf.resize(itsSolInt);
      } else {
        itsBuf.resize(1);
      }

      string modestr = parset.getString (prefix + "caltype");
      itsMode = stringToCalType(modestr);
      uint defaultNChan = 0;
      if (itsMode == TECANDPHASE || itsMode == TEC) {
        defaultNChan = 1;
      }
      itsNChan = parset.getInt(prefix + "nchan", defaultNChan);
      assert(itsMode!=TECSCREEN);
    }

    GainCal::~GainCal()
    {}

    GainCal::CalType GainCal::stringToCalType(const string &modestr) {
      if (modestr=="diagonal"||modestr=="complexgain") return COMPLEXGAIN;
      else if (modestr=="scalarcomplexgain") return SCALARCOMPLEXGAIN;
      else if (modestr=="fulljones") return FULLJONES;
      else if (modestr=="phaseonly") return PHASEONLY;
      else if (modestr=="scalarphase") return SCALARPHASE;
      else if (modestr=="amplitudeonly") return AMPLITUDEONLY;
      else if (modestr=="scalaramplitude") return SCALARAMPLITUDE;
      else if (modestr=="tecandphase") return TECANDPHASE;
      else if (modestr=="tec") return TEC;
      else if (modestr=="tecscreen") return TECSCREEN;
      else if (modestr=="rotation+diagonal") return ROTATIONANDDIAGONAL;
      else if (modestr=="rotation") return ROTATION;
      throw Exception("Unknown mode: " + modestr);
    }

    string GainCal::calTypeToString(GainCal::CalType caltype) {
      switch(caltype)
      {
        case COMPLEXGAIN: return "complexgain";
        case SCALARCOMPLEXGAIN: return "scalarcomplexgain";
        case FULLJONES: return "fulljones";
        case PHASEONLY: return "phaseonly";
        case SCALARPHASE: return "scalarphase";
        case AMPLITUDEONLY: return "amplitudeonly";
        case SCALARAMPLITUDE: return "scalaramplitude";
        case TECANDPHASE: return "tecandphase";
        case TEC: return "tec";
        case TECSCREEN: return "tecscreen";
        case ROTATION: return "rotation";
        case ROTATIONANDDIAGONAL: return "rotation+diagonal";
        default: throw Exception("Unknown caltype: " + std::to_string(caltype));
      }
    }

    void GainCal::setAntennaUsed() {
      Matrix<bool> selbl(itsBaselineSelection.apply (info()));
      uint nBl=info().getAnt1().size();
      itsAntennaUsed.resize(info().antennaNames().size());
      itsAntennaUsed=false;
      for (uint bl=0; bl<nBl; ++bl) {
        if (selbl(info().getAnt1()[bl], info().getAnt2()[bl])) {
          itsAntennaUsed[info().getAnt1()[bl]] = true;
          itsAntennaUsed[info().getAnt2()[bl]] = true;
        }
      }
    }

    void GainCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();

      itsUVWFlagStep.updateInfo(infoIn);

      if (itsUseModelColumn) {
#ifdef HAVE_LOFAR_BEAM
        if (itsApplyBeamToModelColumn) {
          itsApplyBeamStep.updateInfo(infoIn);
        }
#endif
      } else {
        itsPredictStep->updateInfo(infoIn);
      }
      if (itsApplySolution) {
        info().setWriteData();
        info().setWriteFlags();
      }

      if (itsSolInt==0) {
        itsSolInt=info().ntime();
      }
      if (itsTimeSlotsPerParmUpdate==0) {
        itsTimeSlotsPerParmUpdate = info().ntime();
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

      itsSelectedBL = itsBaselineSelection.applyVec(info());
      setAntennaUsed();

      // Compute average frequency for every freqcell
      itsFreqData.resize(itsNFreqCells);
      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        double meanfreq=0;
        uint chmin=itsNChan*freqCell;
        uint chmax=min(info().nchan(), chmin+itsNChan);

        meanfreq = std::accumulate(info().chanFreqs().data()+chmin,
                                   info().chanFreqs().data()+chmax, 0.0);

        itsFreqData[freqCell] = meanfreq / (chmax-chmin);
      }

      // Initialize phase fitters, set their frequency data
      if (itsMode==TECANDPHASE || itsMode==TEC) {
        itsTECSols.reserve(itsTimeSlotsPerParmUpdate);

        itsPhaseFitters.reserve(itsNFreqCells); // TODO: could be numthreads instead

        uint nSt=info().antennaUsed().size();
        for (uint st=0; st<nSt; ++st) {
          itsPhaseFitters.push_back(std::unique_ptr<PhaseFitter>(new PhaseFitter(itsNFreqCells)));
          double* nu = itsPhaseFitters[st]->FrequencyData();
          for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
            nu[freqCell] = itsFreqData[freqCell];
          }
        }
      }

      iS.reserve(itsNFreqCells);
      uint chMax = itsNChan;
      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        if ((freqCell+1)*itsNChan>info().nchan()) { // Last cell can be smaller
          chMax-=((freqCell+1)*itsNChan)%info().nchan();
        }

        StefCal::StefCalMode smode;
        switch (itsMode)
        {
        case COMPLEXGAIN: smode = StefCal::DEFAULT; break;
        case FULLJONES: smode = StefCal::FULLJONES; break;
        case SCALARPHASE:
        case PHASEONLY:
        case TEC:
        case TECANDPHASE: smode = StefCal::PHASEONLY; break;
        case AMPLITUDEONLY:
        case SCALARAMPLITUDE: smode = StefCal::AMPLITUDEONLY; break;
        default: throw Exception("Unhandled mode");
        }

        iS.push_back(StefCal(itsSolInt, chMax, smode, scalarMode(itsMode),
                             itsTolerance, info().antennaUsed().size(),
                             itsDetectStalling, itsDebugLevel));
      }

      itsFlagCounter.init(getInfo());

      itsChunkStartTime = info().startTime();

      if (itsDebugLevel>0) {
        assert(getInfo().nThreads()==1);
        assert(itsTimeSlotsPerParmUpdate >= info().ntime());
        itsAllSolutions.resize(IPosition(6,
                               iS[0].numCorrelations(),
                               info().antennaUsed().size(),
                               (itsMode==TEC||itsMode==TECANDPHASE?2:1),
                               itsNFreqCells,
                               itsMaxIter,
                               info().ntime()
                              ));
      }
    }

    void GainCal::show (std::ostream& os) const
    {
      os << "GainCal " << itsName << endl;
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
      os << endl;
      os << "  solint:              " << itsSolInt <<endl;
      os << "  nchan:               " << itsNChan <<endl;
      os << "  max iter:            " << itsMaxIter << endl;
      os << "  tolerance:           " << itsTolerance << endl;
      os << "  caltype:             " << calTypeToString(itsMode) << endl;
      os << "  apply solution:      " << boolalpha << itsApplySolution << endl;
      os << "  propagate solutions: " << boolalpha << itsPropagateSolutions << endl;
      if (!itsUseH5Parm) {
        os << "  timeslotsperparmupdate: " << itsTimeSlotsPerParmUpdate << endl;
      }
      os << "  detect stalling:     " << boolalpha << itsDetectStalling << endl;
      os << "  use model column:    " << boolalpha << itsUseModelColumn << endl;
      itsBaselineSelection.show (os);
      if (!itsUseModelColumn) {
        itsPredictStep->show(os);
      }
#ifdef HAVE_LOFAR_BEAM
      else if (itsApplyBeamToModelColumn) {
        itsApplyBeamStep.show(os);
      }
#endif
      itsUVWFlagStep.show(os);
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

      if (itsMode == TEC || itsMode == TECANDPHASE) {
        os << "          ";
        FlagCounter::showPerc1 (os, itsTimerPhaseFit.getElapsed(), totaltime);
        os << " of it spent in fitting phases" << endl;
      }

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerWrite.getElapsed(), totaltime);
      os << " of it spent in writing gain solutions to disk" << endl;

      os << "        ";
      os <<"Converged: "<<itsConverged<<", stalled: "<<itsStalled<<", non converged: "<<itsNonconverged<<", failed: "<<itsFailed<<endl;
      os << "        ";
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
        // Need to keep a copy of all solint buffers in this step
        bufIndex=itsStepInSolInt;
        itsBuf[bufIndex].copy(bufin);
      } else {
        // We'll read the necessary info from the buffer and pass it on
        itsBuf[bufIndex].referenceFilled (bufin);
      }
      itsInput->fetchUVW(bufin, itsBuf[bufIndex], itsTimer);
      itsInput->fetchWeights(bufin, itsBuf[bufIndex], itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBuf[bufIndex], itsTimer);

      // UVW flagging happens on a copy of the buffer, so these flags are not written
      itsUVWFlagStep.process(itsBuf[bufIndex]);

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
#ifdef HAVE_LOFAR_BEAM
        if (itsApplyBeamToModelColumn) { // TODO: double check this
          // Temporarily put model data in data column for applybeam step
          // ApplyBeam step will copy the buffer so no harm is done
          itsBuf[bufIndex].getData()=itsModelData;
          itsApplyBeamStep.process(itsBuf[bufIndex]);
          //Put original data back in data column
          itsBuf[bufIndex].getData()=dataCube;
        }
#endif
      } else { // Predict
        itsPredictStep->process(itsBuf[bufIndex]);
      }

      itsTimerPredict.stop();

      itsTimerFill.start();

      if (itsStepInSolInt==0) {
        // Start new solution interval

        for (uint freqCell=0; freqCell<itsNFreqCells; freqCell++) {
          iS[freqCell].clearStationFlagged();
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

      if (!itsUseH5Parm && (itsStepInParmUpdate == itsTimeSlotsPerParmUpdate)) {
        writeSolutionsParmDB(itsChunkStartTime);
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
          uint antA = info().antennaMap()[info().getAnt1()[bl]];
          uint antB = info().antennaMap()[info().getAnt2()[bl]];
          uint freqCell = chan / itsNChan;
          if (nCr>2) {
            ApplyCal::applyFull( &invsol(0, antA, freqCell),
                       &invsol(0, antB, freqCell),
                       &data[bl * 4 * nchan + chan * 4 ],
                       &weight[bl * 4 * nchan + chan * 4 ], // Not passing weights, any pointer should do
                       &flag[  bl * 4 * nchan + chan * 4 ],
                       bl, chan, false, itsFlagCounter); // Update weights is disabled here
          }
          else if (scalarMode(itsMode)) {
            ApplyCal::applyScalar( &invsol(0, antA, freqCell),
                       &invsol(0, antB, freqCell),
                       &data[bl * 4 * nchan + chan * 4 ],
                       &weight[bl * 4 * nchan + chan * 4 ], // Not passing weights, any pointer should do
                       &flag[  bl * 4 * nchan + chan * 4 ],
                       bl, chan, false, itsFlagCounter); // Update weights is disabled here
          } else {
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
    void GainCal::fillMatrices (casacore::Complex* model, casacore::Complex* data, float* weight,
                                const casacore::Bool* flag) {
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = info().ncorr();
      assert(nCr==4 || nCr==2 || nCr==1);

      for (uint ch=0;ch<nCh;++ch) {
        for (uint bl=0;bl<nBl;++bl) {
          if (itsSelectedBL[bl]) {
            int ant1=info().antennaMap()[info().getAnt1()[bl]];
            int ant2=info().antennaMap()[info().getAnt2()[bl]];
            assert(ant1>=0 && ant2>=0);
            if (ant1==ant2 ||
                iS[ch/itsNChan].getStationFlagged()[ant1] ||
                iS[ch/itsNChan].getStationFlagged()[ant2] ||
                flag[bl*nCr*nCh+ch*nCr]) { // Only check flag of cr==0
              continue;
            }

            if (itsMode==TEC || itsMode==TECANDPHASE) {
              iS[ch/itsNChan].incrementWeight(weight[bl*nCr*nCh+ch*nCr]);
            }

            for (uint cr=0;cr<nCr;++cr) {
              // The nCrDiv is there such that for nCr==2 the visibilities end up at (0,0) for cr==0, (1,1) for cr==1
              uint nCrDiv = (nCr==4?2:1);
              iS[ch/itsNChan].getVis() (IPosition(6,ant1,cr/nCrDiv,itsStepInSolInt,ch%itsNChan,cr%2,ant2)) =
                  DComplex(data [bl*nCr*nCh+ch*nCr+cr]) *
                  DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
              iS[ch/itsNChan].getMVis()(IPosition(6,ant1,cr/nCrDiv,itsStepInSolInt,ch%itsNChan,cr%2,ant2)) =
                  DComplex(model[bl*nCr*nCh+ch*nCr+cr]) *
                  DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));

              // conjugate transpose
              iS[ch/itsNChan].getVis() (IPosition(6,ant2,cr%2,itsStepInSolInt,ch%itsNChan,cr/nCrDiv,ant1)) =
                  DComplex(conj(data [bl*nCr*nCh+ch*nCr+cr])) *
                  DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
              iS[ch/itsNChan].getMVis()(IPosition(6,ant2,cr%2,itsStepInSolInt,ch%itsNChan,cr/nCrDiv,ant1)) =
                  DComplex(conj(model[bl*nCr*nCh+ch*nCr+cr] )) *
                  DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
            }
          }
        }
      }
    }

    bool GainCal::scalarMode(CalType caltype) {
      return (caltype==TECANDPHASE || caltype==TEC || caltype==SCALARPHASE ||
              caltype==SCALARAMPLITUDE);
    }

    bool GainCal::diagonalMode(CalType caltype) {
      return (caltype==COMPLEXGAIN || caltype==PHASEONLY ||
              caltype==AMPLITUDEONLY);
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

      casacore::Matrix<double> tecsol(itsMode==TECANDPHASE?2:1,
                                  info().antennaUsed().size(), 0);

      vector<StefCal::Status> converged(itsNFreqCells,StefCal::NOTCONVERGED);
      for (;iter<itsMaxIter;++iter) {
        bool allConverged=true;
        ParallelFor<size_t> loop(getInfo().nThreads());
        loop.Run(0, itsNFreqCells, [&](size_t freqCell, size_t /*thread*/) {
          // Do another step when stalled and not all converged
          if (converged[freqCell]!=StefCal::CONVERGED)
          {
            converged[freqCell] = iS[freqCell].doStep(iter);
            // Only continue if there are steps worth continuing
            // (so not converged, failed or stalled)
            if (converged[freqCell]==StefCal::NOTCONVERGED) {
              allConverged = false;
            }
          }
        });

        if (itsDebugLevel>0) {
          for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
            Matrix<DComplex> fullSolution = iS[freqCell].getSolution(false);
            std::copy(fullSolution.begin(),
                      fullSolution.end(),
                      &(itsAllSolutions(IPosition(6, 0,
                                                     0,
                                                     0,
                                                     freqCell,
                                                     iter,
                                                     itsStepInParmUpdate
                                                 ))));
          }
        }

        if (itsMode==TEC || itsMode==TECANDPHASE) {
          itsTimerSolve.stop();
          itsTimerPhaseFit.start();
          casacore::Matrix<casacore::DComplex> sols_f(itsNFreqCells, info().antennaUsed().size());

          uint nSt = info().antennaUsed().size();

          // TODO: set phase reference so something smarter that station 0
          for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
            casacore::Matrix<casacore::DComplex> sol = iS[freqCell].getSolution(false);
            if (iS[freqCell].getStationFlagged()[0]) {
              // If reference station flagged, flag whole channel
              for (uint st=0; st<info().antennaUsed().size(); ++st) {
                iS[freqCell].getStationFlagged()[st] = true;
              }
            } else {
              for (uint st=0; st<info().antennaUsed().size(); ++st) {
                sols_f(freqCell, st) = sol(st, 0)/sol(0, 0);
                assert(isFinite(sols_f(freqCell, st)));
              }
            }
          }

          ParallelFor<size_t> loop(getInfo().nThreads());
          loop.Run(0, nSt, [&](size_t st, size_t /*thread*/) {
            uint numpoints=0;
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
                  assert(isFinite(phases[freqCell]));
                }
                assert(iS[freqCell].getWeight()>0);
                weights[freqCell] = iS[freqCell].getWeight();
                numpoints++;
              }
            }

            for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
              assert(isFinite(phases[freqCell]));
            }

            if (numpoints>1) { // TODO: limit should be higher
              //cout<<"tecsol(0,"<<st<<")="<<tecsol(0,st)<<", tecsol(1,"<<st<<")="<<tecsol(1,st)<<endl;
              if (itsMode==TECANDPHASE) {
                itsPhaseFitters[st]->FitDataToTEC2Model(tecsol(0, st), tecsol(1,st));
              } else { // itsMode==TEC
                itsPhaseFitters[st]->FitDataToTEC1Model(tecsol(0, st));
              }
              // Update solution in stefcal
              for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
                assert(isFinite(phases[freqCell]));
                iS[freqCell].getSolution(false)(st, 0) = polar(1., phases[freqCell]);
              }
            } else {
              tecsol(0, st) = 0; //std::numeric_limits<double>::quiet_NaN();
              if (itsMode==TECANDPHASE) {
                tecsol(1, st) = 0; //std::numeric_limits<double>::quiet_NaN();
              }
            }

            if (itsDebugLevel>0) {
              for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
                Matrix<DComplex> fullSolution = iS[freqCell].getSolution(false);
                std::copy(fullSolution.begin(),
                          fullSolution.end(),
                          &(itsAllSolutions(IPosition(6, 0,
                                                         0,
                                                         1,
                                                         freqCell,
                                                         iter,
                                                         itsStepInParmUpdate
                                                     ))));
              }
            }
          });
          itsTimerPhaseFit.stop();
          itsTimerSolve.start();
        }

        if (allConverged) {
          break;
        }

      } // End niter

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        switch (converged[freqCell]) {
        case StefCal::CONVERGED: {itsConverged++; itsNIter[0]+=iter; break;}
        case StefCal::STALLED: {itsStalled++; itsNIter[1]+=iter; break;}
        case StefCal::NOTCONVERGED: {itsNonconverged++; itsNIter[2]+=iter; break;}
        case StefCal::FAILED: {itsFailed++; itsNIter[3]+=iter; break;}
        default:
          throw Exception("Unknown converged status");
        }
      }

      // Stefcal terminated (either by maxiter or by converging)

      Cube<DComplex> sol(iS[0].numCorrelations(), info().antennaUsed().size(), itsNFreqCells);

      uint transpose[2][4] = { { 0, 1, 0, 0 }, { 0, 2, 1, 3 } };

      uint nSt = info().antennaUsed().size();

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        casacore::Matrix<casacore::DComplex> tmpsol = iS[freqCell].getSolution(true);

        for (uint st=0; st<nSt; st++) {
          for (uint cr=0; cr<iS[0].nCr(); ++cr) {
            uint crt=transpose[iS[0].numCorrelations()/4][cr];  // Conjugate transpose ! (only for numCorrelations = 4)
            sol(crt, st, freqCell) = conj(tmpsol(st, cr));        // Conjugate transpose
            if (itsMode==COMPLEXGAIN || itsMode==PHASEONLY || itsMode==AMPLITUDEONLY) {
              sol(crt+1, st, freqCell) = conj(tmpsol(st+nSt,cr)); // Conjugate transpose
            }
          }
        }
      }
      itsSols.push_back(sol);
      if (itsMode==TEC || itsMode==TECANDPHASE) {
        itsTECSols.push_back(tecsol);
      }

      itsTimerSolve.stop();
    } // End stefcal()


    void GainCal::initParmDB() {
      itsParmDB = std::shared_ptr<BBS::ParmDB>
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
      if (itsMode==PHASEONLY || itsMode==SCALARPHASE) {
        itsParmDB->getDefValues(parmset, "Gain:0:0:Ampl");
        if (parmset.empty()) {
          ParmValueSet pvset(ParmValue(1.0));
          itsParmDB->putDefValue("Gain:0:0:Ampl",pvset);
          itsParmDB->putDefValue("Gain:1:1:Ampl",pvset);
        }
      }

      // Write out default phases
      if (itsMode==AMPLITUDEONLY || itsMode==SCALARAMPLITUDE) {
        itsParmDB->getDefValues(parmset, "Gain:0:0:Phase");
        if (parmset.empty()) {
          ParmValueSet pvset(ParmValue(0.0));
          itsParmDB->putDefValue("Gain:0:0:Phase",pvset);
          itsParmDB->putDefValue("Gain:1:1:Phase",pvset);
        }
      }

      // Write out default gains
      if (itsMode==COMPLEXGAIN || itsMode==FULLJONES) {
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
      if (itsMode==SCALARPHASE) {
        name=string("CommonScalarPhase:");
      } else if (itsMode==SCALARAMPLITUDE) {
        name=string("CommonScalarAmplitude:");
      } else if (itsMode==TEC || itsMode==TECANDPHASE) {
        name=string("TEC:");
      }
      else {
        name=string("Gain:");
      }

      return name;
    }

    void GainCal::writeSolutionsH5Parm(double) {
      itsTimer.start();
      itsTimerWrite.start();

      H5Parm h5parm(itsParmDBName, true);

      // Fill antenna info in H5Parm, need to convert from casa types to std types
      std::vector<std::string> allAntennaNames(info().antennaNames().size());
      std::vector<std::vector<double> > antennaPos(info().antennaPos().size());
      for (uint i=0; i<info().antennaNames().size(); ++i) {
        allAntennaNames[i]=info().antennaNames()[i];
        casacore::Quantum<casacore::Vector<double> > pos = info().antennaPos()[i].get("m");
        antennaPos[i].resize(3);
        antennaPos[i][0] = pos.getValue()[0];
        antennaPos[i][1] = pos.getValue()[1];
        antennaPos[i][2] = pos.getValue()[2];
      }

      h5parm.addAntennas(allAntennaNames, antennaPos);

      vector<pair<double, double> > pointingPosition(1);
      MDirection phasecenter = info().phaseCenter();
      pointingPosition[0].first = phasecenter.getValue().get()[0];
      pointingPosition[0].second = phasecenter.getValue().get()[1];
      vector<string> pointingName(1, "POINTING");

      h5parm.addSources(pointingName, pointingPosition);

      uint nPol;
      vector<string> polarizations;
      if (scalarMode(itsMode)) {
        nPol = 1;
      } else if (diagonalMode(itsMode)) {
        nPol = 2;
        polarizations.push_back("XX");
        polarizations.push_back("YY");
      } else {
        assert(itsMode==FULLJONES);
        polarizations.push_back("XX");
        polarizations.push_back("XY");
        polarizations.push_back("YX");
        polarizations.push_back("YY");
        nPol = 4;
      }

      // Construct time axis
      uint nSolTimes = (info().ntime()+itsSolInt-1)/itsSolInt;
      vector<double> solTimes(nSolTimes);
      assert(nSolTimes==itsSols.size());
      double starttime=info().startTime();
      for (uint t=0; t<nSolTimes; ++t) {
        solTimes[t] = starttime+(t+0.5)*info().timeInterval()*itsSolInt;
      }

      // Construct frequency axis
      uint nSolFreqs;
      if (itsMode==TEC || itsMode==TECANDPHASE) {
        nSolFreqs = 1;
      } else {
        nSolFreqs = itsNFreqCells;
      }

      vector<H5Parm::AxisInfo> axes;
      axes.push_back(H5Parm::AxisInfo("time", itsSols.size()));
      axes.push_back(H5Parm::AxisInfo("freq", nSolFreqs));
      axes.push_back(H5Parm::AxisInfo("ant", info().antennaUsed().size()));
      if (nPol>1) {
        axes.push_back(H5Parm::AxisInfo("pol", nPol));
      }

      vector<H5Parm::SolTab> soltabs = makeSolTab(h5parm, itsMode, axes);

      std::vector<std::string> antennaUsedNames;
      for (uint st = 0; st<info().antennaUsed().size(); ++st) {
          antennaUsedNames.push_back(info().antennaNames()[info().antennaUsed()[st]]);
      }

      vector<H5Parm::SolTab>::iterator soltabiter = soltabs.begin();
      for (; soltabiter != soltabs.end(); ++soltabiter) {
        (*soltabiter).setAntennas(antennaUsedNames);
        if (nPol>1) {
          (*soltabiter).setPolarizations(polarizations);
        }
        if (itsMode==TEC || itsMode==TECANDPHASE) {
          // Set channel to frequency of middle channel
          // TODO: fix this for nchan
          vector<double> oneFreq(1);
          oneFreq[0] = info().chanFreqs()[info().nchan()/2];
          (*soltabiter).setFreqs(oneFreq);
        } else {
          (*soltabiter).setFreqs(itsFreqData);
        }
        (*soltabiter).setTimes(solTimes);
      }

      // Put solutions in a contiguous piece of memory
      string historyString = "CREATE by DPPP\n" +
                  DPPPVersion::AsString() + "\n" +
                  "step " + itsName + " in parset: \n" + itsParsetString;

      if (itsMode==TEC || itsMode==TECANDPHASE) {
        vector<double> tecsols(nSolFreqs*antennaUsedNames.size()*nSolTimes*nPol);
        vector<double> weights(nSolFreqs*antennaUsedNames.size()*nSolTimes*nPol, 1.);
        vector<double> phasesols;
        if (itsMode==TECANDPHASE) {
          phasesols.resize(nSolFreqs*antennaUsedNames.size()*nSolTimes*nPol);
        }
        size_t i=0;
        for (uint time=0; time<nSolTimes; ++time) {
          for (uint freqCell=0; freqCell<nSolFreqs; ++freqCell) {
            for (uint ant=0; ant<info().antennaUsed().size(); ++ant) {
              for (uint pol=0; pol<nPol; ++pol) {
                assert(!itsTECSols[time].empty());
                tecsols[i] = itsTECSols[time](0, ant) / 8.44797245e9;
                if (!std::isfinite(tecsols[i])) {
                  weights[i] = 0.;
                }
                if (itsMode==TECANDPHASE) {
                  phasesols[i] = -itsTECSols[time](0, ant);
                }
                ++i;
              }
            }
          }
        }
        soltabs[0].setValues(tecsols, weights, historyString);
        if (itsMode==TECANDPHASE) {
          soltabs[1].setValues(phasesols, weights, historyString);
        }
      } else {
        vector<DComplex> sols(nSolFreqs*antennaUsedNames.size()*nSolTimes*nPol);
        vector<double> weights(nSolFreqs*antennaUsedNames.size()*nSolTimes*nPol, 1.);
        size_t i=0;
        for (uint time=0; time<nSolTimes; ++time) {
          for (uint freqCell=0; freqCell<nSolFreqs; ++freqCell) {
            for (uint ant=0; ant<info().antennaUsed().size(); ++ant) {
              for (uint pol=0; pol<nPol; ++pol) {
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

        if (itsMode!=AMPLITUDEONLY) {
          soltabs[0].setComplexValues(sols, weights, false, historyString);
        } else {
          soltabs[0].setComplexValues(sols, weights, true, historyString);
        }
        if (soltabs.size()>1) {
          // Also write amplitudes
          soltabs[1].setComplexValues(sols, weights, true, historyString);
        }
      }

      itsTimerWrite.stop();
      itsTimer.stop();
    }

    vector<H5Parm::SolTab> GainCal::makeSolTab(H5Parm& h5parm, CalType caltype,
                                               vector<H5Parm::AxisInfo>& axes) {
      uint numsols = 1;
      // For [scalar]complexgain, store two soltabs: phase and amplitude
      if (caltype == GainCal::COMPLEXGAIN ||
          caltype == GainCal::SCALARCOMPLEXGAIN ||
          caltype == GainCal::TECANDPHASE ||
          caltype == GainCal::FULLJONES) {
        numsols = 2;
      }
      vector<H5Parm::SolTab> soltabs;
      for (uint solnum=0; solnum<numsols; ++solnum) {
        string solTabName;
        H5Parm::SolTab soltab;
        switch (caltype) {
          case GainCal::SCALARPHASE:
          case GainCal::PHASEONLY:
            solTabName = "phase000";
            soltab = h5parm.createSolTab(solTabName, "phase", axes);
            break;
          case GainCal::SCALARCOMPLEXGAIN:
          case GainCal::COMPLEXGAIN:
          case GainCal::FULLJONES:
            if (solnum==0) {
              solTabName = "phase000";
              soltab = h5parm.createSolTab(solTabName, "phase", axes);
            } else {
              solTabName = "amplitude000";
              soltab = h5parm.createSolTab(solTabName, "amplitude", axes);
            }
            break;
          case GainCal::SCALARAMPLITUDE:
          case GainCal::AMPLITUDEONLY:
            solTabName = "amplitude000";
            soltab = h5parm.createSolTab(solTabName, "amplitude", axes);
            break;
          case GainCal::TEC:
          case GainCal::TECANDPHASE:
            if (solnum==0) {
              solTabName = "tec000";
              soltab = h5parm.createSolTab(solTabName, "tec", axes);
            } else {
              solTabName = "phase000";
              soltab = h5parm.createSolTab(solTabName, "phase", axes);
            }
            break;
          default:
            throw Exception("Unhandled mode in writing H5Parm output: " + calTypeToString(caltype));
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
      } // End initialization of parmdb

      uint ntime=itsSols.size();
      uint nchan, nfreqs;
      if (itsMode==TEC || itsMode==TECANDPHASE) {
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

      // Get end time of the current chunk. For the last chunk, this
      // is chopped off at the end of the MS (only if solint > 1)
      double endTime = min(startTime + ntime * info().timeInterval() * itsSolInt,
                           info().startTime() + info().ntime() * info().timeInterval());

      // Make time axis (can be non regular for last chunk if solint > 1)
      vector<double> lowtimes(ntime), hightimes(ntime);
      for (uint t=0; t<ntime; ++t) {
        lowtimes[t]  = startTime + info().timeInterval() * itsSolInt * t;
        hightimes[t] = min(startTime + info().timeInterval() * itsSolInt * (t+1),
                           endTime);
      }
      BBS::Axis::ShPtr timeAxis = Axis::makeAxis(lowtimes, hightimes);

      BBS::Axis::ShPtr freqAxis(
          new BBS::RegularAxis(
              getInfo().chanFreqs()[0] - freqWidth * 0.5,
              freqWidth*nchan,
              nfreqs));
      BBS::Grid solGrid(freqAxis, timeAxis);

      // Construct domain grid for the current chunk
      BBS::Axis::ShPtr tdomAxis(
          new BBS::RegularAxis(
              startTime,
              endTime - startTime,
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

      uint nSt=info().antennaUsed().size();

      for (size_t st=0; st<nSt; ++st) {
        // Do not write NaN solutions for stations that were not used
        if (!itsAntennaUsed[info().antennaUsed()[st]]) {
          // itsAntennaUsed is indexed with real antenna numbers, so antennaUsed() is needed
          continue;
        }
        for (int pol=0; pol<4; ++pol) { // For 0101
          if (scalarMode(itsMode) && pol>0) {
            continue;
          } else if (diagonalMode(itsMode) && (pol==1||pol==2)) {
            continue;
          }
          int realimmax; // For tecandphase, this functions as dummy between tec and commonscalarphase
          if (itsMode==PHASEONLY || itsMode==SCALARPHASE ||
              itsMode==AMPLITUDEONLY || itsMode==SCALARAMPLITUDE || itsMode==TEC) {
            realimmax=1;
          } else {
            realimmax=2;
          }
          for (int realim=0; realim<realimmax; ++realim) { // For real and imaginary
            string name = parmName();

            if (itsMode!=SCALARPHASE && itsMode!=SCALARAMPLITUDE) {
              name+=str0101[pol];
              if (itsMode==PHASEONLY) {
                name=name+"Phase:";
              } else if (itsMode==AMPLITUDEONLY) {
                name=name+"Ampl:";
              } else {
                name=name+strri[realim];
              }
            }
            if (itsMode==TECANDPHASE || itsMode==TEC) {
              if (realim==0) {
                name="TEC:";
              } else {
                name="CommonScalarPhase:";
              }
            }

            name+=info().antennaNames()[info().antennaUsed()[st]];

            // Collect its solutions for all times and frequency cells in a single array.
            for (uint ts=0; ts<ntime; ++ts) {
              for (uint freqCell=0; freqCell<nfreqs; ++freqCell) {
                if (itsMode==FULLJONES) {
                  if (realim==0) {
                    values(freqCell, ts) = real(itsSols[ts](pol,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](pol,st,freqCell));
                  }
                } else if (itsMode==COMPLEXGAIN) {
                  if (realim==0) {
                    values(freqCell, ts) = real(itsSols[ts](pol/3,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](pol/3,st,freqCell));
                  }
                } else if (itsMode==SCALARPHASE || itsMode==PHASEONLY) {
                  values(freqCell, ts) = arg(itsSols[ts](pol/3,st,freqCell));
                } else if (itsMode==SCALARAMPLITUDE || itsMode==AMPLITUDEONLY) {
                  values(freqCell, ts) = abs(itsSols[ts](pol/3,st,freqCell));
                } else if (itsMode==TEC || itsMode==TECANDPHASE) {
                  if (realim==0) {
                    values(freqCell, ts) = itsTECSols[ts](realim,st) / 8.44797245e9;
                  } else {
                    values(freqCell, ts) = -itsTECSols[ts](realim,st); // TODO: why is there a minus here?
                  }
                }
                else {
                  throw Exception("Unhandled mode");
                }
              }
            }
            BBS::ParmValue::ShPtr pv(new BBS::ParmValue());
            pv->setScalars (solGrid, values);

            BBS::ParmValueSet pvs(domainGrid,
                                  vector<BBS::ParmValue::ShPtr>(1, pv));
            std::map<std::string,int>::const_iterator pit = itsParmIdMap.find(name);

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
        if (itsUseH5Parm) {
          writeSolutionsH5Parm(itsChunkStartTime);
        } else {
          writeSolutionsParmDB(itsChunkStartTime);
        }
        if (itsDebugLevel>0) {
          H5::H5File hdf5file = H5::H5File("debug.h5", H5F_ACC_TRUNC);
          vector<hsize_t> dims(6);
          for (uint i=0; i<6; ++i) {
            dims[i] = itsAllSolutions.shape()[5-i];
          }
          H5::DataSpace dataspace(dims.size(), &(dims[0]), NULL);
          H5::CompType complex_data_type(sizeof(DComplex));
          complex_data_type.insertMember( "r", 0, H5::PredType::IEEE_F64LE);
          complex_data_type.insertMember( "i", sizeof(double), H5::PredType::IEEE_F64LE);
          H5::DataSet dataset = hdf5file.createDataSet("val",
                                                       complex_data_type,
                                                       dataspace);
          dataset.write(itsAllSolutions.data(), complex_data_type);
          hdf5file.close();
        }
      }

      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
