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
        itsBaselines     (),
        itsMaxIter       (parset.getInt (prefix + "maxiter", 50)),
        itsTolerance     (parset.getDouble (prefix + "tolerance", 1.e-5)),
        itsPropagateSolutions
                         (parset.getBool(prefix + "propagatesolutions", false)),
        itsSolInt        (parset.getInt(prefix + "solint", 1)),
        itsNChan         (parset.getInt(prefix + "nchan", 0)),
        itsNFreqCells    (0),
        itsMinBLperAnt   (parset.getInt(prefix + "minblperant", 4)),
        itsTimeSlotsPerParmUpdate
                         (parset.getInt(prefix + "timeslotsperparmupdate", 500)),
        itsConverged     (0),
        itsNonconverged  (0),
        itsStalled       (0),
        itsTimeStep      (0),
        itsChunkStartTime(0),
        itsNTimes        (0)
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
      ASSERT(itsMode=="diagonal" || itsMode=="phaseonly" ||
             itsMode=="fulljones" || itsMode=="scalarphase" ||
             itsMode=="amplitudeonly" || itsMode=="scalaramplitude");
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
      info().setWriteData();
      info().setWriteFlags();

      for (uint i=0; i<nBl; ++i) {
        itsBaselines.push_back (Baseline(info().getAnt1()[i],
                                         info().getAnt2()[i]));
      }

      if (itsSolInt==0) {
        itsSolInt=info().ntime();
      }

      itsSols.reserve(itsTimeSlotsPerParmUpdate);

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

      itsAntennaUsedNames.resize(info().antennaUsed().size());
      for (int ant=0, nAnts=info().antennaUsed().size(); ant<nAnts; ++ant) {
        itsAntennaUsedNames[ant]=info().antennaNames()[info().antennaUsed()[ant]];
      }

      iS.reserve(itsNFreqCells);
      uint chMax = itsNChan;
      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        if ((freqCell+1)*itsNChan>info().nchan()) { // Last cell can be smaller
          chMax-=((freqCell+1)*itsNChan)%info().nchan();
        }
        iS.push_back(StefCal(itsSolInt, chMax, itsMode, itsTolerance,
                             info().antennaNames().size(), itsDetectStalling,
                             itsDebugLevel));
      }

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
      os << "  solint:             " << itsSolInt <<endl;
      os << "  nchan:              " << itsNChan <<endl;
      os << "  max iter:           " << itsMaxIter << endl;
      os << "  tolerance:          " << itsTolerance << endl;
      os << "  mode:               " << itsMode << endl;
      os << "  detect stalling:    " << boolalpha << itsDetectStalling << endl;
      os << "  use model column:   " << boolalpha << itsUseModelColumn << endl;
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

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerWrite.getElapsed(), totaltime);
      os << " of it spent in writing gain solutions to disk" << endl;

      os << "          ";
      os <<"Converged: "<<itsConverged<<", stalled: "<<itsStalled<<", non converged: "<<itsNonconverged<<endl;
    }

    bool GainCal::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuf.referenceFilled (bufin);
      itsInput->fetchUVW(bufin, itsBuf, itsTimer);
      itsInput->fetchWeights(bufin, itsBuf, itsTimer);
      itsInput->fetchFullResFlags(bufin, itsBuf, itsTimer);

      Cube<Complex> dataCube=itsBuf.getData();
      Complex* data=dataCube.data();
      float* weight = itsBuf.getWeights().data();
      const Bool* flag=itsBuf.getFlags().data();

      // Simulate.
      //
      // Model visibilities for each direction of interest will be computed
      // and stored.

      itsTimerPredict.start();

      if (itsUseModelColumn) {
        itsInput->getModelData (itsBuf.getRowNrs(), itsModelData);
        if (itsApplyBeamToModelColumn) {
          // Temporarily put model data in data column for applybeam step
          // ApplyBeam step will copy the buffer so no harm is done
          itsBuf.getData()=itsModelData;
          itsApplyBeamStep.process(itsBuf);
          //Put original data back in data column
          itsBuf.getData()=dataCube;
        }
      } else { // Predict
        itsPredictStep.process(itsBuf);
      }

      itsTimerPredict.stop();

      itsTimerFill.start();

      if (itsNTimes==0) {
        // Start new solution interval

        for (uint freqCell=0; freqCell<itsNFreqCells; freqCell++) {
          uint nSt=setAntennaMaps(flag, freqCell);
          //cout<<"nSt="<<nSt<<endl;
          iS[freqCell].resetVis(nSt);
        }
      }

      // Store data in the stefcal object
      if (itsUseModelColumn && !itsApplyBeamToModelColumn) {
        fillMatrices(itsModelData.data(),data,weight,flag);
      } else {
        fillMatrices(itsResultStep->get().getData().data(),data,weight,flag);
      }
      itsTimerFill.stop();

      if (itsNTimes==itsSolInt-1) {
        // Solve past solution interval
        stefcal();
        itsTimeStep++;
        itsNTimes=0;
      } else {
        itsNTimes++;
      }

      itsTimer.stop();

      if (itsTimeStep == itsTimeSlotsPerParmUpdate) {
        writeSolutions(itsChunkStartTime);
        itsChunkStartTime += itsSolInt * itsTimeSlotsPerParmUpdate * info().timeInterval();
        itsSols.clear();
        itsTimeStep = 0;
      }
      getNextStep()->process(itsBuf);
      return false;
    }

    // Fills itsVis and itsMVis as matrices with all 00 polarizations in the
    // top left, all 11 polarizations in the bottom right, etc.
    void GainCal::fillMatrices (casa::Complex* model, casa::Complex* data, float* weight,
                                const casa::Bool* flag) {
      const size_t nBl = info().nbaselines();
      const size_t nCh = info().nchan();
      const size_t nCr = 4;

      for (uint ch=0;ch<nCh;++ch) {
        for (uint bl=0;bl<nBl;++bl) {
          int ant1=iS[ch/itsNChan].getAntMap()[info().getAnt1()[bl]];
          int ant2=iS[ch/itsNChan].getAntMap()[info().getAnt2()[bl]];
          if (ant1==ant2 || ant1==-1 || ant2 == -1 || flag[bl*nCr*nCh+ch*nCr]) { // Only check flag of cr==0
            continue;
          }

          for (uint cr=0;cr<nCr;++cr) {
            iS[ch/itsNChan].getVis() (IPosition(6,ant1,cr/2,itsNTimes,ch%itsNChan,cr%2,ant2)) =
                DComplex(data [bl*nCr*nCh+ch*nCr+cr]) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
            iS[ch/itsNChan].getMVis()(IPosition(6,ant1,cr/2,itsNTimes,ch%itsNChan,cr%2,ant2)) =
                DComplex(model[bl*nCr*nCh+ch*nCr+cr]) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));

            // conjugate transpose
            iS[ch/itsNChan].getVis() (IPosition(6,ant2,cr%2,itsNTimes,ch%itsNChan,cr/2,ant1)) =
                DComplex(conj(data [bl*nCr*nCh+ch*nCr+cr])) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
            iS[ch/itsNChan].getMVis()(IPosition(6,ant2,cr%2,itsNTimes,ch%itsNChan,cr/2,ant1)) =
                DComplex(conj(model[bl*nCr*nCh+ch*nCr+cr] )) *
                DComplex(sqrt(weight[bl*nCr*nCh+ch*nCr+cr]));
          }
        }
      }
    }

    uint GainCal::setAntennaMaps (const Bool* flag, uint freqCell) {
      uint nCr=info().ncorr();
      uint nCh=info().nchan();
      uint nBl=info().nbaselines();

      casa::Vector<casa::uInt> dataPerAntenna; // nAnt
      dataPerAntenna.resize(info().antennaNames().size());
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
            if (!flag[bl*nCr*nCh + ch*nCr + cr]) {
              dataPerAntenna(ant1)++;
              dataPerAntenna(ant2)++;
            }
          }
        }
      }

      uint nSt=0;

      for (uint ant=0; ant<info().antennaNames().size(); ++ant) {
        if (dataPerAntenna(ant)>nCr*itsMinBLperAnt) {
          iS[freqCell].getAntMap()[ant]=nSt++; // Index in stefcal numbering
        } else {
          iS[freqCell].getAntMap()[ant]=-1; // Not enough data
        }
      }

      return nSt;
    }

    void GainCal::stefcal () {
      itsTimerSolve.start();

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        iS[freqCell].init();
      }

      uint iter=0;

      std::vector<StefCal::Status> converged(itsNFreqCells,StefCal::NOTCONVERGED);
      for (;iter<itsMaxIter;++iter) {
        bool allConverged=true;
        for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
          if (converged[freqCell]==StefCal::CONVERGED) { // Do another step when stalled and not all converged
            continue;
          }
          converged[freqCell] = iS[freqCell].doStep(iter);
          if (converged[freqCell]==StefCal::NOTCONVERGED) {
            allConverged = false;
          }
          if (allConverged) {
            break;
          }
        }       } // End niter

      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        switch (converged[freqCell]) {
        case StefCal::CONVERGED: {itsConverged++; break;}
        case StefCal::STALLED: {itsStalled++; break;}
        case StefCal::NOTCONVERGED: {itsStalled++; break;}
        default:
          THROW(Exception, "Unknown converged status");
        }
      }

      // Stefcal terminated (either by maxiter or by converging)
      // Let's save G...
      Cube<DComplex> allg(iS[0].numCorrelations(), info().antennaNames().size(), itsNFreqCells);

      uint transpose[2][4] = { { 0, 1, 0, 0 }, { 0, 2, 1, 3 } };
      for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
        //cout<<endl<<"freqCell="<<freqCell<<", timeCell="<<itsTStep/itsSolInt<<", tstep="<<itsTStep<<endl;
        casa::Matrix<casa::DComplex> sol = iS[freqCell].getSolution();
        for (uint st=0; st<info().antennaNames().size(); st++) {
          for (uint cr=0; cr<iS[0].numCorrelations(); ++cr) {
            uint crt=transpose[iS[0].numCorrelations()/4][cr]; // Conjugate transpose ! (only for numCorrelations = 4)
            allg(crt, st, freqCell) = conj(sol(st, cr));
          }
        }
      }
      itsSols.push_back(allg);

      itsTimerSolve.stop();
    }

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
      string name=(itsMode=="commonscalarphase"?"CommonScalarPhase:*":"Gain:*");
      if (!itsParmDB->getNames(name).empty()) {
        DPLOG_WARN_STR ("Solutions for "<<name<<" already in "<<itsParmDBName
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

        itsParmDB->deleteValues(name, BBS::Box(
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

    void GainCal::writeSolutions(double startTime) {
      itsTimer.start();
      itsTimerWrite.start();

      // Open the ParmDB at the first write.
      // In that way the instrumentmodel ParmDB can be in the MS directory.
      if (! itsParmDB) {
        initParmDB();
      } // End initialization of parmdb

      uint ntime=itsSols.size();

      // Construct solution grid for the current chunk
      double freqWidth = getInfo().chanWidths()[0];
      if (getInfo().chanFreqs().size()>1) { // Handle data with evenly spaced gaps between channels
        freqWidth = info().chanFreqs()[1]-info().chanFreqs()[0];
      }
      BBS::Axis::ShPtr freqAxis(
          new BBS::RegularAxis(
              getInfo().chanFreqs()[0] - freqWidth * 0.5,
              freqWidth*itsNChan,
              itsNFreqCells));
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
      Matrix<double> values(itsNFreqCells, ntime);

      DComplex sol;

      uint nSt=info().antennaNames().size();
      for (size_t st=0; st<nSt; ++st) {
        uint seqnr = 0; // To take care of real and imaginary part
        string suffix(info().antennaNames()[st]);

        for (int pol=0; pol<4; ++pol) { // For 0101
          if ((itsMode=="diagonal" || itsMode=="phaseonly" ||
               itsMode=="amplitudeonly") && (pol==1||pol==2)) {
            continue;
          }
          if ((itsMode=="scalarphase" || itsMode=="scalaramplitude" ) && pol>0) {
            continue;
          }
          int realimmax;
          if (itsMode=="phaseonly" || itsMode=="scalarphase" ||
              itsMode=="amplitudeonly" || itsMode=="scalaramplitude") {
            realimmax=1;
          } else {
            realimmax=2;
          }
          for (int realim=0; realim<realimmax; ++realim) { // For real and imaginary
            string name;

            if (itsMode=="scalarphase") {
              name=string("CommonScalarPhase:")+suffix;
            } else if (itsMode=="scalaramplitude") {
              name=string("CommonScalarAmplitude:")+suffix;
            } else {
              name=string("Gain:") + str0101[pol];
              if (itsMode=="phaseonly") {
                name=name+"Phase:";
              } else if (itsMode=="amplitudeonly") {
                name=name+"Ampl:";
              } else {
                name=name+strri[realim];
              }
              name+=suffix;
            }
            // Collect its solutions for all times and frequency cells in a single array.
            for (uint ts=0; ts<ntime; ++ts) {
              for (uint freqCell=0; freqCell<itsNFreqCells; ++freqCell) {
                if (itsMode=="fulljones") {
                  if (seqnr%2==0) {
                    values(freqCell, ts) = real(itsSols[ts](seqnr/2,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](seqnr/2,st,freqCell));
                  }
                } else if (itsMode=="diagonal") {
                  if (seqnr%2==0) {
                    values(freqCell, ts) = real(itsSols[ts](pol/3,st,freqCell));
                  } else {
                    values(freqCell, ts) = imag(itsSols[ts](pol/3,st,freqCell));
                  }
                } else if (itsMode=="scalarphase" || itsMode=="phaseonly") {
                  values(freqCell, ts) = arg(itsSols[ts](pol/3,st,freqCell));
                }
              }
            }
            seqnr++;
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
      if (itsNTimes!=0) {
        stefcal();
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
