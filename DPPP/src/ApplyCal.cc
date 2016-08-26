//# ApplyCal.cc: DPPP step class to apply a calibration correction to the data
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
//# $Id: ApplyCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/MSReader.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/OS/File.h>
#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

/// Look at BBSKernel MeasurementExprLOFARUtil.cc and Apply.cc

namespace LOFAR {
  namespace DPPP {

    ApplyCal::ApplyCal (DPInput* input,
                        const ParameterSet& parset,
                        const string& prefix,
                        bool substep)
      : itsInput       (input),
        itsName        (prefix),
        itsParmDBName  (parset.getString (prefix + "parmdb")),
        itsCorrectType (toLower(parset.getString (prefix + "correction", "gain"))),
        itsTimeSlotsPerParmUpdate (parset.getInt (prefix +
            "timeslotsperparmupdate", 500)),
        itsSigmaMMSE   (parset.getDouble (prefix + "MMSE.Sigma", 0)),
        itsUpdateWeights (parset.getBool (prefix + "updateweights", false)),
        itsCount       (0),
        itsTimeStep    (0),
        itsNCorr       (0),
        itsTimeInterval (-1),
        itsLastTime    (-1),
        itsUseAP       (false)
    {
      ASSERT (!itsParmDBName.empty());
      if (substep) {
        itsInvert=false;
      } else {
        itsInvert=parset.getBool (prefix + "invert", true);
      }
      if (itsCorrectType=="fulljones" && itsUpdateWeights) {
        ASSERTSTR (itsInvert, "Updating weights has not been implemented for invert=false and fulljones");
      }
    }

    ApplyCal::ApplyCal()
    {}

    ApplyCal::~ApplyCal()
    {}

    void ApplyCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
      info().setWriteFlags();
      if (itsUpdateWeights) {
        info().setWriteWeights();
      }
      itsTimeInterval = infoIn.timeInterval();
      itsNCorr = infoIn.ncorr();

      ASSERT(itsNCorr==4);

      itsParmDB.reset(new BBS::ParmFacade(itsParmDBName));

      // Detect if full jones solutions are present
      if ((itsCorrectType == "gain" || itsCorrectType=="fulljones") &&
          (itsParmDB->getNames("Gain:0:1:*").size() +
           itsParmDB->getDefNames("Gain:0:1:*").size() >0 )) {
        itsCorrectType="fulljones";
      }

      // Detect if solutions are saved as Real/Imag or Ampl/Phase 
      if (itsCorrectType == "gain" || itsCorrectType == "fulljones" ){
        if (!itsParmDB->getNames("Gain:0:0:Real*").empty()) {
          // Values with :Real present
          itsUseAP = false; 
        } else if (!itsParmDB->getNames("Gain:0:0:Ampl*").empty() || 
                   !itsParmDB->getNames("Phase:0:0:Ampl*").empty()) {
          // Values with :Ampl present
          itsUseAP = true;
        } else if (!itsParmDB->getDefNames("Gain:0:0:Real*").empty()) {
          // Defvalues with :Real present
          itsUseAP = false;
        } else if (!itsParmDB->getDefNames("Gain:0:0:Ampl*").empty() ||
                   !itsParmDB->getDefNames("Gain:0:0:Phase*").empty()) {
          // Defvalues with :Ampl present
          itsUseAP = true;
        } else {
          THROW (Exception, "No gains found in parmdb "+itsParmDBName);
        }
      }

      if (itsCorrectType == "gain") {
        if (itsUseAP) {
          itsParmExprs.push_back("Gain:0:0:Ampl");
          itsParmExprs.push_back("Gain:0:0:Phase");
          itsParmExprs.push_back("Gain:1:1:Ampl");
          itsParmExprs.push_back("Gain:1:1:Phase");
        } else {
          itsParmExprs.push_back("Gain:0:0:Real");
          itsParmExprs.push_back("Gain:0:0:Imag");
          itsParmExprs.push_back("Gain:1:1:Real");
          itsParmExprs.push_back("Gain:1:1:Imag");
        }
      } else if (itsCorrectType == "fulljones") {
        if (itsUseAP) {
          itsParmExprs.push_back("Gain:0:0:Ampl");
          itsParmExprs.push_back("Gain:0:0:Phase");
          itsParmExprs.push_back("Gain:0:1:Ampl");
          itsParmExprs.push_back("Gain:0:1:Phase");
          itsParmExprs.push_back("Gain:1:0:Ampl");
          itsParmExprs.push_back("Gain:1:0:Phase");
          itsParmExprs.push_back("Gain:1:1:Ampl");
          itsParmExprs.push_back("Gain:1:1:Phase");
        } else {
          itsParmExprs.push_back("Gain:0:0:Real");
          itsParmExprs.push_back("Gain:0:0:Imag");
          itsParmExprs.push_back("Gain:0:1:Real");
          itsParmExprs.push_back("Gain:0:1:Imag");
          itsParmExprs.push_back("Gain:1:0:Real");
          itsParmExprs.push_back("Gain:1:0:Imag");
          itsParmExprs.push_back("Gain:1:1:Real");
          itsParmExprs.push_back("Gain:1:1:Imag");
        }
      } else if (itsCorrectType == "tec") {
        if (itsParmDB->getNames("TEC:0:*").empty() &&
                    itsParmDB->getDefNames("TEC:0:*").empty() ) {
          itsParmExprs.push_back("TEC");
        }
        else {
          itsParmExprs.push_back("TEC:0");
          itsParmExprs.push_back("TEC:1");
        }
      } else if (itsCorrectType == "clock") {
        if (itsParmDB->getNames("Clock:0:*").empty() &&
            itsParmDB->getDefNames("Clock:0:*").empty() ) {
          itsParmExprs.push_back("Clock");
        }
        else {
          itsParmExprs.push_back("Clock:0");
          itsParmExprs.push_back("Clock:1");
        }
      } else if (itsCorrectType == "commonrotationangle") {
        itsParmExprs.push_back("{Common,}RotationAngle");
      } else if (itsCorrectType == "commonscalarphase") {
        itsParmExprs.push_back("{Common,}ScalarPhase");
      } else if (itsCorrectType == "rotationmeasure") {
        itsParmExprs.push_back("RotationMeasure");
      } else if (itsCorrectType == "commonscalaramplitude") {
        itsParmExprs.push_back("{Common,}ScalarAmplitude");
      }
      else {
        THROW (Exception, "Correction type " + itsCorrectType +
                         " is unknown");
      }

      initDataArrays();
      itsFlagCounter.init(getInfo());

      // Check that channels are evenly spaced
      if (info().nchan()>1) {
        Vector<Double> upFreq = info().chanFreqs()(
                                  Slicer(IPosition(1,1),
                                         IPosition(1,info().nchan()-1)));
        Vector<Double> lowFreq = info().chanFreqs()(
                                  Slicer(IPosition(1,0),
                                         IPosition(1,info().nchan()-1)));
        Double freqstep0=upFreq(0)-lowFreq(0);
        // Compare up to 1kHz accuracy
        bool regularChannels=allNearAbs(upFreq-lowFreq, freqstep0, 1.e3) &&
                             allNearAbs(info().chanWidths(),
                                        info().chanWidths()(0), 1.e3);
        ASSERTSTR(regularChannels, 
                  "ApplyCal requires evenly spaced channels.");
      }
    }

    void ApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << std::endl;
      os << "  parmdb:         " << itsParmDBName << endl;
      os << "  correction:     " << itsCorrectType << endl;
      if (itsCorrectType=="gain" || itsCorrectType=="fulljones") {
        os << "    Ampl/Phase:   " << boolalpha << itsUseAP << endl;
      }
      os << "  update weights: " << boolalpha << itsUpdateWeights << endl;
      os << "  sigmaMMSE:      " << itsSigmaMMSE << endl;
      os << "  invert:         " << boolalpha << itsInvert <<endl;
      os << "  timeSlotsPerParmUpdate: " << itsTimeSlotsPerParmUpdate <<endl;
    }

    void ApplyCal::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " ApplyCal " << itsName << endl;
    }

    bool ApplyCal::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      double bufStartTime = bufin.getTime() - 0.5*itsTimeInterval;

      if (bufin.getTime() > itsLastTime) {
        updateParms(bufStartTime);
        itsTimeStep=0;
      }
      else {
        itsTimeStep++;
      }

      // Loop through all baselines in the buffer.
      size_t nbl = itsBuffer.getData().shape()[2];

      Complex* data = itsBuffer.getData().data();

      itsInput->fetchWeights (bufin, itsBuffer, itsTimer);
      float* weight = itsBuffer.getWeights().data();

      bool* flag = itsBuffer.getFlags().data();

      size_t nchan = itsBuffer.getData().shape()[1];

#pragma omp parallel for
      for (size_t bl=0; bl<nbl; ++bl) {
        for (size_t chan=0;chan<nchan;chan++) {
          uint timeFreqOffset=(itsTimeStep*info().nchan())+chan;
          uint antA = info().getAnt1()[bl];
          uint antB = info().getAnt2()[bl];
          if (itsParms.shape()[0]>2) {
            applyFull( &itsParms(0, antA, timeFreqOffset),
                       &itsParms(0, antB, timeFreqOffset),
                       &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &flag[  bl * itsNCorr * nchan + chan * itsNCorr ],
                       bl, chan, itsUpdateWeights, itsFlagCounter);
          }
          else {
            applyDiag( &itsParms(0, antA, timeFreqOffset),
                       &itsParms(0, antB, timeFreqOffset),
                       &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &flag[  bl * itsNCorr * nchan + chan * itsNCorr ],
                       bl, chan, itsUpdateWeights, itsFlagCounter);
          }
        }
      }

      itsTimer.stop();
      getNextStep()->process(itsBuffer);

      itsCount++;
      return true;
    }

    void ApplyCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


    void ApplyCal::updateParms (const double bufStartTime)
    {
      int numAnts = info().antennaNames().size();

      // itsParms contains the parameters to a grid, first for all parameters
      // (e.g. Gain:0:0 and Gain:1:1), next all antennas, next over freq * time
      // as returned by ParmDB
      vector<vector<vector<double> > > parmvalues;
      parmvalues.resize(itsParmExprs.size());
      for (size_t i=0;i<parmvalues.size();++i) {
        parmvalues[i].resize(numAnts);
      }

      uint numFreqs         (info().chanFreqs().size());
      double freqInterval  (info().chanWidths()[0]);
      if (numFreqs>1) { // Handle data with evenly spaced gaps between channels
        freqInterval = info().chanFreqs()[1]-info().chanFreqs()[0];
      }
      double minFreq       (info().chanFreqs()[0]-0.5*freqInterval);
      double maxFreq (info().chanFreqs()[numFreqs-1]+0.5*freqInterval);
      itsLastTime = bufStartTime + itsTimeSlotsPerParmUpdate * itsTimeInterval;
      uint numTimes = itsTimeSlotsPerParmUpdate;

      double lastMSTime = info().startTime() + info().ntime() * itsTimeInterval;
      if (itsLastTime > lastMSTime && !nearAbs(itsLastTime, lastMSTime, 1.e-3)) {
        itsLastTime = lastMSTime;
        numTimes = info().ntime() % itsTimeSlotsPerParmUpdate;
      }

      map<string, vector<double> > parmMap;
      map<string, vector<double> >::iterator parmIt;

      uint tfDomainSize=numTimes*numFreqs;

      for (uint parmExprNum = 0; parmExprNum<itsParmExprs.size();++parmExprNum) {
        // parmMap contains parameter values for all antennas
        parmMap = itsParmDB->getValuesMap( itsParmExprs[parmExprNum] + "*",
        minFreq, maxFreq, freqInterval,
        bufStartTime, itsLastTime, itsTimeInterval, true);

        for (int ant = 0; ant < numAnts; ++ant) {
          parmIt = parmMap.find(
                    itsParmExprs[parmExprNum] + ":" + info().antennaNames()[ant]);

          if (parmIt != parmMap.end()) {
            parmvalues[parmExprNum][ant].swap(parmIt->second);
            ASSERT(parmvalues[parmExprNum][ant].size()==tfDomainSize);
          } else {// No value found, try default
            Array<double> defValues;
            double defValue;

            if (itsParmDB->getDefValues(itsParmExprs[parmExprNum] + ":" +
                info().antennaNames()[ant]).size()==1) { // Default for antenna
              itsParmDB->getDefValues(itsParmExprs[parmExprNum] + ":" +
                  info().antennaNames()[ant]).get(0,defValues);
              ASSERT(defValues.size()==1);
              defValue=defValues.data()[0];
            }
            else if (itsParmDB->getDefValues(itsParmExprs[parmExprNum]).size()
                == 1) { //Default value
              itsParmDB->getDefValues(itsParmExprs[parmExprNum]).get(0,defValues);
              ASSERT(defValues.size()==1);
              defValue=defValues.data()[0];
            } else if (itsParmExprs[parmExprNum].substr(0,5)=="Gain:") {
              defValue=0.;
            }
            else {
              THROW (Exception, "No parameter value found for "+
                 itsParmExprs[parmExprNum]+":"+info().antennaNames()[ant]);
            }

            parmvalues[parmExprNum][ant].resize(tfDomainSize);
            for (uint tf=0; tf<tfDomainSize;++tf) {
              parmvalues[parmExprNum][ant][tf]=defValue;
            }
          }
        }
      }

      ASSERT(parmvalues[0][0].size() <= tfDomainSize); // Catches multiple matches

      double freq;

      // Make parameters complex
      for (uint tf=0;tf<tfDomainSize;++tf) {
        for (int ant=0;ant<numAnts;++ant) {

          freq=info().chanFreqs()[tf % numFreqs];

          if (itsCorrectType=="gain") {
            if (itsUseAP) { // Data as Amplitude / Phase
              itsParms(0, ant, tf) = polar(parmvalues[0][ant][tf],
                                           parmvalues[1][ant][tf]);
              itsParms(1, ant, tf) = polar(parmvalues[2][ant][tf],
                                           parmvalues[3][ant][tf]);
            } else { // Data as Real / Imaginary
              itsParms(0, ant, tf) = DComplex(parmvalues[0][ant][tf],
                                              parmvalues[1][ant][tf]);
              itsParms(1, ant, tf) = DComplex(parmvalues[2][ant][tf],
                                              parmvalues[3][ant][tf]);
            }
          }
          else if (itsCorrectType=="fulljones") {
            if (itsUseAP) { // Data as Amplitude / Phase
              itsParms(0, ant, tf) = polar(parmvalues[0][ant][tf],
                                           parmvalues[1][ant][tf]);
              itsParms(1, ant, tf) = polar(parmvalues[2][ant][tf],
                                           parmvalues[3][ant][tf]);
              itsParms(2, ant, tf) = polar(parmvalues[4][ant][tf],
                                           parmvalues[5][ant][tf]);
              itsParms(3, ant, tf) = polar(parmvalues[6][ant][tf],
                                           parmvalues[7][ant][tf]);
            } else { // Data as Real / Imaginary
              itsParms(0, ant, tf) = DComplex(parmvalues[0][ant][tf],
                                              parmvalues[1][ant][tf]);
              itsParms(1, ant, tf) = DComplex(parmvalues[2][ant][tf],
                                              parmvalues[3][ant][tf]);
              itsParms(2, ant, tf) = DComplex(parmvalues[4][ant][tf],
                                              parmvalues[5][ant][tf]);
              itsParms(3, ant, tf) = DComplex(parmvalues[6][ant][tf],
                                              parmvalues[7][ant][tf]);
            }
          }
          else if (itsCorrectType=="tec") {
            itsParms(0, ant, tf)=polar(1.,
                parmvalues[0][ant][tf] * -8.44797245e9 / freq);
            if (itsParmExprs.size() == 1) { // No TEC:0, only TEC:
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[0][ant][tf] * -8.44797245e9 / freq);
            }
            else { // TEC:0 and TEC:1
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[1][ant][tf] * -8.44797245e9 / freq);
            }
          }
          else if (itsCorrectType=="clock") {
            itsParms(0, ant, tf)=polar(1.,
                parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            if (itsParmExprs.size() == 1) { // No Clock:0, only Clock:
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            }
            else { // Clock:0 and Clock:1
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[1][ant][tf] * freq * casa::C::_2pi);
            }
          }
          else if (itsCorrectType=="commonrotationangle") {
            double phi=parmvalues[0][ant][tf];
            if (itsInvert) {
              phi = -phi;
            }
            double sinv=sin(phi);
            double cosv=cos(phi);
            itsParms(0, ant, tf) =  cosv;
            itsParms(1, ant, tf) = -sinv;
            itsParms(2, ant, tf) =  sinv;
            itsParms(3, ant, tf) =  cosv;
          }
          else if (itsCorrectType=="rotationmeasure") {
            double lambda2 = casa::C::c / freq;
            lambda2 *= lambda2;
            double chi = parmvalues[0][ant][tf] * lambda2;
            if (itsInvert) {
              chi = -chi;
            }
            double sinv = sin(chi);
            double cosv = cos(chi);
            itsParms(0, ant, tf) =  cosv;
            itsParms(1, ant, tf) = -sinv;
            itsParms(2, ant, tf) =  sinv;
            itsParms(3, ant, tf) =  cosv;
          }
          else if (itsCorrectType=="commonscalarphase") {
            itsParms(0, ant, tf) = polar(1., parmvalues[0][ant][tf]);
            itsParms(1, ant, tf) = polar(1., parmvalues[0][ant][tf]);
          }
          else if (itsCorrectType=="commonscalaramplitude") {
            itsParms(0, ant, tf) = parmvalues[0][ant][tf];
            itsParms(1, ant, tf) = parmvalues[0][ant][tf];
          }

          // Invert (rotationmeasure and commonrotationangle are already inverted)
          if (itsInvert && itsParms.shape()[0]==2) {
            itsParms(0, ant, tf) = 1./itsParms(0, ant, tf);
            itsParms(1, ant, tf) = 1./itsParms(1, ant, tf);
          }
          else if (itsInvert && itsCorrectType=="fulljones") {
            invert(&itsParms(0, ant, tf),itsSigmaMMSE);
          }
        }
      }
    }

    void ApplyCal::initDataArrays() {
      uint numAnts=info().antennaNames().size();
      uint tfDomainSize=itsTimeSlotsPerParmUpdate*info().chanFreqs().size();

      uint numParms;
      if (itsCorrectType=="fulljones" || 
          itsCorrectType=="commonrotationangle" || 
          itsCorrectType=="rotationmeasure") {
        numParms = 4;
      }
      else {
        numParms = 2;
      }

      itsParms.resize(numParms, numAnts, tfDomainSize);
    }

    void ApplyCal::applyDiag (const DComplex* gainA, const DComplex* gainB,
                              Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool updateWeights,
                              FlagCounter& flagCounter) {
      // If parameter is NaN or inf, do not apply anything and flag the data
      if (! (isFinite(gainA[0].real()) && isFinite(gainA[0].imag()) &&
             isFinite(gainB[0].real()) && isFinite(gainB[0].imag()) &&
             isFinite(gainA[1].real()) && isFinite(gainA[1].imag()) &&
             isFinite(gainB[1].real()) && isFinite(gainB[1].imag())) ) {
        // Only update flagcounter for first correlation
        if (!flag[0]) {
          flagCounter.incrChannel(chan);
          flagCounter.incrBaseline(bl);
        }
        for (uint corr=0; corr<4; ++corr) {
          flag[corr]=true;
        }
        return;
      }

      vis[0] *= gainA[0] * conj(gainB[0]);
      vis[1] *= gainA[0] * conj(gainB[1]);
      vis[2] *= gainA[1] * conj(gainB[0]);
      vis[3] *= gainA[1] * conj(gainB[1]);

      if (updateWeights) {
        weight[0] /= norm(gainA[0]) * norm(gainB[0]);
        weight[1] /= norm(gainA[0]) * norm(gainB[1]);
        weight[2] /= norm(gainA[1]) * norm(gainB[0]);
        weight[3] /= norm(gainA[1]) * norm(gainB[1]);
      }
    }

    // Inverts complex 2x2 input matrix
    void ApplyCal::invert (DComplex* v, double sigmaMMSE)
    {
      // Add the variance of the nuisance term to the elements on the diagonal.
      const double variance = sigmaMMSE * sigmaMMSE;
      DComplex v0 = v[0] + variance;
      DComplex v3 = v[3] + variance;
      // Compute inverse in the usual way.
      DComplex invDet(1.0 / (v0 * v3 - v[1] * v[2]));
      v[0] = v3 * invDet;
      v[2] = v[2] * -invDet;
      v[1] = v[1] * -invDet;
      v[3] = v0 * invDet;
    }

    void ApplyCal::applyFull (const DComplex* gainA, const DComplex* gainB,
                              Complex* vis, float* weight, bool* flag,
                              uint bl, uint chan, bool doUpdateWeights,
                              FlagCounter& flagCounter) {
      DComplex gainAxvis[4];

      // If parameter is NaN or inf, do not apply anything and flag the data
      bool anyinfnan = false;
      for (uint corr=0; corr<4; ++corr) {
        if (! (isFinite(gainA[corr].real()) && isFinite(gainA[corr].imag()) &&
               isFinite(gainB[corr].real()) && isFinite(gainB[corr].imag())) ) {
          anyinfnan = true;
          break;
        }
      }
      if (anyinfnan) {
        // Only update flag counter for first correlation
        if (!flag[0]) {
          flagCounter.incrChannel(chan);
          flagCounter.incrBaseline(bl);
        }
        for (uint corr=0; corr<4; ++corr) {
          flag[corr]=true;
        }
        return;
      }

      // gainAxvis = gainA * vis
      for (uint row=0;row<2;++row) {
        for (uint col=0;col<2;++col) {
          gainAxvis[2*row+col]=gainA[2*row+0] * DComplex(vis[2*0+col]) +
                               gainA[2*row+1] * DComplex(vis[2*1+col]);
        }
      }

      // vis = gainAxvis * gainB^H
      for (uint row=0;row<2;++row) {
        for (uint col=0;col<2;++col) {
          vis[2*row+col]=gainAxvis[2*row+0] * conj(gainB[2*col+0])+
                         gainAxvis[2*row+1] * conj(gainB[2*col+1]);
        }
      }

      if (doUpdateWeights) {
        applyWeights(gainA, gainB, weight);
      }
    }

    void ApplyCal::applyWeights(const DComplex* gainA,
                                const DComplex* gainB,
                                float* weight) {
      float cov[4], normGainA[4], normGainB[4];
      for (uint i=0;i<4;++i) {
        cov[i]=1./weight[i];
        normGainA[i]=norm(gainA[i]);
        normGainB[i]=norm(gainB[i]);
      }

      weight[0]=cov[0]*(normGainA[0]*normGainB[0])
                     +cov[1]*(normGainA[0]*normGainB[1])
                     +cov[2]*(normGainA[1]*normGainB[0])
                     +cov[3]*(normGainA[1]*normGainB[1]);
      weight[0]=1./weight[0];

      weight[1]=cov[0]*(normGainA[0]*normGainB[2])
                     +cov[1]*(normGainA[0]*normGainB[3])
                     +cov[2]*(normGainA[1]*normGainB[2])
                     +cov[3]*(normGainA[1]*normGainB[3]);
      weight[1]=1./weight[1];

      weight[2]=cov[0]*(normGainA[2]*normGainB[0])
                     +cov[1]*(normGainA[2]*normGainB[1])
                     +cov[2]*(normGainA[3]*normGainB[0])
                     +cov[3]*(normGainA[3]*normGainB[1]);
      weight[2]=1./weight[2];

      weight[3]=cov[0]*(normGainA[2]*normGainB[2])
                     +cov[1]*(normGainA[2]*normGainB[3])
                     +cov[2]*(normGainA[3]*normGainB[2])
                     +cov[3]*(normGainA[3]*normGainB[3]);
      weight[3]=1./weight[3];
    }

    void ApplyCal::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by ApplyCal " << itsName;
      os << endl << "=======================" << endl;
      itsFlagCounter.showBaseline (os, itsCount);
      itsFlagCounter.showChannel  (os, itsCount);
    }

  } //# end namespace
}
