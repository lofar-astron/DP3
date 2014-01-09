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
                        const string& prefix)
      : itsInput       (input),
        itsName        (prefix),
        itsParmDBName  (parset.getString (prefix + "parmdb")),
        itsCorrectType (toLower(parset.getString (prefix + "correction"))),
        itsTimeSlotsPerParmUpdate (parset.getInt (prefix +
            "timeslotsperparmupdate", 500)),
        itsSigmaMMSE   (parset.getDouble (prefix + "MMSE.Sigma", 0)),
        itsUpdateWeights (parset.getBool (prefix + "updateweights", false)),
        itsTimeStep    (0),
        itsNCorr       (0),
        itsTimeInterval (-1),
        itsLastTime    (-1),
        itsUseAP       (false)
    {
      ASSERT (!itsParmDBName.empty());
    }

    ApplyCal::~ApplyCal()
    {}

    void ApplyCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();
      if (itsUpdateWeights) {
        info().setNeedWrite(info().needWrite() | DPInfo::NeedWriteWeight);
      }
      itsTimeInterval = infoIn.timeInterval();
      itsNCorr = infoIn.ncorr();

      ASSERT(itsNCorr==4);

      itsParmDB.reset(new BBS::ParmFacade(itsParmDBName));

      // Handle the correction type.

      if ((itsCorrectType == "gain" || itsCorrectType=="fullgain") &&
          (itsParmDB->getNames("Gain:0:1:*").size() +
           itsParmDB->getDefNames("Gain:0:1:*").size() >0 )) {
        itsCorrectType="fullgain";
      }

      if (itsCorrectType == "gain") {
        itsUseAP = (itsParmDB->getNames("Gain:0:0:Real*").empty() &&
            itsParmDB->getDefNames("Gain:0:0:Real*").empty());
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
      } else if (itsCorrectType == "fullgain") {
        itsUseAP = (itsParmDB->getNames("Gain:0:0:Real*").empty() &&
            itsParmDB->getDefNames("Gain:0:0:Real*").empty());
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
        itsParmExprs.push_back("TEC");
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
        itsParmExprs.push_back("CommonRotationAngle");
      } else if (itsCorrectType == "commonscalarphase") {
        itsParmExprs.push_back("CommonScalarPhase");
      }
      else {
        THROW (Exception, "Correction type " + itsCorrectType +
                         " is unknown");
      }

      initDataArrays();
      itsFlagCounter.init(getInfo());
    }

    void ApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << std::endl;
      os << "  parmdb:         " << itsParmDBName << endl;
      os << "  correction:     " << itsCorrectType << endl;
      os << "  sigmaMMSE:      " << itsSigmaMMSE << endl;
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
      DPBuffer buf(bufin);
      buf.getData().unique();
      RefRows rowNrs(buf.getRowNrs());

      double bufStartTime = buf.getTime() - 0.5*itsTimeInterval;

      if (buf.getTime() > itsLastTime) {
        updateParms(bufStartTime);
        itsTimeStep=0;
      }
      else {
        itsTimeStep++;
      }

      // Loop through all baselines in the buffer.
      size_t nbl = bufin.getData().shape()[2];

      Complex* data = buf.getData().data();

      if (itsUpdateWeights) {
        buf.setWeights(itsInput->fetchWeights (buf, rowNrs, itsTimer));
      }
      float* weight = buf.getWeights().data();

      size_t nchan = buf.getData().shape()[1];

#pragma omp parallel for
      for (size_t bl=0; bl<nbl; ++bl) {
        for (size_t chan=0;chan<nchan;chan++) {
          if (itsParms.size()>2) {
            applyFull( &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                info().getAnt1()[bl], info().getAnt2()[bl], chan, itsTimeStep);
          }
          else {
            applyDiag( &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                info().getAnt1()[bl], info().getAnt2()[bl], chan, itsTimeStep);
          }
        }
      }

      MSReader::flagInfNaN(buf.getData(),buf.getFlags(),itsFlagCounter);

      itsTimer.stop();
      getNextStep()->process(buf);
      return false;
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
      // (e.g. Gain:0:0 and Gain:1:1), next all antennas, next over frec * time
      // as returned by ParmDB
      vector<vector<vector<double> > > parmvalues;
      parmvalues.resize(itsParmExprs.size());
      for (size_t i=0;i<parmvalues.size();++i) {
        parmvalues[i].resize(numAnts);
      }

      uint numFreqs         (info().chanFreqs().size());
      double freqInterval  (info().chanWidths()[0]);
      double minFreq       (info().chanFreqs()[0]-0.5*freqInterval);
      double maxFreq (info().chanFreqs()[numFreqs-1]+0.5*freqInterval);

      itsLastTime = std::min(
          bufStartTime + itsTimeSlotsPerParmUpdate * itsTimeInterval,
          info().startTime() + info().ntime() * itsTimeInterval);

      map<string, vector<double> > parmMap;
      map<string, vector<double> >::iterator parmIt;

      uint tfDomainSize=itsTimeSlotsPerParmUpdate*numFreqs;

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
              itsParms[0][ant][tf] = polar(parmvalues[0][ant][tf],
                               parmvalues[1][ant][tf]);
              itsParms[1][ant][tf] = polar(parmvalues[2][ant][tf],
                               parmvalues[3][ant][tf]);
            } else { // Data as Real / Imaginary
              itsParms[0][ant][tf] = DComplex(parmvalues[0][ant][tf],
                                 parmvalues[1][ant][tf]);
              itsParms[1][ant][tf] = DComplex(parmvalues[2][ant][tf],
                                 parmvalues[3][ant][tf]);
            }
          }
          else if (itsCorrectType=="fullgain") {
            if (itsUseAP) { // Data as Amplitude / Phase
              itsParms[0][ant][tf] = polar(parmvalues[0][ant][tf],
                               parmvalues[1][ant][tf]);
              itsParms[1][ant][tf] = polar(parmvalues[2][ant][tf],
                               parmvalues[3][ant][tf]);
              itsParms[2][ant][tf] = polar(parmvalues[4][ant][tf],
                               parmvalues[5][ant][tf]);
              itsParms[3][ant][tf] = polar(parmvalues[6][ant][tf],
                               parmvalues[7][ant][tf]);
            } else { // Data as Real / Imaginary
              itsParms[0][ant][tf] = DComplex(parmvalues[0][ant][tf],
                                 parmvalues[1][ant][tf]);
              itsParms[1][ant][tf] = DComplex(parmvalues[2][ant][tf],
                                 parmvalues[3][ant][tf]);
              itsParms[2][ant][tf] = DComplex(parmvalues[4][ant][tf],
                                 parmvalues[5][ant][tf]);
              itsParms[3][ant][tf] = DComplex(parmvalues[6][ant][tf],
                                 parmvalues[7][ant][tf]);
            }
          }
          else if (itsCorrectType=="tec") {
            itsParms[0][ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * -8.44797245e9 / freq);
            itsParms[1][ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * -8.44797245e9 / freq);
          }
          else if (itsCorrectType=="clock") {
            itsParms[0][ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            if (itsParmExprs.size() == 1) { // No Clock:0, only Clock:
              itsParms[1][ant][tf]=polar(1.,
                  parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            }
            else { // Clock:0 and Clock:1
              itsParms[1][ant][tf]=polar(1.,
                  parmvalues[1][ant][tf] * freq * casa::C::_2pi);
            }
          }
          else if (itsCorrectType=="commonrotationangle") {
            itsParms[0][ant][tf] =  cos(parmvalues[0][ant][tf]);
            itsParms[1][ant][tf] = -sin(parmvalues[0][ant][tf]);
            itsParms[2][ant][tf] =  sin(parmvalues[0][ant][tf]);
            itsParms[3][ant][tf] =  cos(parmvalues[0][ant][tf]);
          }
          else if (itsCorrectType=="commonscalarphase") {
            itsParms[0][ant][tf] = polar(1., parmvalues[0][ant][tf]);
            itsParms[1][ant][tf] = polar(1., parmvalues[0][ant][tf]);
          }
        }
      }
    }

    void ApplyCal::initDataArrays() {
      uint numAnts=info().antennaNames().size();
      uint tfDomainSize=itsTimeSlotsPerParmUpdate*info().chanFreqs().size();

      uint numParms;
      if (itsCorrectType=="fullgain" || itsCorrectType=="commonrotationangle") {
        numParms = 4;
      }
      else {
        numParms = 2;
      }

      itsParms.resize(numParms);

      for (uint parmNum=0;parmNum<numParms;parmNum++) {
        itsParms[parmNum].resize(numAnts);
        for (uint ant=0;ant<numAnts;++ant) {
            itsParms[parmNum][ant].resize(tfDomainSize);
        }
      }
    }

    void ApplyCal::applyDiag (Complex* vis, float* weight, int antA,
        int antB, int chan, int time) {
      int timeFreqOffset=(time*info().nchan())+chan;

      DComplex diag0A = itsParms[0][antA][timeFreqOffset];
      DComplex diag1A = itsParms[1][antA][timeFreqOffset];
      DComplex diag0B = itsParms[0][antB][timeFreqOffset];
      DComplex diag1B = itsParms[1][antB][timeFreqOffset];

      vis[0] /= diag0A * conj(diag0B);
      vis[1] /= diag0A * conj(diag1B);
      vis[2] /= diag1A * conj(diag0B);
      vis[3] /= diag1A * conj(diag1B);

      if (itsUpdateWeights) {
        weight[0] *= norm(diag0A) * norm(diag0B);
        weight[1] *= norm(diag0A) * norm(diag1B);
        weight[2] *= norm(diag1A) * norm(diag0B);
        weight[3] *= norm(diag1A) * norm(diag1B);
      }
    }

    // Inverts complex 2x2 input matrix
    void ApplyCal::invert (DComplex* v, double sigmaMMSE) const
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

    void ApplyCal::applyFull (Complex* vis, float* weight, int antA,
        int antB, int chan, int time) {
      int timeFreqOffset=(time*info().nchan())+chan;
      DComplex gainA[4];
      DComplex gainB[4];

      gainA[0] = itsParms[0][antA][timeFreqOffset];
      gainA[1] = itsParms[1][antA][timeFreqOffset];
      gainA[2] = itsParms[2][antA][timeFreqOffset];
      gainA[3] = itsParms[3][antA][timeFreqOffset];

      gainB[0] = itsParms[0][antB][timeFreqOffset];
      gainB[1] = itsParms[1][antB][timeFreqOffset];
      gainB[2] = itsParms[2][antB][timeFreqOffset];
      gainB[3] = itsParms[3][antB][timeFreqOffset];

      DComplex gainAxvis[4];
      invert(gainA,itsSigmaMMSE);
      invert(gainB,itsSigmaMMSE);

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

      // The code below does the same as the combination of BBS + python script
      // covariance2weight.py (cookbook), except it stores weights per freq.
      // The diagonal of covariance matrix is transferred to the weights.
      // Note that the real covariance (mixing of noise terms after which they
      // are not independent anymore) is not stored.
      // The input covariance matrix C is assumed to be diagonal with elements
      // w_i (the weights), the result the diagonal of
      // (gainA kronecker gainB^H).C.(gainA kronecker gainB^H)^H
      if (itsUpdateWeights) {
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
    }

  } //# end namespace
}
