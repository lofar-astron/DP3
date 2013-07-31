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
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
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
        itsCorrectType (parset.getString (prefix + "correction")),
        itsTimeSlotsPerParmUpdate (parset.getInt (prefix +
            "timeslotsperparmupdate", 100)),
        itsTimeStep    (0),
        itsNCorr       (0),
        itsTimeInterval (-1),
        itsLastTime    (-1),
        itsUseAP       (false)
    {
      ASSERT (!itsParmDBName.empty());
      // Possible corrections one (or more?) of:
      //   Gain (real/imag or ampl/phase), RM, TEC, Clock, Bandpass
    }

    ApplyCal::~ApplyCal()
    {}

    void ApplyCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();
      itsTimeInterval = infoIn.timeInterval();
      itsNCorr = infoIn.ncorr();

      ASSERT(itsNCorr==4);

      itsParmDB.reset(new BBS::ParmFacade(itsParmDBName));

      // Handle the correction type.

      string corrType = toLower(itsCorrectType);


      if (corrType == "gain") {
        itsUseAP = (itsParmDB->getNames("Gain:*:Real:*").empty() &&
            itsParmDB->getDefNames("Gain:*:Real:*").empty());
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
      } else if (corrType == "tec") {
        itsParmExprs.push_back("TEC");
      } else if (corrType == "clock") {
        if (itsParmDB->getNames("Clock:0:*").empty() &&
            itsParmDB->getDefNames("Clock:0:*").empty() ) {
          itsParmExprs.push_back("Clock");
        }
        else {
          itsParmExprs.push_back("Clock:0");
          itsParmExprs.push_back("Clock:1");
        }
      } else {
        THROW (Exception, "Correction type " + itsCorrectType +
                         " is unknown");
      }

      initDataArrays();
    }

    void ApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << std::endl;
      os << "  parmdb:         " << itsParmDBName << endl;
      os << "  correction:     " << itsCorrectType << endl;
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

      float* weight = buf.getWeights().data();

      size_t nchan = buf.getData().shape()[1];

#pragma omp parallel for
      for (size_t bl=0; bl<nbl; ++bl) {
        for (size_t chan=0;chan<nchan;chan++) {
          if (itsCorrectType=="gain") {
            applyGain( &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                info().getAnt1()[bl], info().getAnt2()[bl], chan, itsTimeStep);
          }
          else if (itsCorrectType=="tec" || itsCorrectType=="clock") {
            applyPhase( &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                info().getAnt1()[bl], info().getAnt2()[bl], chan, itsTimeStep);
          }
        }
      }

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

      itsLastTime = bufStartTime + itsTimeSlotsPerParmUpdate * itsTimeInterval;

      map<string, vector<double> > parmMap;
      map<string, vector<double> >::iterator parmIt;

      uint tfDomainSize=itsTimeSlotsPerParmUpdate*numFreqs;

      for (uint parmNum = 0; parmNum<itsParmExprs.size();++parmNum) {
        // parmMap contains parameter values for all antennas
        parmMap = itsParmDB->getValuesMap( itsParmExprs[parmNum] + "*",
        minFreq, maxFreq, freqInterval,
        bufStartTime, itsLastTime, itsTimeInterval, true);

        for (int ant = 0; ant < numAnts; ++ant) {
          parmIt = parmMap.find(
                    itsParmExprs[parmNum] + ":" + info().antennaNames()[ant]);

          if (parmIt != parmMap.end()) {
            parmvalues[parmNum][ant].swap(parmIt->second);
          } else {// No value found, try default
            Array<double> defValues;
            double defValue;

            if (itsParmDB->getDefValues(itsParmExprs[parmNum] + ":" +
                info().antennaNames()[ant]).size()==1) { // Default for antenna
              itsParmDB->getDefValues(itsParmExprs[parmNum] + ":" +
                  info().antennaNames()[ant]).get(0,defValues);
              ASSERT(defValues.size()==1);
              defValue=defValues.data()[0];
            }
            else if (itsParmDB->getDefValues(itsParmExprs[parmNum]).size()
                == 1) { //Default value
              //TODO: not including * in the pattern above may be too strict
              itsParmDB->getDefValues(itsParmExprs[parmNum]).get(0,defValues);
              ASSERT(defValues.size()==1);
              defValue=defValues.data()[0];
            }
            else {
              THROW (Exception, "No parameter value found for "+
                 itsParmExprs[parmNum]+":"+info().antennaNames()[ant]);
            }

            parmvalues[parmNum][ant].resize(tfDomainSize);
            for (uint tf=0; tf<tfDomainSize;++tf) {
              parmvalues[parmNum][ant][tf]=defValue;
            }
          }
        }
      }

      ASSERT(tfDomainSize==parmvalues[0][0].size());

      double freq;

      // Make parameters complex
      for (uint tf=0;tf<tfDomainSize;++tf) {
        for (int ant=0;ant<numAnts;++ant) {

          freq=info().chanFreqs()[tf % numFreqs];

          if (itsCorrectType=="gain") {
            if (itsUseAP) { // Data as Amplitude / Phase
              itsParms0[ant][tf] = polar(parmvalues[0][ant][tf],
                               parmvalues[1][ant][tf]);
              itsParms1[ant][tf] = polar(parmvalues[2][ant][tf],
                               parmvalues[3][ant][tf]);
            } else { // Data as Real / Imaginary
              itsParms0[ant][tf] = DComplex(parmvalues[0][ant][tf],
                                 parmvalues[1][ant][tf]);
              itsParms1[ant][tf] = DComplex(parmvalues[2][ant][tf],
                                 parmvalues[3][ant][tf]);
            }
          }
          else if (itsCorrectType=="tec") {
            itsParms0[ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * -8.44797245e9 / freq);
            itsParms1[ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * -8.44797245e9 / freq);
          }
          else if (itsCorrectType=="clock") {
            itsParms0[ant][tf]=polar(1.,
                parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            if (itsParmExprs.size() == 1) {
              itsParms1[ant][tf]=polar(1.,
                  parmvalues[0][ant][tf] * freq * casa::C::_2pi);
            }
            else {
              itsParms1[ant][tf]=polar(1.,
                  parmvalues[1][ant][tf] * freq * casa::C::_2pi);
            }
          }
        }
      }
    }

    void ApplyCal::initDataArrays() {
      uint numAnts=info().antennaNames().size();
      uint tfDomainSize=itsTimeSlotsPerParmUpdate*info().chanFreqs().size();

      itsParms0.resize(numAnts);
      itsParms1.resize(numAnts);
      for (uint ant=0;ant<numAnts;++ant) {
          itsParms0[ant].resize(tfDomainSize);
          itsParms1[ant].resize(tfDomainSize);
      }
    }

    void ApplyCal::applyGain (Complex* vis, float* weight, int antA,
        int antB, int chan, int time) {
      int timeFreqOffset=(time*info().nchan())+chan;

      DComplex gain00A = itsParms0[antA][timeFreqOffset];
      DComplex gain11A = itsParms1[antA][timeFreqOffset];
      DComplex gain00B = itsParms0[antB][timeFreqOffset];
      DComplex gain11B = itsParms1[antB][timeFreqOffset];

      vis[0] /= gain00A * conj(gain00B);
      vis[1] /= gain00A * conj(gain11B);
      vis[2] /= gain11A * conj(gain00B);
      vis[3] /= gain11A * conj(gain11B);

      //cout<<weight[0]<<endl;
      //weight[0]*= real(gain00A) * real(gain00A) * real(gain00B) * real(gain00B);
      //weight[1]*= real(gain00A) * real(gain00A) * real(gain11B) * real(gain11B);
      //weight[2]*= real(gain11A) * real(gain11A) * real(gain00B) * real(gain00B);
      //weight[3]*= real(gain11A) * real(gain11A) * real(gain11B) * real(gain11B);
    }

    void ApplyCal::applyPhase(Complex* vis, int antA, int antB,
        int chan, int time) {
      int timeFreqOffset=(time*info().nchan())+chan;
      for (size_t i=0;i<itsNCorr;i++) {
        vis[i] /= itsParms0[antA][timeFreqOffset]/
            itsParms0[antB][timeFreqOffset];
      }
    }

  } //# end namespace
}
