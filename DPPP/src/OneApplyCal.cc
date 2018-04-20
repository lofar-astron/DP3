//# OneApplyCal.cc: DPPP step class to apply a calibration correction to the data
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
//# $Id: OneApplyCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/OneApplyCal.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/MSReader.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/OS/File.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iomanip>

using namespace casacore;
using namespace LOFAR::BBS;

/// Look at BBSKernel MeasurementExprLOFARUtil.cc and Apply.cc

namespace LOFAR {
  namespace DPPP {

    OneApplyCal::OneApplyCal (DPInput* input,
                        const ParameterSet& parset,
                        const string& prefix,
                        const string& defaultPrefix,
                        bool substep,
                        string predictDirection
                        )
      : itsInput       (input),
        itsName        (prefix),
        itsParmDBName  (
            parset.isDefined(prefix+"parmdb") ?
            parset.getString(prefix + "parmdb") :
            parset.getString(defaultPrefix + "parmdb")),
        itsUseH5Parm   (itsParmDBName.find(".h5") != string::npos),
        itsSigmaMMSE   (
            parset.isDefined(prefix + "MMSE.Sigma") ?
            parset.getDouble(prefix + "MMSE.Sigma") :
            parset.getDouble(defaultPrefix + "MMSE.Sigma", 0.)),
        itsUpdateWeights (
            parset.isDefined(prefix + "updateweights") ?
            parset.getBool (prefix + "updateweights") :
            parset.getBool (defaultPrefix + "updateweights", false)),
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
        itsInvert = (parset.isDefined(prefix + "invert") ?
                     parset.getBool (prefix + "invert") :
                     parset.getBool (defaultPrefix + "invert", true));
      }

      if (itsUseH5Parm) {
        itsTimeSlotsPerParmUpdate = 0;
        string directionStr;
        directionStr = (parset.isDefined(prefix + "direction") ?
                        parset.getString(prefix + "direction") :
                        parset.getString(defaultPrefix + "direction",
                          predictDirection));
        itsH5Parm = H5Parm(itsParmDBName);

        itsSolTabName = (parset.isDefined(prefix + "correction") ?
                         parset.getString(prefix + "correction") :
                         parset.getString(defaultPrefix + "correction"));

        itsSolTab = itsH5Parm.getSolTab(itsSolTabName);
        itsCorrectType = stringToCorrectType(itsSolTab.getType());
        if (itsCorrectType==PHASE && nPol("")==1) {
          itsCorrectType = SCALARPHASE;
        }
        if (itsCorrectType==AMPLITUDE && nPol("")==1) {
          itsCorrectType = SCALARAMPLITUDE;
        }
        itsDirection = 0;
        if (directionStr=="") {
          ASSERT(!itsSolTab.hasAxis("dir") || itsSolTab.getAxis("dir").size==1);
          // If there is only one direction, silently assume it is the right one
        } else if (itsSolTab.hasAxis("dir") && itsSolTab.getAxis("dir").size>1) {
          itsDirection = itsSolTab.getDirIndex(directionStr);
        }
      } else {
        itsTimeSlotsPerParmUpdate =
            parset.isDefined(prefix + "timeslotsperparmupdate") ?
            parset.getInt (prefix + "timeslotsperparmupdate") :
            parset.getInt (defaultPrefix + "timeslotsperparmupdate", 500);
        string correctTypeStr = toLower(
            parset.isDefined(prefix + "correction") ?
            parset.getString(prefix + "correction") :
            parset.getString (defaultPrefix + "correction", "gain"));
        itsCorrectType = stringToCorrectType(correctTypeStr);
      }

      if (itsCorrectType==FULLJONES && itsUpdateWeights) {
        ASSERTSTR (itsInvert, "Updating weights has not been implemented for invert=false and fulljones");
      }
    }

    string OneApplyCal::correctTypeToString(CorrectType ct) {
      if (ct==GAIN) return "gain";
      if (ct==FULLJONES) return "fulljones";
      if (ct==TEC) return "tec";
      if (ct==CLOCK) return "clock";
      if (ct==SCALARPHASE) return "scalarphase";
      if (ct==SCALARAMPLITUDE) return "scalaramplitude";
      if (ct==ROTATIONANGLE) return "rotationangle";
      if (ct==ROTATIONMEASURE) return "rotationmeasure";
      if (ct==PHASE) return "phase";
      if (ct==AMPLITUDE) return "amplitude";
      THROW(Exception, "Unknown correction type: "<<ct);
      return "";
    }

    OneApplyCal::CorrectType OneApplyCal::stringToCorrectType(const string& ctStr) {
      if (ctStr=="gain") return GAIN;
      if (ctStr=="fulljones") return FULLJONES;
      if (ctStr=="tec") return TEC;
      if (ctStr=="clock") return CLOCK;
      if (ctStr=="scalarphase" || ctStr=="commonscalarphase") return SCALARPHASE;
      if (ctStr=="scalaramplitude" || ctStr=="commonscalaramplitude") return SCALARAMPLITUDE;
      if (ctStr=="phase") return PHASE;
      if (ctStr=="amplitude") return AMPLITUDE;
      if (ctStr=="rotationangle" || ctStr=="commonrotationangle" || ctStr=="rotation") return ROTATIONANGLE;
      if (ctStr=="rotationmeasure") return ROTATIONMEASURE;
      THROW(Exception, "Unknown correction type: "<<ctStr);
      return GAIN;
    }

    OneApplyCal::~OneApplyCal()
    {}

    void OneApplyCal::updateInfo (const DPInfo& infoIn)
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

      if (itsUseH5Parm) {
          itsTimeSlotsPerParmUpdate = info().ntime();
      } else { // Use ParmDB
        itsParmDB.reset(new BBS::ParmFacade(itsParmDBName));
      }

      // Detect if full jones solutions are present
      if (!itsUseH5Parm &&
          (itsCorrectType == GAIN || itsCorrectType==FULLJONES) &&
          (itsParmDB->getNames("Gain:0:1:*").size() +
           itsParmDB->getDefNames("Gain:0:1:*").size() >0 )) {
        itsCorrectType=FULLJONES;
      }

      // Detect if solutions are saved as Real/Imag or Ampl/Phase
      if (itsCorrectType == GAIN || itsCorrectType == FULLJONES ){
        if (itsUseH5Parm) {
          // H5Parm uses amplitude / phase by definition
          itsUseAP = true;
        } else {
          // Determine from values present in parmdb what to use
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
      }

      if (itsCorrectType == GAIN) {
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
      } else if (itsCorrectType == FULLJONES) {
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
      } else if (itsCorrectType == TEC) {
        if (nPol("TEC")==1) {
          itsParmExprs.push_back("TEC");
        }
        else {
          itsParmExprs.push_back("TEC:0");
          itsParmExprs.push_back("TEC:1");
        }
      } else if (itsCorrectType == CLOCK) {
        if (nPol("Clock")==1) {
          itsParmExprs.push_back("Clock");
        }
        else {
          itsParmExprs.push_back("Clock:0");
          itsParmExprs.push_back("Clock:1");
        }
      } else if (itsCorrectType == ROTATIONANGLE) {
        itsParmExprs.push_back("{Common,}RotationAngle");
      } else if (itsCorrectType == SCALARPHASE) {
        itsParmExprs.push_back("{Common,}ScalarPhase");
      } else if (itsCorrectType == ROTATIONMEASURE) {
        itsParmExprs.push_back("RotationMeasure");
      } else if (itsCorrectType == SCALARAMPLITUDE) {
        itsParmExprs.push_back("{Common,}ScalarAmplitude");
      } else if (itsCorrectType == PHASE) {
        ASSERT(itsUseH5Parm);
        itsParmExprs.push_back("Phase:0");
        itsParmExprs.push_back("Phase:1");
      } else if (itsCorrectType == AMPLITUDE) {
        ASSERT(itsUseH5Parm);
        itsParmExprs.push_back("Amplitude:0");
        itsParmExprs.push_back("Amplitude:1");
      } else {
        THROW (Exception, "Correction type "<<
                          correctTypeToString(itsCorrectType)<<" is unknown");
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

    void OneApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << std::endl;
      if (itsUseH5Parm) {
        os << "  H5Parm:         " << itsParmDBName << endl;
        os << "    SolTab:       " << itsSolTabName << endl;
      } else {
        os << "  parmdb:         " << itsParmDBName << endl;
      }
      os << "  correction:     " << correctTypeToString(itsCorrectType) << endl;
      if (itsCorrectType==GAIN || itsCorrectType==FULLJONES) {
        os << "    Ampl/Phase:   " << boolalpha << itsUseAP << endl;
      }
      os << "  update weights: " << boolalpha << itsUpdateWeights << endl;
      os << "  invert:         " << boolalpha << itsInvert <<endl;
      if (itsInvert) {
      os << "    sigmaMMSE:    " << itsSigmaMMSE << endl;
      }
      os << "  timeSlotsPerParmUpdate: " << itsTimeSlotsPerParmUpdate <<endl;
    }

    void OneApplyCal::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " OneApplyCal " << itsName << endl;
    }

    bool OneApplyCal::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);

      if (bufin.getTime() > itsLastTime) {
        updateParms(bufin.getTime());
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
            ApplyCal::applyFull( &itsParms(0, antA, timeFreqOffset),
                       &itsParms(0, antB, timeFreqOffset),
                       &data[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &weight[bl * itsNCorr * nchan + chan * itsNCorr ],
                       &flag[  bl * itsNCorr * nchan + chan * itsNCorr ],
                       bl, chan, itsUpdateWeights, itsFlagCounter);
          }
          else {
            ApplyCal::applyDiag( &itsParms(0, antA, timeFreqOffset),
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

    void OneApplyCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


    void OneApplyCal::updateParms (const double bufStartTime)
    {
      uint numAnts = info().antennaNames().size();

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
      itsLastTime = bufStartTime - 0.5*itsTimeInterval + 
                    itsTimeSlotsPerParmUpdate * itsTimeInterval;
      uint numTimes = itsTimeSlotsPerParmUpdate;

      double lastMSTime = info().startTime() + info().ntime() * itsTimeInterval;
      if (itsLastTime > lastMSTime && !nearAbs(itsLastTime, lastMSTime, 1.e-3)) {
        itsLastTime = lastMSTime;
        numTimes = info().ntime() % itsTimeSlotsPerParmUpdate;
      }

      map<string, vector<double> > parmMap;
      map<string, vector<double> >::iterator parmIt;

      uint tfDomainSize=numTimes*numFreqs;

      // Fill parmvalues here, get raw data from H5Parm or ParmDB
      if (itsUseH5Parm) {
#pragma omp critical(updateH5ParmValues)
{
        // TODO: understand polarization etc.
        ASSERT(itsParmExprs.size()==1 || itsParmExprs.size()==2);
        hsize_t startTime = 0;
        if (itsSolTab.hasAxis("time")) {
          startTime = itsSolTab.getTimeIndex(bufStartTime);
        }
        hsize_t startFreq = 0;
        if (itsSolTab.hasAxis("freq")) {
          startFreq = itsSolTab.getFreqIndex(info().chanFreqs()[0]);
        }
        uint freqUpsampleFactor = numFreqs;

        double h5freqinterval = 0.;
        if (itsSolTab.hasAxis("freq") && itsSolTab.getAxis("freq").size > 1) {
          h5freqinterval = itsSolTab.getFreqInterval();
          ASSERT(h5freqinterval>0);
          freqUpsampleFactor = h5freqinterval/info().chanWidths()[0] + 0.5; // Round;
          ASSERTSTR(near(h5freqinterval, freqUpsampleFactor*info().chanWidths()[0],1.e-5),
                    "H5Parm freq interval ("<<h5freqinterval<<") is not an integer " <<
                    "multiple of MS freq interval ("<<info().chanWidths()[0]<<")");
          if (freqUpsampleFactor > numFreqs) {
            freqUpsampleFactor = numFreqs;
          }
        }

        uint timeUpsampleFactor = numTimes;
        if (itsSolTab.hasAxis("time") && itsSolTab.getAxis("time").size > 1) {
          double h5timeInterval = itsSolTab.getTimeInterval();
          timeUpsampleFactor = h5timeInterval/itsTimeInterval+0.5; // Round
          ASSERTSTR(near(h5timeInterval,timeUpsampleFactor*itsTimeInterval,1.e-5),
                    "H5Parm time interval ("<<h5timeInterval<<") is not an integer " <<
                    "multiple of MS time interval ("<<itsTimeInterval<<")");
          if (timeUpsampleFactor > numTimes) {
            timeUpsampleFactor = numTimes;
          }
        }

        // Figure out whether time or frequency is first axis
        bool freqvariesfastest = true;
        if (itsSolTab.hasAxis("freq") && itsSolTab.hasAxis("time") &&
            itsSolTab.getAxisIndex("freq") < itsSolTab.getAxisIndex("time")) {
          freqvariesfastest = false;
        }
        ASSERT(freqvariesfastest);

        // Take the ceiling of numTimes/timeUpsampleFactor, same for freq
        uint numTimesInH5Parm = (numTimes+timeUpsampleFactor-1)/timeUpsampleFactor;
        uint numFreqsInH5Parm = (numFreqs+freqUpsampleFactor-1)/freqUpsampleFactor;

        // Check that frequencies match
        if (itsSolTab.hasAxis("freq") && itsSolTab.getAxis("freq").size > 1) {
          vector<double> h5parmfreqs = itsSolTab.getRealAxis("freq");
          for (uint f=0; f<info().nchan(); ++f) {
            ASSERT(nearAbs(info().chanFreqs()[f],
                           h5parmfreqs[startFreq + f/freqUpsampleFactor],
                           h5freqinterval*0.501));
          }
        }

        for (uint ant = 0; ant < numAnts; ++ant) {
          for (uint pol=0; pol<itsParmExprs.size(); ++pol) {
            vector<double> rawsols, rawweights;
            rawsols = itsSolTab.getValues(info().antennaNames()[ant],
                                        startTime, numTimesInH5Parm, 1,
                                        startFreq, numFreqsInH5Parm, 1, pol, itsDirection);

            rawweights = itsSolTab.getWeights(info().antennaNames()[ant],
                                        startTime, numTimesInH5Parm, 1,
                                        startFreq, numFreqsInH5Parm, 1, pol, itsDirection);

            parmvalues[pol][ant].resize(tfDomainSize);

            size_t tf=0;
            for (uint t=0; t<numTimesInH5Parm; ++t) {
              for (uint ti=0; ti<timeUpsampleFactor; ++ti) {
                for (uint f=0; f<numFreqsInH5Parm; ++f) {
                  for (uint fi=0; fi<freqUpsampleFactor; ++fi) {
                    if (tf<tfDomainSize) {
                      if (rawweights[t*numFreqsInH5Parm + f]>0) {
                        parmvalues[pol][ant][tf++] = rawsols[t*numFreqsInH5Parm + f];
                      } else {
                        parmvalues[pol][ant][tf++] = std::numeric_limits<double>::quiet_NaN();
                      }
                    }
                  }
                }
              }
            }
            ASSERT(tf==tfDomainSize);
          }
        }
} // End pragma omp critical
      } else { // Use ParmDB
        for (uint parmExprNum = 0; parmExprNum<itsParmExprs.size();++parmExprNum) {
          // parmMap contains parameter values for all antennas
          parmMap = itsParmDB->getValuesMap( itsParmExprs[parmExprNum] + "{:phase_center,}*",
                                 minFreq, maxFreq, freqInterval,
                                 bufStartTime-0.5*itsTimeInterval, itsLastTime,
                                 itsTimeInterval, true);

          string parmExpr = itsParmExprs[parmExprNum];

          // Resolve {Common,}Bla to CommonBla or Bla
          if (!parmMap.empty() &&
              parmExpr.find("{Common,}") != string::npos) {
            // Take the name of the first parm, e.g. Bla:CS001, and remove the antenna name
            uint colonPos = (parmMap.begin()->first).find(":");
            parmExpr = (parmMap.begin()->first).substr(0, colonPos);
          }
          
          string name_postfix = "";
          // Remove :phase_center postfix
          if (!parmMap.empty()) {
            // Take the name of the first parm, e.g. Bla:CS001, and remove the antenna name
            // If necessary, append :phase_center
            if ((parmMap.begin()->first).find(":phase_center") != string::npos) {
              name_postfix = ":phase_center";
            }
          }

          for (uint ant = 0; ant < numAnts; ++ant) {
            parmIt = parmMap.find(parmExpr + ":" + string(info().antennaNames()[ant])
                                   + name_postfix);

            if (parmIt != parmMap.end()) {
              parmvalues[parmExprNum][ant].swap(parmIt->second);
              ASSERT(parmvalues[parmExprNum][ant].size()==tfDomainSize);
            } else {// No value found, try default
              Array<double> defValues;
              double defValue;

              if (itsParmDB->getDefValues(parmExpr + ":" +
                  string(info().antennaNames()[ant]) + name_postfix).size()==1) { // Default for antenna
                itsParmDB->getDefValues(string(itsParmExprs[parmExprNum]) + ":" +
                  string(info().antennaNames()[ant]) + name_postfix).get(0,defValues);
                ASSERT(defValues.size()==1);
                defValue=defValues.data()[0];
              }
              else if (itsParmDB->getDefValues(parmExpr).size() == 1) { //Default value
                itsParmDB->getDefValues(parmExpr).get(0,defValues);
                ASSERT(defValues.size()==1);
                defValue=defValues.data()[0];
              } else if (parmExpr.substr(0,5)=="Gain:") {
                defValue=0.;
              }
              else {
                THROW (Exception, "No parameter value found for "+
                   parmExpr + ":" + string(info().antennaNames()[ant]) + name_postfix);
              }

              parmvalues[parmExprNum][ant].resize(tfDomainSize);
              for (uint tf=0; tf<tfDomainSize;++tf) {
                parmvalues[parmExprNum][ant][tf]=defValue;
              }
            }
          }
        }
      }

      ASSERT(parmvalues[0][0].size() <= tfDomainSize); // Catches multiple matches

      double freq;

      // Make parameters complex
      for (uint tf=0;tf<tfDomainSize;++tf) {
        for (uint ant=0;ant<numAnts;++ant) {

          freq=info().chanFreqs()[tf % numFreqs];

          if (itsCorrectType==GAIN) {
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
          else if (itsCorrectType==FULLJONES) {
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
          else if (itsCorrectType==TEC) {
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
          else if (itsCorrectType==CLOCK) {
            itsParms(0, ant, tf)=polar(1.,
                parmvalues[0][ant][tf] * freq * casacore::C::_2pi);
            if (itsParmExprs.size() == 1) { // No Clock:0, only Clock:
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[0][ant][tf] * freq * casacore::C::_2pi);
            }
            else { // Clock:0 and Clock:1
              itsParms(1, ant, tf)=polar(1.,
                  parmvalues[1][ant][tf] * freq * casacore::C::_2pi);
            }
          }
          else if (itsCorrectType==ROTATIONANGLE) {
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
          else if (itsCorrectType==ROTATIONMEASURE) {
            double lambda2 = casacore::C::c / freq;
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
          else if (itsCorrectType==PHASE || itsCorrectType==SCALARPHASE) {
            itsParms(0, ant, tf) = polar(1., parmvalues[0][ant][tf]);
            if (itsCorrectType==SCALARPHASE) { // Same value for x and y
              itsParms(1, ant, tf) = polar(1., parmvalues[0][ant][tf]);
            } else { // Different value for x and y
              itsParms(1, ant, tf) = polar(1., parmvalues[1][ant][tf]);
            }
          }
          else if (itsCorrectType==AMPLITUDE || itsCorrectType==SCALARAMPLITUDE) {
            itsParms(0, ant, tf) = parmvalues[0][ant][tf];
            if (itsCorrectType==SCALARAMPLITUDE) { // Same value for x and y
              itsParms(1, ant, tf) = parmvalues[0][ant][tf];
            } else { // Different value for x and y
              itsParms(1, ant, tf) = parmvalues[1][ant][tf];
            }
          }

          // Invert
          if (itsInvert) {
            if (itsParms.shape()[0]==2) {
            itsParms(0, ant, tf) = 1./itsParms(0, ant, tf);
            itsParms(1, ant, tf) = 1./itsParms(1, ant, tf);
            } else if (itsCorrectType==FULLJONES) {
              ApplyCal::invert(&itsParms(0, ant, tf),itsSigmaMMSE);
            } else {
              ASSERT (itsCorrectType==ROTATIONMEASURE || itsCorrectType==ROTATIONANGLE);
              // rotationmeasure and commonrotationangle are already inverted above
            }
          }
        }
      }
    }

    uint OneApplyCal::nPol(const string& parmName) {
      if (itsUseH5Parm) {
        if (!itsSolTab.hasAxis("pol")) {
          return 1;
        } else {
          return itsSolTab.getAxis("pol").size;
        }
      } else { // Use ParmDB
        if (itsParmDB->getNames(parmName+":0:*").empty() &&
                    itsParmDB->getDefNames(parmName+":0:*").empty() ) {
          return 1;
        } else {
          return 2;
        }
      }
    }

    void OneApplyCal::initDataArrays() {
      uint numAnts=info().antennaNames().size();
      uint tfDomainSize=itsTimeSlotsPerParmUpdate*info().chanFreqs().size();

      uint numParms;
      if (itsCorrectType==FULLJONES ||
          itsCorrectType==ROTATIONANGLE ||
          itsCorrectType==ROTATIONMEASURE) {
        numParms = 4;
      }
      else {
        numParms = 2;
      }

      itsParms.resize(numParms, numAnts, tfDomainSize);
    }

    void OneApplyCal::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by OneApplyCal " << itsName;
      os << endl << "=======================" << endl;
      itsFlagCounter.showBaseline (os, itsCount);
      itsFlagCounter.showChannel  (os, itsCount);
    }

  } //# end namespace
}
