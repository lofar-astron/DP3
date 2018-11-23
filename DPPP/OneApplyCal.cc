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

#include "OneApplyCal.h"

#include "ApplyCal.h"
#include "Exceptions.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "MSReader.h"

#include "../Common/ParallelFor.h"
#include "../Common/ParameterSet.h"
#include "../Common/StringUtil.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/OS/File.h>

#include <iostream>
#include <limits>
#include <algorithm>
#include <iomanip>

#include <boost/algorithm/string/case_conv.hpp>

using namespace casacore;
using namespace DP3::BBS;

/// Look at BBSKernel MeasurementExprLOFARUtil.cc and Apply.cc

namespace DP3 {
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
        itsSolSetName  (
            parset.isDefined(prefix + "solset") ?
            parset.getString(prefix + "solset") :
            parset.getString(defaultPrefix + "solset", "")),
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

      assert (!itsParmDBName.empty());

      if (substep) {
        itsInvert=false;
      } else {
        itsInvert = (parset.isDefined(prefix + "invert") ?
                     parset.getBool (prefix + "invert") :
                     parset.getBool (defaultPrefix + "invert", true));
      }

      if (itsUseH5Parm) {
        string interpolationStr = (parset.isDefined(prefix + "interpolation") ?
                            parset.getString(prefix + "interpolation") :
                            parset.getString(prefix + "interpolation", "nearest"));
        if (interpolationStr == "nearest") {
          itsInterpolationType = InterpolationType::NEAREST;
        } else if (interpolationStr == "linear") {
          itsInterpolationType = InterpolationType::LINEAR;
        } else {
          throw std::runtime_error("Unsupported interpolation mode: " + interpolationStr);
        }
        itsTimeSlotsPerParmUpdate = 0;
        string directionStr;
        directionStr = (parset.isDefined(prefix + "direction") ?
                        parset.getString(prefix + "direction") :
                        parset.getString(defaultPrefix + "direction",
                          predictDirection));
        itsH5Parm = H5Parm(itsParmDBName, false, false, itsSolSetName);

        itsSolTabName = (parset.isDefined(prefix + "correction") ?
                         parset.getString(prefix + "correction") :
                         parset.getString(defaultPrefix + "correction"));
        if(itsSolTabName == "fulljones")
        {
          std::vector<string> solTabs = parset.getStringVector(prefix + "soltab", std::vector<string>{"amplitude000", "phase000"});
          if(solTabs.size() != 2)
            throw std::runtime_error("The soltab parameter requires two soltabs for fulljones correction (amplitude and phase)");
          itsSolTab = itsH5Parm.getSolTab(solTabs[0]);
          itsSolTab2 = itsH5Parm.getSolTab(solTabs[1]);
          itsSolTabName = solTabs[0] + ", " + solTabs[1]; // this is only so that show() shows these tables
          itsCorrectType = FULLJONES;
        }
        else {
          itsSolTab = itsH5Parm.getSolTab(itsSolTabName);
          itsCorrectType = stringToCorrectType(itsSolTab.getType());
        }
        if (itsCorrectType==PHASE && nPol("")==1) {
          itsCorrectType = SCALARPHASE;
        }
        if (itsCorrectType==AMPLITUDE && nPol("")==1) {
          itsCorrectType = SCALARAMPLITUDE;
        }
        itsDirection = 0;
        if (directionStr=="") {
          assert(!itsSolTab.hasAxis("dir") || itsSolTab.getAxis("dir").size==1);
          // If there is only one direction, silently assume it is the right one
        } else if (itsSolTab.hasAxis("dir") && itsSolTab.getAxis("dir").size>1) {
          itsDirection = itsSolTab.getDirIndex(directionStr);
        }
      } else {
        itsTimeSlotsPerParmUpdate =
            parset.isDefined(prefix + "timeslotsperparmupdate") ?
            parset.getInt (prefix + "timeslotsperparmupdate") :
            parset.getInt (defaultPrefix + "timeslotsperparmupdate", 500);
        string correctTypeStr = boost::to_lower_copy(
            parset.isDefined(prefix + "correction") ?
            parset.getString(prefix + "correction") :
            parset.getString (defaultPrefix + "correction", "gain"));
        itsCorrectType = stringToCorrectType(correctTypeStr);
      }

      if (itsCorrectType==FULLJONES && itsUpdateWeights) {
        if (!itsInvert)
          throw Exception("Updating weights has not been implemented for invert=false and fulljones");
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
      throw Exception("Unknown correction type: " + std::to_string(ct));
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
      throw Exception("Unknown correction type: " + ctStr);
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

      assert(itsNCorr==4);

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
            throw Exception("No gains found in parmdb "+itsParmDBName);
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
        assert(itsUseH5Parm);
        itsParmExprs.push_back("Phase:0");
        itsParmExprs.push_back("Phase:1");
      } else if (itsCorrectType == AMPLITUDE) {
        assert(itsUseH5Parm);
        itsParmExprs.push_back("Amplitude:0");
        itsParmExprs.push_back("Amplitude:1");
      } else {
        throw Exception("Correction type " +
                          correctTypeToString(itsCorrectType) + " is unknown");
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

        if (!itsUseH5Parm) {
          if(!regularChannels)
            throw Exception(
                    "ApplyCal with parmdb requires evenly spaced channels.");
        }
      }
    }

    void OneApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << '\n';
      if (itsUseH5Parm) {
        os << "  H5Parm:         " << itsParmDBName << '\n';
        os << "    SolSet:       " << itsH5Parm.getSolSetName() << '\n';
        os << "    SolTab:       " << itsSolTabName << '\n';
        os << "  Direction:      " << itsDirection << '\n';
        os << "  Interpolation:  " << (itsInterpolationType==InterpolationType::NEAREST?"nearest":"linear") << '\n';
      } else {
        os << "  Parmdb:         " << itsParmDBName << '\n';
      }
      os << "  Correction:     " << correctTypeToString(itsCorrectType) << '\n';
      if (itsCorrectType==GAIN || itsCorrectType==FULLJONES) {
        os << "    Ampl/Phase:   " << boolalpha << itsUseAP << '\n';
      }
      os << "  Update weights: " << boolalpha << itsUpdateWeights << '\n';
      os << "  Invert:         " << boolalpha << itsInvert <<'\n';
      if (itsInvert) {
      os << "    SigmaMMSE:    " << itsSigmaMMSE << '\n';
      }
      os << "  TimeSlotsPerParmUpdate: " << itsTimeSlotsPerParmUpdate <<'\n';
    }

    void OneApplyCal::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " OneApplyCal " << itsName << '\n';
    }

    bool OneApplyCal::process (const DPBuffer& bufin, std::mutex* hdf5Mutex)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);

      if (bufin.getTime() > itsLastTime) {
        updateParms(bufin.getTime(), hdf5Mutex);
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

      ParallelFor<size_t> loop(getInfo().nThreads());
      loop.Run(0, nbl, [&](size_t bl, size_t /*thread*/) {
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
      });

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

    void OneApplyCal::applyFlags(vector<double>& values,
                                 const vector<double>& weights) {
      assert(values.size() == weights.size());
      vector<double>::iterator values_it = values.begin();
      vector<double>::const_iterator weights_it = weights.begin();

      for (; values_it != values.end(); ++values_it) {
        if (*weights_it == 0.) {
          *values_it = std::numeric_limits<float>::quiet_NaN();
        }
        weights_it++;
      } 
    }

    void OneApplyCal::updateParms (const double bufStartTime, std::mutex* hdf5Mutex)
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

      std::map<std::string, std::vector<double> > parmMap;
      std::map<std::string, std::vector<double> >::iterator parmIt;

      uint tfDomainSize=numTimes*numFreqs;

      // Fill parmvalues here, get raw data from H5Parm or ParmDB
      if (itsUseH5Parm) {
        std::unique_lock<std::mutex> lock;
        if(hdf5Mutex != nullptr)
          lock = std::unique_lock<std::mutex>(*hdf5Mutex);
        
        // TODO: understand polarization etc.
        //  assert(itsParmExprs.size()==1 || itsParmExprs.size()==2);

        // Figure out whether time or frequency is first axis
        bool freqvariesfastest = true;
        if (itsSolTab.hasAxis("freq") && itsSolTab.hasAxis("time") &&
            itsSolTab.getAxisIndex("freq") < itsSolTab.getAxisIndex("time")) {
          freqvariesfastest = false;
        }
        assert(freqvariesfastest);

        vector<double> times(info().ntime());
        for (uint t=0; t<times.size(); ++t) {
          // time centroids
          times[t] = info().startTime() + (t+0.5) * info().timeInterval();
        }
        vector<double> freqs(info().chanFreqs().size());
        for (uint ch=0; ch<info().chanFreqs().size(); ++ch) {
          freqs[ch] = info().chanFreqs()[ch];
        }

        vector<double> weights;
        for (uint ant = 0; ant < numAnts; ++ant) {
          if(itsCorrectType == FULLJONES)
          {
            for (uint pol=0; pol<4; ++pol) {
              // Place amplitude in even and phase in odd elements
              parmvalues[pol*2][ant] = itsSolTab.getValuesOrWeights("val",
                info().antennaNames()[ant],
                times, freqs,
                pol, itsDirection, itsInterpolationType==InterpolationType::NEAREST);
              weights = itsSolTab.getValuesOrWeights("weight",
                info().antennaNames()[ant], times, freqs, pol, itsDirection,
                itsInterpolationType==InterpolationType::NEAREST);
              applyFlags(parmvalues[pol*2][ant], weights);
              parmvalues[pol*2+1][ant] = itsSolTab2.getValuesOrWeights("val",
                info().antennaNames()[ant],
                times, freqs,
                pol, itsDirection, itsInterpolationType==InterpolationType::NEAREST);
              weights = itsSolTab2.getValuesOrWeights("weight",
                info().antennaNames()[ant], times, freqs, pol, itsDirection,
                itsInterpolationType==InterpolationType::NEAREST);
              applyFlags(parmvalues[pol*2+1][ant], weights);
            }
          }
          else {
            for (uint pol=0; pol<itsParmExprs.size(); ++pol) {
              parmvalues[pol][ant] = itsSolTab.getValuesOrWeights("val",
                info().antennaNames()[ant],
                times, freqs,
                pol, itsDirection, itsInterpolationType==InterpolationType::NEAREST);
              weights = itsSolTab.getValuesOrWeights("weight",
                info().antennaNames()[ant], times, freqs, pol, itsDirection,
                itsInterpolationType==InterpolationType::NEAREST);
              applyFlags(parmvalues[pol][ant], weights);
            }
          }
        }
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
              assert(parmvalues[parmExprNum][ant].size()==tfDomainSize);
            } else {// No value found, try default
              Array<double> defValues;
              double defValue;

              if (itsParmDB->getDefValues(parmExpr + ":" +
                  string(info().antennaNames()[ant]) + name_postfix).size()==1) { // Default for antenna
                itsParmDB->getDefValues(string(itsParmExprs[parmExprNum]) + ":" +
                  string(info().antennaNames()[ant]) + name_postfix).get(0,defValues);
                assert(defValues.size()==1);
                defValue=defValues.data()[0];
              }
              else if (itsParmDB->getDefValues(parmExpr).size() == 1) { //Default value
                itsParmDB->getDefValues(parmExpr).get(0,defValues);
                assert(defValues.size()==1);
                defValue=defValues.data()[0];
              } else if (parmExpr.substr(0,5)=="Gain:") {
                defValue=0.;
              }
              else {
                throw Exception("No parameter value found for "+
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

      assert(parmvalues[0][0].size() <= tfDomainSize); // Catches multiple matches

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
              assert (itsCorrectType==ROTATIONMEASURE || itsCorrectType==ROTATIONANGLE);
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
      os << "\nFlags set by OneApplyCal " << itsName;
      os << "\n=======================\n";
      itsFlagCounter.showBaseline (os, itsCount);
      itsFlagCounter.showChannel  (os, itsCount);
    }

  } //# end namespace
}
