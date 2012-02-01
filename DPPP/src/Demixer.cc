//# Demixer.cc: DPPP step class to subtract A-team sources
//# Copyright (C) 2011
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
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/Demixer.h>
#include <DPPP/Averager.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <ParmDB/Axis.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <Common/OpenMP.h>
#include <BBSKernel/MeasurementAIPS.h>

#include <scimath/Mathematics/MatrixMathLA.h>
#include <casa/Arrays/MatrixIter.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <casa/Quanta/MVAngle.h>

#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    Demixer::Demixer (DPInput* input,
                      const ParSet& parset, const string& prefix)
      : itsInput        (input),
        itsName         (prefix),
        itsTarget       (parset.getString(prefix+"target")),
        itsSources      (parset.getStringVector (prefix+"sources")),
//        itsExtraSources (parset.getStringVector (prefix+"sources")),
        itsJointSolve   (parset.getBool  (prefix+"jointsolve", true)),
        itsNChanAvg     (parset.getUint  (prefix+"freqstep", 1)),
        itsNTimeAvg     (parset.getUint  (prefix+"timestep", 1)),
        itsResChanAvg   (parset.getUint  (prefix+"avgfreqstep", itsNChanAvg)),
        itsResTimeAvg   (parset.getUint  (prefix+"avgtimestep", itsNTimeAvg)),
        itsNTimeChunk   (parset.getUint  (prefix+"ntimechunk", 0)),
        itsNTimeIn      (0),
        itsNTimeOut     (0),

        itsBaselineMask(True),
        itsCorrelationMask(True)
    {
      cout << "PREFIX: " << prefix << endl;

      // Default nr of time chunks is maximum number of threads.
      if (itsNTimeChunk == 0) {
        itsNTimeChunk = OpenMP::maxThreads();
      }
      cout << "NTIMECHUNK: " << itsNTimeChunk << endl;

      // Check that time and freq windows fit nicely.
      ASSERTSTR ((itsNTimeChunk * itsNTimeAvg) % itsResTimeAvg == 0,
                 "time window should fit averaging integrally");

      // JVZ: Added to keep track of time info needed for BBS grid.
      itsTimeCenters.reserve(itsNTimeChunk);
      itsTimeWidths.reserve(itsNTimeChunk);

      // Collect all source names.
      itsNrDir = itsSources.size() + itsExtraSources.size() + 1;
      itsAllSources.reserve (itsNrDir);
      itsAllSources.insert (itsAllSources.end(),
                            itsSources.begin(), itsSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsExtraSources.begin(), itsExtraSources.end());
//      itsAllSources.push_back("target"); //// probably not needed
      itsAllSources.push_back(itsTarget);

      // Size buffers.
      itsFactors.resize (itsNTimeChunk);
//      itsBuf.resize (itsNTimeChunk * itsNTimeAvg);
      itsPhaseShifts.reserve (itsNrDir-1);
      itsFirstSteps.reserve  (itsNrDir);
      itsAvgResults.reserve  (itsNrDir);
//      itsBBSExpr.reserve     (itsNrDir);

      // Create the steps for the sources to be removed.
      // Demixing consists of the following steps:
      // - phaseshift data to each demix source
      // - average data in each direction, also for original phasecenter.
      // - determine demix factors for all directions
      // - use BBS to predict and solve in each direction. It is possible to
      //   predict more directions than to solve (for strong sources in field).
      // - use BBS to subtract the solved sources using the demix factors.
      //   The averaging used here can be smaller than used when solving.
      for (uint i=0; i<itsNrDir-1; ++i) {
        // First make the phaseshift and average steps for each demix source.
        // The resultstep gets the result.
        // The phasecenter can be given in a parameter. Its name is the default.
        // Note the PhaseShift knows about source names CygA, etc.
        itsPhaseShifts.push_back (new PhaseShift
                                  (input, parset,
                                   prefix + itsAllSources[i] + '.',
                                   itsAllSources[i]));
        DPStep::ShPtr step1 (itsPhaseShifts[i]);
        itsFirstSteps.push_back (step1);
        DPStep::ShPtr step2 (new Averager(input, parset, prefix));
        step1->setNextStep (step2);
        MultiResultStep* step3 = new MultiResultStep(itsNTimeChunk);
        step2->setNextStep (DPStep::ShPtr(step3));
        // There is a single demix factor step which needs to get all results.
        itsAvgResults.push_back (step3);
      }
      // Now create the step to average the data themselves.
      DPStep::ShPtr targetAvg(new Averager(input, prefix,
                                           itsNChanAvg, itsNTimeAvg));
      itsFirstSteps.push_back (targetAvg);
      MultiResultStep* targetAvgRes = new MultiResultStep(itsNTimeChunk);
      targetAvg->setNextStep (DPStep::ShPtr(targetAvgRes));
      itsAvgResults.push_back (targetAvgRes);

      // Open ParmDB and SourceDB.
      try {
        itsSourceDB = boost::shared_ptr<SourceDB>
          (new SourceDB(ParmDBMeta("casa", "sky")));
        ParmManager::instance().initCategory(SKY, itsSourceDB->getParmDB());
      } catch (Exception &e) {
        THROW(Exception, "Failed to open sky model parameter database: "
              << "sky");
      }

      try {
        ParmManager::instance().initCategory(INSTRUMENT,
                                             ParmDB(ParmDBMeta("casa",
                                                               "instrument")));
      } catch (Exception &e) {
        THROW(Exception, "Failed to open instrument model parameter database: "
              << "instrument");
      }

      // Create Instrument instance using information present in DPInput.
      size_t nStations = input->antennaNames().size();

//      vector<Station::Ptr> stations;
//      stations.reserve(nStations);
//      for (size_t i = 0; i < nStations; ++i) {
//        // Get station name and ITRF position.
//        casa::MPosition position = MPosition::Convert(input->antennaPos()[i],
//                                                      MPosition::ITRF)();
//        // Store station information.
//        stations.push_back(Station::Ptr(new Station(input->antennaNames()(i),
//                                                    position)));
//      }
//      MPosition position = MPosition::Convert(input->arrayPos(),
//                                              MPosition::ITRF)();
//      Instrument::Ptr instrument(new Instrument("LOFAR", position,
//        stations.begin(), stations.end()));

      MeasurementAIPS ___bla(input->msName());
      Instrument::Ptr instrument = ___bla.instrument();


      // Get directions and make sure they are in J2000.
      MDirection refDelay = MDirection::Convert(input->delayCenter(),
                                               MDirection::J2000)();
      MDirection refTile  = MDirection::Convert(input->tileBeamDir(),
                                               MDirection::J2000)();

      // Construct frequency axis (needs channel width information).
      double chanWidth = input->chanWidths()[0];
      ASSERT(allEQ(input->chanWidths(), chanWidth));
      itsFreqAxisAvg = Axis::ShPtr(new RegularAxis(input->chanFreqs()(0)
        - 0.5 * chanWidth, chanWidth, input->chanFreqs().size()));

//      double factor = input->chanFreqs().size() / itsNChanAvg;
//      if(input->chanFreqs().size() % itsNChanAvg > 0)
//      {
//        ++factor;
//      }

//      LOG_DEBUG_STR("factor: " << factor);
      itsFreqAxisAvg = itsFreqAxisAvg->compress(itsNChanAvg);

      ASSERT(input->getAnt1().size() == input->getAnt2().size());
      for(size_t i = 0; i < input->getAnt1().size(); ++i)
      {
        unsigned int ant1 = input->getAnt1()[i];
        unsigned int ant2 = input->getAnt2()[i];
        itsBaselines.append(baseline_t(ant1, ant2));
      }

      // Deselect auto-correlations.
      for(size_t i = 0; i < nStations; ++i)
      {
          itsBaselineMask.clear(i, i);
      }

      ASSERT(input->ncorr() == 4);
      itsCorrelations.append(Correlation::XX);
      itsCorrelations.append(Correlation::XY);
      itsCorrelations.append(Correlation::YX);
      itsCorrelations.append(Correlation::YY);

      ModelConfig config;
//      config.setDirectionalGain();
      config.setGain();
      BeamConfig beamConfig(BeamConfig::DEFAULT, false,
        casa::Path("$LOFARROOT/share"));
      config.setBeamConfig(beamConfig);
      config.setCache();

      vector<string> incl, excl;
      incl.push_back("Gain:*");
//      incl.push_back("DirectionalGain:*");

      // TODO: Need to know reference frequency!!!
      double refFreq = itsFreqAxisAvg->center(itsFreqAxisAvg->size() / 2);
      for(size_t i = 0; i < itsAllSources.size(); ++i)
      {
        try
        {
          // TODO: Find sane way to derive phase center.
//          MDirection refPhase = MDirection::makeMDirection(itsAllSources[i]);

//          MDirection refPhase;
//          string keyName = prefix+itsAllSources[i]+".phasecenter";
//          if(parset.parameterSet().isDefined(keyName))
//          {
//            refPhase = handleCenter(parset.getStringVector(prefix+itsAllSources[i]+".phasecenter"));
//          }
//          else
//          {
//            refPhase = handleCenter(vector<string>(1, itsAllSources[i]));
//          }

          MDirection refPhase = MDirection::Convert(input->phaseCenter(),
            MDirection::J2000)();

          cout << "PHASE CENTER BBS: " << i << " " << refPhase << endl;
//          config.setSources(vector<string>(1, itsAllSources[i]));
          config.setSources(vector<string>(1, "SB000*"));
          MeasurementExpr::Ptr model(new MeasurementExprLOFAR(*itsSourceDB,
            BufferMap(), config, instrument, itsBaselines, refFreq, refPhase,
            refDelay, refTile));
          itsModels.push_back(model);
        }
        catch(Exception &e)
        {
          THROW(Exception, "Unable to construct model expression for source: "
            << itsAllSources[i] << " (" << e.what() << ")");
        }

        ParmGroup parms = ParmManager::instance().makeSubset(incl, excl,
          itsModels.back()->parms());
        itsModelParms.push_back(parms);
        itsParms.insert(parms.begin(), parms.end());
      }

//      itsParms = ParmManager::instance().makeSubset(incl, excl);

      SolverOptions lsqOptions;
//      lsqOptions.maxIter = 200;
//      lsqOptions.epsValue = 1e-8;
//      lsqOptions.epsDerivative = 1e-8;
//      lsqOptions.colFactor = 1e-6;
//      lsqOptions.lmFactor = 1e-3;
//      lsqOptions.balancedEq = false;
//      lsqOptions.useSVD = true;

      lsqOptions.maxIter = 40;
      lsqOptions.epsValue = 1e-9;
      lsqOptions.epsDerivative = 1e-9;
      lsqOptions.colFactor = 1e-9;
      lsqOptions.lmFactor = 1.0;
      lsqOptions.balancedEq = false;
      lsqOptions.useSVD = true;

      itsOptions = EstimateOptions(EstimateOptions::COMPLEX,
        EstimateOptions::L2, false, 1, false, ~flag_t(0), flag_t(4),
        lsqOptions);

    }

    Demixer::~Demixer()
    {
    }

    MDirection Demixer::handleCenter(const vector<string> &center) const
    {
      // A case-insensitive name can be given for a moving source (e.g. SUN)
      // or a known source (e.g. CygA).
      if (center.size() == 1) {
        return MDirection::makeMDirection (center[0]);
      }
      // The phase center must be given in J2000 as two values (ra,dec).
      // In the future time dependent frames can be done as in UVWFlagger.
      ASSERTSTR (center.size() == 2,
                 "2 values must be given in PhaseShift phasecenter");
      ///ASSERTSTR (center.size() < 4,
      ///"Up to 3 values can be given in UVWFlagger phasecenter");
      MDirection phaseCenter;
      if (center.size() == 1) {
        string str = toUpper(center[0]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid source type"
                   " in UVWFlagger phasecenter");
        return MDirection(tp);
      }
      Quantity q0, q1;
      ASSERTSTR (MVAngle::read (q0, center[0]),
                 center[0] << " is an invalid RA or longitude"
                 " in UVWFlagger phasecenter");
      ASSERTSTR (MVAngle::read (q1, center[1]),
                 center[1] << " is an invalid DEC or latitude"
                 " in UVWFlagger phasecenter");
      MDirection::Types type = MDirection::J2000;
      if (center.size() > 2) {
        string str = toUpper(center[2]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid direction type in UVWFlagger"
                   " in UVWFlagger phasecenter");
      }
      return MDirection(q0, q1, type);
    }

    void Demixer::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();

      itsTimeInterval = info.timeInterval();
      cout << "itsTimeInterval: " << itsTimeInterval << endl;

      itsNrChanIn = info.nchan();
      itsNrBl     = info.nbaselines();
      itsNrCorr   = info.ncorr();
      itsFactorBuf.resize (IPosition(4, itsNrBl, itsNrChanIn, itsNrCorr,
                                     itsNrDir*(itsNrDir-1)/2));
      // Let the internal steps update their data.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infocp;
      for (uint i=0; i<itsFirstSteps.size(); ++i) {
        infocp = info;
        DPStep::ShPtr step = itsFirstSteps[i];
        while (step) {
          step->updateInfo (infocp);
          step = step->getNextStep();
        }
        // Create the BBSexpression.
//        itsBBSExpr.push_back (BBSExpr::ShPtr(new BBSExpr(*itsInput, infocp,
//                                                         itsAllSources[i])));
//        itsModels.push_back (itsBBSExpr[i]->getModel());
      }
      // Keep the averaged time interval.
      itsTimeIntervalAvg = infocp.timeInterval();
      // Update the info of this object.
      info.update (itsResChanAvg, itsResTimeAvg);
      itsNrChanOut = info.nchan();
      itsTimeIntervalRes = info.timeInterval();
    }

    void Demixer::show (std::ostream& os) const
    {
      os << "Demixer " << itsName << std::endl;
      os << "  target:         " << itsTarget << std::endl;
      os << "  sources:        " << itsSources << std::endl;
      os << "  extrasources:   " << itsExtraSources << std::endl;
      os << "  jointsolve:     " << itsJointSolve << std::endl;
      os << "  freqstep:       " << itsNChanAvg << std::endl;
      os << "  timestep:       " << itsNTimeAvg << std::endl;
      os << "  avgfreqstep:    " << itsResChanAvg << std::endl;
      os << "  avgtimestep:    " << itsResTimeAvg << std::endl;
    }

    void Demixer::showTimings (std::ostream& os, double duration) const
    {
      const double self = itsTimer.getElapsed();

      os << "  ";
      FlagCounter::showPerc1 (os, self, duration);
      os << " Demixer " << itsName << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerPhaseShift.getElapsed(), self);
      os << " Phase shift" << endl;

      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerDemix.getElapsed(), self);
      os << " Decorrelation factors" << endl;
    }

    bool Demixer::process (const DPBuffer& buf)
    {
      itsTimer.start();

      if(itsNTimeIn == 0)
      {
        itsTimeStart = buf.getTime() - 0.5 * itsTimeInterval;
      }

      ++itsNTimeIn;

      // Make sure all required data arrays are filled in.
      DPBuffer newBuf(buf);
      RefRows refRows(newBuf.getRowNrs());
      if (newBuf.getUVW().empty()) {
        newBuf.setUVW(itsInput->fetchUVW(newBuf, refRows, itsTimer));
      }
      if (newBuf.getWeights().empty()) {
        newBuf.setWeights(itsInput->fetchWeights(newBuf, refRows, itsTimer));
      }
      if (newBuf.getFullResFlags().empty()) {
        newBuf.setFullResFlags(itsInput->fetchFullResFlags(newBuf, refRows,
                                                            itsTimer));
      }

      // Do the initial steps (phaseshift and average).
      itsTimerPhaseShift.start();
///#pragma omp parallel for
      for (int i=0; i<int(itsFirstSteps.size()); ++i) {
        itsFirstSteps[i]->process(newBuf);
      }
      itsTimerPhaseShift.stop();

      // For each itsNTimeAvg times, calculate the
      // phase rotation per direction.
      itsTimerDemix.start();
      addFactors(newBuf);
      if (itsNTimeIn % itsNTimeAvg == 0) {
        // TODO: NB: This call increases itsNTimeOut as a side effect!!!!
        averageFactors();
        itsNTimeIn  = 0;
      }
      itsTimerDemix.stop();

      // Do BBS solve, etc. when sufficient time slots have been collected.
      if (itsNTimeOut == itsNTimeChunk) {
        demix();

        ASSERT(itsNTimeIn == 0);
        itsNTimeOut = 0;
      }

      itsTimer.stop();
      return true;
    }

    void Demixer::finish()
    {
      // Process remaining entries.
      // Let the next steps finish.
      if (itsNTimeIn > 0) {
        itsTimer.start();

        // Finish the initial steps (phaseshift and average).
        itsTimerPhaseShift.start();
        ///#pragma omp parallel for
        for (int i=0; i<int(itsFirstSteps.size()); ++i) {
          itsFirstSteps[i]->finish();
        }
        itsTimerPhaseShift.stop();

        itsTimerDemix.start();
        // TODO: NB: This call increases itsNTimeOut as a side effect!!!!
        averageFactors();
        itsTimerDemix.stop();

        demix();
        itsTimer.stop();
      }

      getNextStep()->finish();
    }

    void Demixer::addFactors (const DPBuffer& newBuf)
    {
///#pragma omp parallel
      {
///#pragma omp for
        uint ncorr  = newBuf.getData().shape()[0];
        uint nchan  = newBuf.getData().shape()[1];
        uint nbl    = newBuf.getData().shape()[2];
        DComplex* factorPtr = itsFactorBuf.data();
        //# If ever in the future a time dependent phase center is used,
        //# the machine must be reset for each new time, thus each new call
        //# to process.
        for (uint i1=0; i1<itsNrDir-1; ++i1) {
          for (uint i0=i1+1; i0<itsNrDir; ++i0) {
            const double* uvw       = newBuf.getUVW().data();
            const bool*   flagPtr   = newBuf.getFlags().data();
            const float*  weightPtr = newBuf.getWeights().data();
            const DComplex* phasor1 = itsPhaseShifts[i1]->getPhasors().data();
            if (i0 == itsNrDir-1) {
              for (uint i=0; i<nbl; ++i) {
                for (uint j=0; j<nchan; ++j) {
                  DComplex factor = conj(*phasor1++);
                  for (uint k=0; k<ncorr; ++k) {
                    if (! *flagPtr) {
                      *factorPtr += factor * double(*weightPtr);
                    }
                    flagPtr++;
                    weightPtr++;
                    factorPtr++;
                  }
                }
                uvw += 3;
              }
            } else {
              const DComplex* phasor0 = itsPhaseShifts[i0]->getPhasors().data();
              for (uint i=0; i<nbl; ++i) {
                for (uint j=0; j<nchan; ++j) {
                  // Probably multiply with conj
                  DComplex factor = *phasor0++ / *phasor1++;
                  for (uint k=0; k<ncorr; ++k) {
                    if (! *flagPtr) {
                      *factorPtr += factor * double(*weightPtr);
                    }
                    flagPtr++;
                    weightPtr++;
                    factorPtr++;
                  }
                }
                uvw += 3;
              }
            }
          }
        }
      }
    }

    void Demixer::averageFactors()
    {
      // The averaged weights are calculated in the Averager, so use those.
      const Cube<float>& weightSums = itsAvgResults[0]->get()[itsNTimeOut].getWeights();
      ASSERT (! weightSums.empty());
      itsFactors[itsNTimeOut].resize (IPosition(5, itsNrDir, itsNrDir,
                                                itsNrCorr, itsNrChanOut,
                                                itsNrBl));
      itsFactors[itsNTimeOut] = DComplex(1,0);
      const DComplex* phin = itsFactorBuf.data();
      for (uint d0=0; d0<itsNrDir; ++d0) {
        for (uint d1=d0+1; d1<itsNrDir; ++d1) {
          DComplex* ph1 = itsFactors[itsNTimeOut].data() + d0*itsNrDir + d1;
          DComplex* ph2 = itsFactors[itsNTimeOut].data() + d1*itsNrDir + d0;
          // Average for all channels and divide by the summed weights.
          const float* weightPtr = weightSums.data();
          for (uint k=0; k<itsNrBl; ++k) {
            for (uint c0=0; c0<itsNrChanOut; ++c0) {
              DComplex sum[4];
              uint nch = std::min(itsNChanAvg, itsNrChanIn-c0*itsNChanAvg);
              for (uint c1=0; c1<nch; ++c1) {
                for (uint j=0; j<itsNrCorr; ++j) {
                  sum[j] += *phin++;
                }
              }
              for (uint j=0; j<itsNrCorr; ++j) {
                *ph1 = sum[j] / double(*weightPtr++);
                *ph2 = conj(*ph1);
                ph1 += itsNrDir*itsNrDir;
                ph2 += itsNrDir*itsNrDir;
              }
            }
          }
        }
      }
      ///      cout << "factor=" <<itsFactors[itsNTimeOut] << endl;
      // Clear the summation buffer.
      itsFactorBuf = DComplex();

      // TODO: Factor this out somehow??
      unsigned int nTime = itsNTimeIn % itsNTimeAvg;
      nTime = (nTime == 0 ? itsNTimeAvg : nTime);
      itsTimeWidths.push_back(nTime * itsTimeInterval);
      itsTimeCenters.push_back(itsTimeStart + 0.5 * itsTimeWidths.back());

      itsNTimeOut++;
    }

    void Demixer::demix()
    {
      // Collect buffers for each direction.
      vector<vector<DPBuffer> > buffers;
      for(uint i = 0; i < itsAvgResults.size(); ++i)
      {
        buffers.push_back(itsAvgResults[i]->get());
        itsAvgResults[i]->clear();
      }

      // Solve for the gains in the various directions.
      for(uint i = 0; i < itsModels.size(); ++i)
      {
        itsModels[i]->setSolvables(itsModelParms[i]);
      }

      // Make time axis based on averaged target visibilities.
      ASSERT(itsTimeCenters.size() == itsNTimeOut
        && itsTimeWidths.size() == itsNTimeOut);

      Axis::ShPtr timeAxis(new OrderedAxis(itsTimeCenters, itsTimeWidths));
      itsTimeCenters.clear();
      itsTimeWidths.clear();

      // Make time axis and grid.
      Grid visGrid(itsFreqAxisAvg, timeAxis);

      // Solve for each time slot over all channels.
      Grid solGrid(itsFreqAxisAvg->compress(itsFreqAxisAvg->size()), timeAxis);

      LOG_DEBUG_STR("SHAPES: " << itsFactors[0].shape() << " " << itsFreqAxisAvg->size() << " " << buffers[0][0].getData().shape());

//      double startTime = itsBuf[0].getTime() - itsInput->timeInterval() * 0.5;
//      Axis::ShPtr timeAxis(new RegularAxis(startTime, itsTimeIntervalAvg,
//                                             itsNTimeOut));
//      Grid grid(itsBBSExpr[0]->getFreqAxis(), timeAxis);

      // Set parameter domain.
      ParmManager::instance().setDomain(solGrid.getBoundingBox());

      // Estimate model parameters.
      estimate(buffers, itsModels, itsFactors, itsBaselines, itsCorrelations,
        itsBaselineMask, itsCorrelationMask, visGrid, solGrid, itsOptions);

      // Flush solutions to disk.
      ParmManager::instance().flush();


      // Subtract the demixed sources.
      // TODO: As soon as we allow a different output resolution for the target
      // field, visGrid needs to be re-derived for the subtract.
      ASSERT(itsNChanAvg == itsResChanAvg && itsNTimeAvg == itsResTimeAvg);
      LOG_DEBUG_STR("subtracting....");

      // Solve for the gains in the various directions.
      for(uint i = 0; i < itsModels.size(); ++i)
      {
        itsModels[i]->clearSolvables();
      }

      vector<unsigned int> directions(itsSources.size());
      for(size_t i = 0; i < itsSources.size(); ++i)
      {
        directions[i] = i;
      }

      const unsigned int target = itsAllSources.size() - 1;

      LOG_DEBUG_STR("target: " << target << " directions: " << directions);

      subtract(buffers.back(), itsModels, itsFactors, itsBaselines,
        itsCorrelations, itsBaselineMask, itsCorrelationMask, visGrid, target,
        directions);

      // Let the next steps process the data.
      // TODO: As soon as we allow a different output resolution for the target
      // field, this has to be adapted.
      ASSERT(itsNChanAvg == itsResChanAvg && itsNTimeAvg == itsResTimeAvg);
      for(uint i = 0; i < buffers.back().size(); ++i)
      {
        getNextStep()->process(buffers.back()[i]);
      }

      // Clear the intermediate buffers.
      // TODO: Could already clear all buffers except for the target field
      // before flushing buffers down the pipline (above).
      for(uint i = 0; i < itsAvgResults.size(); ++i)
      {
        itsAvgResults[i]->clear();
      }
    }

//    void Demixer::subtract()
//    {
//      // Set expressions to not solvable.
//      for (uint i=0; i<itsModels.size(); ++i) {
//        itsModels[i]->clearSolvables();
//      }
//      // Loop through all time windows.
//      for (uint i=0; i<itsNTimeOut; ++i) {
//        // Subtract data for each time window.
//      }
//    }



      // array dim: baseline, time, freq, dir1, dir2
      // NDPPP compares data with predictions in Estimate class
      // Beam Info lezen en aan BBS doorgeven (evt. via BBS class)
      //   Joris maakt zijn read functie public
      // Get RA and DEC of phase center (take care of moving targets).
      /// MeqExpr for FFT/degrid per baseline met apart (gedeeld) degrid object
      ///     MeqMatrix[2][2] degrid (uvCoordStat1, uvCoordStat2, times, freqs)
      /// Joris:
      //  - publiek maken beam info lezen
      //  - estimate functie voor demixing
      /// Ger
      //  - create demixing matrix
      //  - multiple predict result with demix matrix (also derivatives)
      //  - subtract mbv demixing matrix
      //  awimager
      //  - for subbands in awimager use correct reffreq
      //    both for separate windows and for subbands combined in single band
      //    moet parameter worden in ATerm::evaluate
      //  - option to treat channels in spw as subbands
      //  Sven
      //  - on-the-fly degridding
      //  Ronald
      //  - Can Solver solve in block matrices?
      //  Wim
      //  - can solve be parallelized?
      // Cyril:
      // - separate time window for element beam?
      // - gridding element beam because small conv.func?
      // - degridding: apply element beam per time window?
      // - commit the code
      // - Joris: new beam model in imager (in Cyril's branch)
      // - Bas: merge ionosphere in imager (gridding takes most time)
      // - Sanjay: write paper wide-band MSMFS
      // -         paper wide-band A-projection
      //    Cyril: A-projection or LOFAR plus element beam trick
      //    Bas:   ionosphere
      //    Johan/Stefan: verify beam model
      //    George: is MSMFS needed in awimager?

    // awimager: default channels by MFS is channel 0 only.

  } //# end namespace
}
