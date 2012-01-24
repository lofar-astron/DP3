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

#include <scimath/Mathematics/MatrixMathLA.h>
#include <casa/Arrays/MatrixIter.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>

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
        itsSources      (parset.getStringVector (prefix+"sources")),
        itsExtraSources (parset.getStringVector (prefix+"sources")),
        itsJointSolve   (parset.getBool  (prefix+"jointsolve", true)),
        itsNChanAvg     (parset.getUint  (prefix+"freqstep", 1)),
        itsNTimeAvg     (parset.getUint  (prefix+"timestep", 1)),
        itsResChanAvg   (parset.getUint  (prefix+"avgfreqstep", itsNChanAvg)),
        itsResTimeAvg   (parset.getUint  (prefix+"avgtimestep", itsNTimeAvg)),
        itsNTimeChunk   (parset.getUint  (prefix+"ntimechunk", 0)),
        itsNTimeIn      (0),
        itsNTimeOut     (0)
    {
      // Default nr of time chunks is maximum number of threads.
      if (itsNTimeChunk == 0) {
        itsNTimeChunk = OpenMP::maxThreads();
      }
      // Check that time and freq windows fit nicely.
      ASSERTSTR ((itsNTimeChunk * itsNTimeAvg) % itsResTimeAvg == 0,
                 "time window should fit averaging integrally");
      // Collect all source names.
      itsNrDir = itsSources.size() + itsExtraSources.size() + 1;
      itsAllSources.reserve (itsNrDir);
      itsAllSources.insert (itsAllSources.end(),
                            itsSources.begin(), itsSources.end());
      itsAllSources.insert (itsAllSources.end(),
                            itsExtraSources.begin(), itsExtraSources.end());
      /// itsAllSources.push_back ("target"); //// probably not needed
      // Size buffers.
      itsFactors.resize (itsNTimeChunk);
      itsBuf.resize (itsNTimeChunk * itsNTimeAvg);
      itsPhaseShifts.reserve (itsNrDir-1);
      itsFirstSteps.reserve  (itsNrDir);
      itsAvgResults.reserve  (itsNrDir);
      itsBBSExpr.reserve     (itsNrDir);
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
                                   prefix + '.' + itsAllSources[i],
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
    }

    Demixer::~Demixer()
    {}

    void Demixer::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
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
        itsBBSExpr.push_back (BBSExpr::ShPtr(new BBSExpr(*itsInput, infocp, 
                                                         itsAllSources[i])));
        itsModels.push_back (itsBBSExpr[i]->getModel());
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
      double timing = itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, timing, duration);
      os << " Demixer " << itsName << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerPhaseShift.getElapsed(), timing);
      os << " of it spent in phase shifting data" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsTimerDemix.getElapsed(), timing);
      os << " of it spent in calculating demix factors" << endl;
    }

    bool Demixer::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Keep the buffer.
      DPBuffer& newBuf = itsBuf[itsNTimeIn++];
      newBuf = buf;
      // Make sure all required data arrays are filled in.
      RefRows refRows(newBuf.getRowNrs());
      if (newBuf.getUVW().empty()) {
        newBuf.setUVW (itsInput->fetchUVW(newBuf, refRows, itsTimer));
      }
      if (newBuf.getWeights().empty()) {
        newBuf.setWeights (itsInput->fetchWeights(newBuf, refRows, itsTimer));
      }
      if (newBuf.getFullResFlags().empty()) {
        newBuf.setFullResFlags (itsInput->fetchFullResFlags(newBuf, refRows,
                                                            itsTimer));
      }
      // Do the initial steps (phaseshift and average).
      itsTimerPhaseShift.start();
///#pragma omp parallel for
      for (int i=0; i<int(itsFirstSteps.size()); ++i) {
        itsFirstSteps[i]->process (newBuf);
      }
      itsTimerPhaseShift.stop();
      // For each itsNTimeAvg times, calculate the
      // phase rotation per direction.
      itsTimerDemix.start();
      addFactors (newBuf);
      if (itsNTimeIn % itsNTimeAvg == 0) {
        averageFactors();
      }
      itsTimerDemix.stop();
      itsTimer.stop();
      // Do BBS solve, etc. when sufficient time slots have been collected.
      if (itsNTimeOut == itsNTimeChunk) {
        cout << "process time chunks" << endl;
        demix();
        itsNTimeIn  = 0;
        itsNTimeOut = 0;
      }
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
        averageFactors();
        itsTimerDemix.stop();
        itsTimer.stop();
        cout << "final process time chunks" << endl;
        demix();
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
      itsNTimeOut++;
    }

    void Demixer::demix()
    {
      // Solve for the gains in the various directions.
      for (uint i=0; i<itsModels.size(); ++i) {
        /// itsModels[i]->setSolvables();  //What to put in ParmGroup?
      }
      // Make time axis and grid.
      double startTime = itsBuf[0].getTime() - itsInput->timeInterval() * 0.5;
      Axis::ShPtr timeAxis (new RegularAxis (startTime, itsTimeIntervalAvg,
                                             itsNTimeOut));
      Grid grid(itsBBSExpr[0]->getFreqAxis(), timeAxis);
      // estimate (dpbuffers, exprs, grid, baselineMask, ...);
      // 
      // Subtract the demixed sources.
      subtract();
      // Clear the input buffers (to cut in memory usage).
      for (uint i=0; i<itsNTimeIn; ++i) {
        itsBuf[i].clear();
      }
      // Let the next steps process the data.
      for (uint i=0; i<itsNTimeOut*itsNTimeAvg / itsResTimeAvg; ++i) {
        ///        getNextStep()->process (buf2);
      }
      // Clear the result buffers.
      for (uint i=0; i<itsAvgResults.size(); ++i) {
        itsAvgResults[i]->clear();
      }
    }

    void Demixer::subtract()
    {
      // Set expressions to not solvable.
      for (uint i=0; i<itsModels.size(); ++i) {
        itsModels[i]->clearSolvables();
      }
      // Loop through all time windows.
      for (uint i=0; i<itsNTimeOut; ++i) {
        // Subtract data for each time window.
      }
    }



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
