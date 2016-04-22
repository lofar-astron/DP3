//# AOFlaggerStep.cc: DPPP step class to flag data based on rficonsole
//# Copyright (C) 2010
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
//# $Id: AOFlaggerStep.cc 31423 2015-04-03 14:06:21Z dijkema $
//#
//# @author Andre Offringa, Ger van Diepen

#include <lofar_config.h>
#include <DPPP_AOFlag/AOFlaggerStep.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <Common/OpenMP.h>

#include <casa/OS/HostInfo.h>
#include <casa/OS/File.h>

#include <aoflagger.h>

#include <iostream>
#include <algorithm>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    AOFlaggerStep::AOFlaggerStep (DPInput* input,
                                  const ParameterSet& parset,
                                  const string& prefix)
      : itsName        (prefix),
        itsBufIndex    (0),
        itsNTimes      (0),
        itsMemoryNeeded(0),
        itsFlagCounter (input->msName(), parset, prefix+"count."),
        itsMoveTime    (0),
        itsFlagTime    (0),
        itsQualTime    (0),
        itsRfiStats    ()
    {
      itsStrategyName = parset.getString (prefix+"strategy", string());
      itsWindowSize   = parset.getUint   (prefix+"timewindow", 0);
      itsMemory       = parset.getUint   (prefix+"memorymax", 0);
      itsMemoryPerc   = parset.getUint   (prefix+"memoryperc", 0);
      itsOverlap      = parset.getUint   (prefix+"overlapmax", 0);
      // Also look for keyword overlap for backward compatibility.
      if (itsOverlap == 0) {
        itsOverlap    = parset.getUint   (prefix+"overlap", 0);
      }
      itsOverlapPerc  = parset.getDouble (prefix+"overlapperc", -1);
      itsPulsarMode   = parset.getBool   (prefix+"pulsar", false);
      itsPedantic     = parset.getBool   (prefix+"pedantic", false);
      itsDoAutoCorr   = parset.getBool   (prefix+"autocorr", true);
      itsDoRfiStats   = parset.getBool   (prefix+"keepstatistics", true);
    }

    AOFlaggerStep::~AOFlaggerStep()
    {}

    DPStep::ShPtr AOFlaggerStep::makeStep (DPInput* input,
                                           const ParameterSet& parset,
                                           const std::string& prefix)
    {
      return DPStep::ShPtr(new AOFlaggerStep(input, parset, prefix));
    }

    void AOFlaggerStep::show (std::ostream& os) const
    {
      os << "AOFlaggerStep " << itsName << std::endl;
      os << "  strategy:       " << itsStrategyName << std::endl;
      os << "  timewindow:     " << itsWindowSize << std::endl;
      os << "  overlap:        " << itsOverlap << std::endl;
      os << "  pulsar:         " << itsPulsarMode << std::endl;
      os << "  pedantic:       " << itsPedantic << std::endl;
      os << "  keepstatistics: " << itsDoRfiStats << std::endl;
      os << "  autocorr:       " << itsDoAutoCorr << std::endl;
      os << "  nthreads (omp)  " << OpenMP::maxThreads() << std::endl;
      os << "  max memory used ";
      formatBytes(os, itsMemoryNeeded);
      os << std::endl;
    }

    void AOFlaggerStep::formatBytes(std::ostream& os, double bytes) {
      int exp=0;
      while (bytes >= 1024 && exp<5) {
        bytes/=1024;
        exp++;
      }

      uint origPrec=os.precision();
      os.precision(1);

      if (exp==0) {
        os<<fixed<<bytes<<" "<<"B";
      } else {
        os<<fixed<<bytes<<" "<<"KMGTPE"[exp-1]<<"B";
      }

      os.precision(origPrec);
    }

    void AOFlaggerStep::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteFlags();
      // Get nr of threads.
      uint nthread = OpenMP::maxThreads();
      // Determine available memory.
      double availMemory = HostInfo::memoryTotal() * 1024.;
      // Determine how much memory can be used.
      double memoryMax = itsMemory * 1024*1024*1024;
      double memory    = memoryMax;
      if (itsMemoryPerc > 0) {
        memory = itsMemoryPerc * availMemory / 100.;
        if (memoryMax > 0  &&  memory > memoryMax) {
          memory = memoryMax;
        }
      } else if (itsMemory <= 0) {
        // Nothing given, so use available memory on this machine.
        // Set 50% (max 2 GB) aside for other purposes.
        memory = availMemory - std::min(0.5 * availMemory, 2.*1024*1024*1024);
      }
      // Determine how much buffer space is needed per time slot.
      // The flagger needs 3 extra work buffers (data+flags) per thread.
      double timeSize = (sizeof(Complex) + sizeof(bool)) *
        (infoIn.nbaselines() + 3*nthread) * infoIn.nchan() * infoIn.ncorr();
      // If no overlap percentage is given, set it to 1%.
      if (itsOverlapPerc < 0  &&  itsOverlap == 0) {
        itsOverlapPerc = 1;
      }
      // If no time window given, determine it from the available memory.
      if (itsWindowSize == 0) {
        double nt = memory / timeSize;
        if (itsOverlapPerc > 0) {
          // Determine the overlap (add 0.5 for rounding).
          // If itsOverLap is also given, it is the maximum.
          double tw = nt / (1 + 2*itsOverlapPerc/100);
          uint overlap = uint(itsOverlapPerc*tw/100 + 0.5);
          if (itsOverlap == 0  ||  overlap < itsOverlap) {
            itsOverlap = overlap;
          }
        }
        itsWindowSize = uint(std::max(1., nt-2*itsOverlap));
        // Make the window size divide the nr of times nicely (if known).
        // In that way we cannot have a very small last window.
        if (infoIn.ntime() > 0) {
          uint nwindow = 1 + (infoIn.ntime() - 1) / itsWindowSize;
          itsWindowSize = 1 + (infoIn.ntime() - 1) / nwindow;
          if (itsOverlapPerc > 0) {
            uint overlap = uint(itsOverlapPerc*itsWindowSize/100 + 0.5);
            if (overlap < itsOverlap) {
              itsOverlap = overlap;
            }
          }
        }
      }
      if (itsOverlap == 0) {
        itsOverlap = uint(itsOverlapPerc*itsWindowSize/100);
      }
      // Check if it all fits in memory.
      itsMemoryNeeded = (itsWindowSize + 2*itsOverlap) * timeSize;
      ASSERTSTR (itsMemoryNeeded < availMemory,
                 "Timewindow " << itsWindowSize
                 << " and/or overlap " << itsOverlap
                 << ' ' << memory
                 << " too large for available memory " << availMemory);
      // Size the buffer (need overlap on both sides).
      itsBuf.resize (itsWindowSize + 2*itsOverlap);
      // Initialize the flag counters.
      itsFlagCounter.init (getInfo());
      itsFreqs = infoIn.chanFreqs();
      // Fill the strategy (used by all threads)
      // (A thread does not need a private strategy; it is
      // safe to share one among different threads.)
      fillStrategy();
    }

    void AOFlaggerStep::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by AOFlaggerStep " << itsName;
      os << endl << "===========================" << endl;
      itsFlagCounter.showBaseline (os, itsNTimes);
      itsFlagCounter.showChannel  (os, itsNTimes);
      itsFlagCounter.showCorrelation (os, itsNTimes);
    }

    void AOFlaggerStep::showTimings (std::ostream& os, double duration) const
    {
      double flagDur = itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, flagDur, duration);
      os << " AOFlaggerStep " << itsName << endl;
      os << "          ";
      // move time and flag time are sum of all threads.
      // Scale them to a single elapsed time.
      double factor = (itsComputeTimer.getElapsed() /
                       (itsMoveTime + itsFlagTime + itsQualTime));
      FlagCounter::showPerc1 (os, itsMoveTime*factor, flagDur);
      os << " of it spent in shuffling data" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, itsFlagTime*factor, flagDur);
      os << " of it spent in calculating flags" << endl;
      if (itsDoRfiStats) {
        os << "          ";
        FlagCounter::showPerc1 (os, itsQualTime*factor +  itsQualityTimer.getElapsed(),
                                flagDur);
        os << " of it spent in making quality statistics" << endl;
      }
    }

    // Alternative strategy is to flag in windows
    //  0 ..  n+2m
    //  n .. 2n+2m
    // 2n .. 3n+2m  etc.
    // and also update the flags in the overlaps
    bool AOFlaggerStep::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Accumulate in the time window until the window and overlap are full. 
      itsNTimes++;
      itsBuf[itsBufIndex].copy (buf);
      ++itsBufIndex;
      if (itsBufIndex == itsWindowSize+2*itsOverlap) {
        flag(2*itsOverlap);
      }
      itsTimer.stop();
      return true;
    }

    void AOFlaggerStep::finish()
    {
      cerr << "  " << itsBufIndex << " time slots to finish in AOFlaggerStep ..."
           << endl;
      itsTimer.start();
      // Set window size to all entries left.
      itsWindowSize = itsBufIndex;
      if (itsWindowSize > 0) {
        // Flag the remaining time slots (without right overlap).
        flag (0);
      }
      itsBuf.clear();
      itsTimer.stop();
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

    void AOFlaggerStep::addToMS (const string& msName)
    {
      itsTimer.start();
      if (itsDoRfiStats) {
        itsQualityTimer.start();
	itsAOFlagger.WriteStatistics(*itsRfiStats, msName);
        itsQualityTimer.stop();
      }
      itsTimer.stop();
      getPrevStep()->addToMS (msName);
    }

    void AOFlaggerStep::flag (uint rightOverlap)
    {
      // Get the sizes of the axes.
      // Note: OpenMP 2.5 needs signed iteration variables.
      int  nrbl   = itsBuf[0].getData().shape()[2];
      uint ncorr  = itsBuf[0].getData().shape()[0];
      ASSERTSTR (ncorr==4, "AOFlaggerStep can only handle all 4 correlations");
      // Get antenna numbers in case applyautocorr is true.
      const Vector<Int>& ant1 = getInfo().getAnt1();
      const Vector<Int>& ant2 = getInfo().getAnt2();
      itsComputeTimer.start();
      // Now flag each baseline for this time window.
      // The baselines can be processed in parallel.
#pragma omp parallel
      {
        // Create thread-private counter object.
        FlagCounter counter;
        counter.init (getInfo());
	
        // Create a statistics object for all polarizations.
	std::vector<double> scanTimes(itsBuf.size());
	for (size_t i=0; i<itsBuf.size(); ++i) {
	  scanTimes[i] = itsBuf[i].getTime();
        }
	aoflagger::QualityStatistics rfiStats =
	  itsAOFlagger.MakeQualityStatistics (scanTimes.data(),
                                              scanTimes.size(),
                                              itsFreqs.data(),
                                              itsFreqs.size(),
                                              4, false);   // no histograms

        // The for loop can be parallellized. This must be done dynamically,
        // because the execution times of iterations can vary.
#pragma omp for schedule(dynamic)
        // GCC-4.3 only supports OpenMP 2.5 that needs signed iteration
        // variables.
        for (int ib=0; ib<nrbl; ++ib) {
          // Do autocorrelations only if told so.
          if (ant1[ib] == ant2[ib]) {
            if (itsDoAutoCorr) {
              flagBaseline (0, itsWindowSize+rightOverlap, 0, ib,
                            counter, rfiStats);
            }
          } else {
            flagBaseline (0, itsWindowSize+rightOverlap, 0, ib,
                          counter, rfiStats);
          }
        } // end of OMP for
#pragma omp critical(aorflagger_updatecounts)
        {
          // Add the counters to the overall object.
          itsFlagCounter.add (counter);
          if (itsDoRfiStats) {
            itsQualityTimer.stop();
            // Add the rfi statistics to the global object.
            if (itsRfiStats == 0) {
              itsRfiStats.reset(new aoflagger::QualityStatistics(rfiStats));
            } else {
              (*itsRfiStats) += rfiStats;
            }
            itsQualityTimer.start();
          }
        }
      } // end of OMP parallel
      itsComputeTimer.stop();
      itsTimer.stop();
      // Let the next step process the buffers.
      // If possible, discard the buffer processed to minimize memory usage.
      for (uint i=0; i<itsWindowSize; ++i) {
        getNextStep()->process (itsBuf[i]);
        ///        itsBuf[i] = DPBuffer();
        ///cout << "cleared buffer " << i << endl;
      }
      itsTimer.start();
      // Shift the buffers still needed to the beginning of the vector.
      // This is a bit easier than keeping a wrapped vector.
      // Note it is a cheap operation, because shallow copies are made.
      for (uint i=0; i<rightOverlap; ++i) {
        itsBuf[i].copy (itsBuf[i+itsWindowSize]);
        ///cout << "moved buffer " <<i+itsWindowSize<<" to "<< i << endl;
      }
      itsBufIndex = rightOverlap;
    }

    void AOFlaggerStep::flagBaseline (uint leftOverlap, uint windowSize,
                                      uint rightOverlap, uint bl,
                                      FlagCounter& counter,
                                      aoflagger::QualityStatistics& rfiStats)
    {
      NSTimer moveTimer, flagTimer, qualTimer;
      moveTimer.start();
      // Get the sizes of the axes.
      uint ntime  = leftOverlap + windowSize + rightOverlap;
      uint nchan  = itsBuf[0].getData().shape()[1];
      uint blsize = nchan * itsBuf[0].getData().shape()[0];
      // Fill the rficonsole buffers and flag.
      // Create the objects for the real and imaginary data of all corr.
      aoflagger::ImageSet imageSet =
	itsAOFlagger.MakeImageSet(ntime, nchan, 8);
      aoflagger::FlagMask origFlags =
	itsAOFlagger.MakeFlagMask(ntime, nchan);
      const uint iStride = imageSet.HorizontalStride();
      const uint fStride = origFlags.HorizontalStride();
      for (uint i=0; i<ntime; ++i) {
        const Complex* data = itsBuf[i].getData().data()  + bl*blsize;
        const bool*   flags = itsBuf[i].getFlags().data() + bl*blsize;
        for (uint j=0; j<nchan; ++j) {
	  for (uint p=0; p!=4; ++p) {
	    imageSet.ImageBuffer(p*2  )[i + j*iStride] = data->real();
	    imageSet.ImageBuffer(p*2+1)[i + j*iStride] = data->imag();
	    data++;
	  }
	  origFlags.Buffer()[i + j*fStride] = *flags;
          flags += 4;
        }
      }
      // Execute the strategy to do the flagging.
      moveTimer.stop();
      flagTimer.start();
      aoflagger::FlagMask rfiMask = itsAOFlagger.Run(*itsStrategy, imageSet);
      flagTimer.stop();
      // Put back the true flags and count newly set flags.
      moveTimer.start();
      for (uint i=leftOverlap; i<windowSize+leftOverlap; ++i) {
        bool* flags = itsBuf[i].getFlags().data() + bl*blsize;
        for (uint j=0; j<nchan; ++j) {
          // Only set if not already set.
          // If any corr is newly set, set all corr.
          if (! flags[0]) {
            bool setFlag = true;
            if (rfiMask.Buffer()[i + j*fStride]) {
              counter.incrCorrelation(0);
              counter.incrCorrelation(1);
              counter.incrCorrelation(2);
              counter.incrCorrelation(3);
            } else {
              setFlag = false;
            }
            if (setFlag) {
              counter.incrBaseline(bl);
              counter.incrChannel(j);
              for (int k=0; k<4; ++k) {
                flags[k] = true;
              }
            }
          }
          flags += 4;
        }
      }
      moveTimer.stop();
      // Update the RFI statistics if needed.
      if (itsDoRfiStats) {
        qualTimer.start();
        addStats (rfiStats, imageSet, rfiMask, origFlags, bl);
        qualTimer.stop();
      }
#pragma omp critical(aorflagger_updatetimers)
      {
        // Add the timings.
        itsMoveTime += moveTimer.getElapsed();
        itsFlagTime += flagTimer.getElapsed();
        itsQualTime += qualTimer.getElapsed();
      } // end of OMP critical
    }

    void AOFlaggerStep::addStats (aoflagger::QualityStatistics& rfiStats,
                                  const aoflagger::ImageSet& values,
                                  const aoflagger::FlagMask& rfiMask, const aoflagger::FlagMask& origMask,
                                  int bl)
    {
      itsAOFlagger.CollectStatistics(rfiStats, values, rfiMask, origMask,
                                     getInfo().getAnt1()[bl], getInfo().getAnt2()[bl]);
    }

    void AOFlaggerStep::fillStrategy ()
    {
      if (! itsStrategyName.empty()) {
        File file(itsStrategyName);
        if (! file.exists()) {
          file = File("$LOFARROOT/share/rfistrategies/" + itsStrategyName);
          if (! file.exists()) {
            THROW (Exception, "Unknown rfistrategy file " << itsStrategyName);
          }
        }
	itsStrategy.reset(new aoflagger::Strategy
                          (itsAOFlagger.LoadStrategy(file.path().absoluteName())));
      } else {
        double centralFrequency = 0.5*(itsFreqs[0] + itsFreqs[itsFreqs.size()-1]);
        double timeRes;
        if (itsBuf.size() >=2 ) {
          timeRes = (itsBuf[itsBuf.size()-1].getTime() -
                     itsBuf[0].getTime()) / (itsBuf.size()-1);
        } else {
          timeRes = 0.0;
        }
        double frequencyRes;
        if (itsFreqs.size() >= 2) {
          frequencyRes = (itsFreqs[itsFreqs.size()-1] -
                          itsFreqs[0]) / (itsFreqs.size()-1);
        } else {
          frequencyRes = 0.0;
        }				
        itsStrategy.reset(new aoflagger::Strategy
                          (itsAOFlagger.MakeStrategy(aoflagger::LOFAR_TELESCOPE,
                                                     aoflagger::StrategyFlags::NONE,
                                                     centralFrequency,
                                                     timeRes,
                                                     frequencyRes)));
      }
    }

  } //# end namespace
}
