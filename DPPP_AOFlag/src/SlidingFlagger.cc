//# SlidingFlagger.cc: DPPP step class to flag data using a sliding window
//# Copyright (C) 2015
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
//# $Id: SlidingFlagger.cc 30697 2015-01-14 12:17:14Z diepen $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP_AOFlag/SlidingFlagger.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <Common/OpenMP.h>

#include <casa/OS/File.h>
#include <casa/OS/HostInfo.h>

#include <AOFlagger/msio/image2d.h>
#include <AOFlagger/msio/mask2d.h>
#include <AOFlagger/msio/timefrequencydata.h>
#include <AOFlagger/strategy/actions/changeresolutionaction.h>
#include <AOFlagger/strategy/actions/combineflagresultsaction.h>
#include <AOFlagger/strategy/actions/foreachcomplexcomponentaction.h>
#include <AOFlagger/strategy/actions/foreachpolarisationaction.h>
#include <AOFlagger/strategy/actions/frequencyselectionaction.h>
#include <AOFlagger/strategy/actions/iterationaction.h>
#include <AOFlagger/strategy/actions/setflaggingaction.h>
#include <AOFlagger/strategy/actions/setimageaction.h>
#include <AOFlagger/strategy/actions/slidingwindowfitaction.h>
#include <AOFlagger/strategy/actions/statisticalflagaction.h>
#include <AOFlagger/strategy/actions/strategyaction.h>
#include <AOFlagger/strategy/actions/sumthresholdaction.h>
#include <AOFlagger/strategy/actions/timeselectionaction.h>
#include <AOFlagger/strategy/control/artifactset.h>
#include <AOFlagger/strategy/control/strategyreader.h>

#include <iostream>
#include <algorithm>

using namespace casa;
using namespace rfiStrategy;

namespace LOFAR {
  namespace DPPP {

    SlidingFlagger::SlidingFlagger (DPInput* input,
                                    const ParameterSet& parset,
                                    const string& prefix)
      : itsName        (prefix),
        itsBufIndex    (0),
        itsNTimes      (0),
        itsNThreads    (OpenMP::maxThreads()),
        itsMemoryNeeded(0),
        itsFlagCounter (input->msName(), parset, prefix+"count.")
    {
      itsStrategyName = parset.getString (prefix+"strategy", string());
      itsWindowSize   = parset.getUint   (prefix+"timewindow", 0);
      itsBufferSize   = parset.getUint   (prefix+"buffersize", 0);
      itsPulsarMode   = parset.getBool   (prefix+"pulsar", false);
      itsPedantic     = parset.getBool   (prefix+"pedantic", false);
      itsDoAutoCorr   = parset.getBool   (prefix+"autocorr", true);
      ASSERT (itsWindowSize > 0  &&  itsBufferSize > 0);
      if (itsBufferSize > itsWindowSize) {
        THROW (Exception, "timewindow < buffersize, so AOFlagger should "
               "be used instead of SlidingFlagger");
      }
    }

    SlidingFlagger::~SlidingFlagger()
    {}

    DPStep::ShPtr SlidingFlagger::makeStep (DPInput* input,
                                            const ParameterSet& parset,
                                            const std::string& prefix)
    {
      return DPStep::ShPtr(new SlidingFlagger(input, parset, prefix));
    }

    void SlidingFlagger::show (std::ostream& os) const
    {
      os << "SlidingFlagger " << itsName << std::endl;
      os << "  strategy:       " << itsStrategyName << std::endl;
      os << "  timewindow:     " << itsWindowSize << std::endl;
      os << "  buffersize:     " << itsBufferSize << std::endl;
      os << "  pulsar:         " << itsPulsarMode << std::endl;
      os << "  pedantic:       " << itsPedantic << std::endl;
      os << "  autocorr:       " << itsDoAutoCorr << std::endl;
      os << "  nthreads (omp)  " << itsNThreads << std::endl;
      os << "  memory used     " << itsMemoryNeeded << std::endl;
    }

    void SlidingFlagger::updateInfo (const DPInfo& infoIn)
    {
      ASSERTSTR (infoIn.ncorr()==4,
                 "SlidingFlagger can only handle all 4 correlations");
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteFlags();
      // Determine available memory.
      double availMemory = HostInfo::memoryTotal() * 1024.;
      // Determine how much buffer space is needed per time slot.
      // The flagger needs 3 extra work buffers (data+flags) per thread.
      double imgSize = 8*sizeof(num_t) * itsWindowSize * infoIn.nchan();
      double bufSize = (infoIn.nbaselines() * infoIn.nchan() * 4 *
                        (sizeof(Complex) + sizeof(bool)));
      // Check if it all fits in memory.
      itsMemoryNeeded = (infoIn.nbaselines() * imgSize +
                         itsBufferSize * bufSize);
      ASSERTSTR (itsMemoryNeeded < availMemory,
                 "timewindow " << itsWindowSize
                 << " and buffersize " << itsBufferSize
                 << ' ' << itsMemoryNeeded
                 << " too large for available memory " << availMemory);
      // Create the AOFlagger objects for each baseline.
      // Fill them with 0 data and True flags.
      ///const Vector<Int>& ant1 = getInfo().getAnt1();
      ///const Vector<Int>& ant2 = getInfo().getAnt2();
      ///int nbl   = ant1.size();
      ///int nchan = getInfo().nchan();
      // Create the rficonsole Image2D objects for each baseline to flag.
      // The axes are freq,time (time varies slowest) to make the shifting
      // in flagBaseline easier. Note that AORFlagger uses time,freq, but
      // the order does not matter for rficonsole.
      ///itsBLData.resize (nbl);
      ///for (int i=0; i<nbl; ++i) {
      ///if (itsDoAutoCorr  ||  ant1[i] != ant2[i]) {
      ///}
      ///}
      // Fill the strategy (used by all threads)
      // (A thread does not need a private strategy; it is
      // safe to share one among different threads.)
      fillStrategy();
      itsFreqs = infoIn.chanFreqs();
      // Create the counters for all possible threads.
      itsTD.resize (itsNThreads);
      for (uint i=0; i<itsNThreads; ++i) {
        itsTD[i].flagCounter.init (getInfo());
      }
      itsBuf.resize (itsBufferSize);
      itsFlagCounter.init (getInfo());
    }

    void SlidingFlagger::showCounts (std::ostream& os) const
    {
      for (size_t i=0; i<itsTD.size(); ++i) {
        itsFlagCounter.add (itsTD[i].flagCounter);
      }
      os << endl << "Flags set by SlidingFlagger " << itsName;
      os << endl << "============================" << endl;
      itsFlagCounter.showBaseline (os, itsNTimes);
      itsFlagCounter.showChannel  (os, itsNTimes);
      itsFlagCounter.showCorrelation (os, itsNTimes);
    }

    void SlidingFlagger::showTimings (std::ostream& os, double duration) const
    {
      double moveTime = 0;
      double flagTime = 0;
      for (size_t i=0; i<itsNThreads; ++i) {
        moveTime += itsTD[i].moveTimer.getElapsed();
        flagTime += itsTD[i].flagTimer.getElapsed();
      }
      double flagDur = itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, flagDur, duration);
      os << " SlidingFlagger " << itsName << endl;
      os << "          ";
      // move time and flag time are sum of all threads.
      // Scale them to a single elapsed time.
      double factor = itsComputeTimer.getElapsed() / (moveTime + flagTime);
      FlagCounter::showPerc1 (os, moveTime*factor, flagDur);
      os << " of it spent in shuffling data" << endl;
      os << "          ";
      FlagCounter::showPerc1 (os, flagTime*factor, flagDur);
      os << " of it spent in calculating flags" << endl;
    }

    bool SlidingFlagger::process (const DPBuffer& buf)
    {
      itsTimer.start();
      itsNTimes++;
      if (itsBufferSize == 1) {
        itsBuf[0].referenceFilled (buf);
        itsBufIndex = 1;
      } else {
        itsBuf[itsBufIndex++].copy (buf);
      }
      if (itsBufIndex == itsBufferSize) {
        flag();
      }
      itsTimer.stop();
      return true;
    }

    void SlidingFlagger::finish()
    {
      cerr << "  " << itsBufIndex << " time slots to finish in SlidingFlagger ..."
           << endl;
      itsTimer.start();
      if (itsBufIndex > 0) {
        // Flag the remaining time slots.
        flag();
      }
      itsTimer.stop();
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

    void SlidingFlagger::flag()
    {
      // Get the sizes of the axes.
      // Note: OpenMP 2.5 needs signed iteration variables.
      int nrbl = itsBuf[0].getData().shape()[2];
      // Get antenna numbers in case applyautocorr is true.
      const Vector<Int>& ant1 = getInfo().getAnt1();
      const Vector<Int>& ant2 = getInfo().getAnt2();
      itsComputeTimer.start();
      // Now flag each baseline for this time window.
      // The baselines can be processed in parallel.
#pragma omp parallel
      {
        // The for loop can be parallellized. This must be done dynamically,
        // because the execution times of iterations can vary.
#pragma omp for schedule(dynamic)
        for (int ib=0; ib<nrbl; ++ib) {
          // Flag baseline only if told so.
          // Do autocorrelations only if told so.
          ///if (itsBLData[ib].realXX != 0) {
          if (ant1[ib] == ant2[ib]) {
            if (itsDoAutoCorr) {
              flagBaseline (ib, itsTD[OpenMP::threadNum()]);
            }
          } else {
            flagBaseline (ib, itsTD[OpenMP::threadNum()]);
          }
        } // end of OMP for
      } // end of OMP parallel
      itsComputeTimer.stop();
      itsTimer.stop();
      // Let the next step process the buffers.
      for (uint i=0; i<itsBufIndex; ++i) {
        getNextStep()->process (itsBuf[i]);
      }
      itsTimer.start();
      itsBufIndex = 0;
    }

    void SlidingFlagger::flagBaseline (uint bl, ThreadData& td)
    {
      td.moveTimer.start();
      // Get the sizes of the axes.
      uint ntime  = itsWindowSize;
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
	  for (uint p=0; p<4; ++p) {
	    imageSet.ImageBuffer(p*2  )[i + j*iStride] = data->real();
	    imageSet.ImageBuffer(p*2+1)[i + j*iStride] = data->imag();
	    data++;
	  }
	  origFlags.Buffer()[i + j*fStride] = *flags;
          flags += 4;
        }
      }
      // Execute the strategy to do the flagging.
      td.moveTimer.stop();
      td.flagTimer.start();
      aoflagger::FlagMask rfiMask = itsAOFlagger.Run(*itsStrategy, imageSet);
      td.flagTimer.stop();
      // Put back the true flags and count newly set flags.
      td.moveTimer.start();
      for (uint i=0; i<itsBufIndex; ++i) {
        bool* flags = itsBuf[i].getFlags().data() + bl*blsize;
        for (uint j=0; j<nchan; ++j) {
          // Only set if not already set.
          // Note that if first corr flag is true, all are true.
          // If any corr is newly set, set all corr.
          if (! flags[0]) {
            bool setFlag = true;
            if (rfiMask.Buffer()[i + j*fStride]) {
              td.flagCounter.incrCorrelation(0);
              td.flagCounter.incrCorrelation(1);
              td.flagCounter.incrCorrelation(2);
              td.flagCounter.incrCorrelation(3);
            } else {
              setFlag = false;
            }
            if (setFlag) {
              td.flagCounter.incrBaseline(bl);
              td.flagCounter.incrChannel(j);
              for (int k=0; k<4; ++k) {
                flags[k] = true;
              }
            }
          }
          flags += 4;
        }
      }
      td.moveTimer.stop();
    }

    void SlidingFlagger::fillStrategy()
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
