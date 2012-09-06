//# AORFlagger.cc: DPPP step class to flag data based on rficonsole
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
//# $Id$
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/AORFlagger.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>

#include <casa/OS/HostInfo.h>
#include <casa/OS/File.h>

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
#include <AOFlagger/strategy/actions/sumthresholdaction.h>
#include <AOFlagger/strategy/actions/timeselectionaction.h>
#include <AOFlagger/strategy/control/artifactset.h>
#include <AOFlagger/strategy/control/strategyreader.h>
#include <AOFlagger/quality/qualitytablesformatter.h>

#include <Common/StreamUtil.h>
#include <Common/LofarLogger.h>
#include <Common/OpenMP.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Containers/Record.h>
#include <casa/Containers/RecordField.h>
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/RecordGram.h>
#include <iostream>
#include <algorithm>

using namespace casa;
using namespace rfiStrategy;

namespace LOFAR {
  namespace DPPP {

    AORFlagger::AORFlagger (DPInput* input,
                            const ParSet& parset, const string& prefix)
      : itsInput       (input),
        itsName        (prefix),
        itsBufIndex    (0),
        itsNTimes      (0),
        itsNTimesToDo  (0),
        itsMemoryNeeded(0),
        itsFlagCounter (input->msName(), parset, prefix+"count."),
        itsMoveTime    (0),
        itsFlagTime    (0),
        itsQualTime    (0),
        itsRfiStats    (4)
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
      // Fill the strategy for all possible threads.
      fillStrategy (itsStrategy);
    }

    AORFlagger::~AORFlagger()
    {}

    void AORFlagger::show (std::ostream& os) const
    {
      os << "AOFlagger " << itsName << std::endl;
      os << "  timewindow:     " << itsWindowSize << std::endl;
      os << "  overlap:        " << itsOverlap << std::endl;
      os << "  pulsar:         " << itsPulsarMode << std::endl;
      os << "  pedantic:       " << itsPedantic << std::endl;
      os << "  keepstatistics: " << itsDoRfiStats << std::endl;
      os << "  autocorr:       " << itsDoAutoCorr << std::endl;
      os << "  nthreads (omp)  " << OpenMP::maxThreads() << std::endl;
      os << "  max memory used " << itsMemoryNeeded << std::endl;
    }

    void AORFlagger::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();
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
      itsNTimesToDo = infoIn.ntime();
      // Size the RFI statistics.
      itsFreqs = infoIn.chanFreqs();
      if (itsDoRfiStats) {
        for (int i=0; i<4; ++i) {
          itsRfiStats.InitializeBand (0, itsFreqs.data(), itsFreqs.size());
        }
      }
    }

    void AORFlagger::showCounts (std::ostream& os) const
    {
      os << endl << "Flags set by AOFlagger " << itsName;
      os << endl << "=======================" << endl;
      itsFlagCounter.showBaseline (os, itsNTimes);
      itsFlagCounter.showChannel  (os, itsNTimes);
      itsFlagCounter.showCorrelation (os, itsNTimes);
    }

    void AORFlagger::showTimings (std::ostream& os, double duration) const
    {
      double flagDur = itsTimer.getElapsed();
      os << "  ";
      FlagCounter::showPerc1 (os, flagDur, duration);
      os << " AOFlagger " << itsName << endl;
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
    bool AORFlagger::process (const DPBuffer& buf)
    {
      itsTimer.start();
      // Accumulate in the time window until the window and overlap are full. 
      itsNTimes++;
      ///      cout<<"inserted at " << itsBufIndex<<endl;
      itsBuf[itsBufIndex++] = buf;
      if (itsBufIndex == itsWindowSize+2*itsOverlap) {
        flag (2*itsOverlap);
      }
      itsTimer.stop();
      return true;
    }

    void AORFlagger::finish()
    {
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

    void AORFlagger::addToMS (const string& msName)
    {
      itsTimer.start();
      if (itsDoRfiStats) {
	itsQualityTimer.start();
        QualityTablesFormatter qualityData(msName);
        itsRfiStats.Save (qualityData);
	itsQualityTimer.stop();
      }
      itsTimer.stop();
    }

    void AORFlagger::flag (uint rightOverlap)
    {
      // Get the sizes of the axes.
      // Note: OpenMP 2.5 needs signed iteration variables.
      int  nrbl   = itsBuf[0].getData().shape()[2];
      uint ncorr  = itsBuf[0].getData().shape()[0];
      ASSERTSTR (ncorr==4, "AOFlagger can only handle all 4 correlations");
      // Get antenna numbers in case applyautocorr is true.
      const Vector<Int>& ant1 = getInfo().getAnt1();
      const Vector<Int>& ant2 = getInfo().getAnt2();
      itsComputeTimer.start();
      // Now flag each baseline for this time window.
      // The baselines can be processed in parallel.
#pragma omp parallel
      {
	// Create thread-private counter object.
        FlagCounter counter (itsFlagCounter);
	// Create thread-private strategy object.
        boost::shared_ptr<Strategy> strategy;
	fillStrategy (strategy);
        // Create a statistics object for all polarizations.
        StatisticsCollection rfiStats(4);
        if (itsDoRfiStats) {
          rfiStats.InitializeBand (0, itsFreqs.data(), itsFreqs.size());
        }
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
                            counter, *strategy, rfiStats);
            }
          } else {
            flagBaseline (0, itsWindowSize+rightOverlap, 0, ib,
                          counter, *strategy, rfiStats);
          }
        } // end of OMP for
#pragma omp critical(aorflagger_updatecounts)
        {
          // Add the counters to the overall object.
          itsFlagCounter.add (counter);
          if (itsDoRfiStats) {
	    itsQualityTimer.stop();
            // Add the rfi statistics to the global object.
            itsRfiStats.Add (rfiStats);
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
        itsBuf[i] = DPBuffer();
        ///cout << "cleared buffer " << i << endl;
      }
      itsTimer.start();
      // Shift the buffers still needed to the beginning of the vector.
      // This is a bit easier than keeping a wrapped vector.
      // Note it is a cheap operation, because shallow copies are made.
      for (uint i=0; i<rightOverlap; ++i) {
        itsBuf[i] = itsBuf[i+itsWindowSize];
        ///cout << "moved buffer " <<i+itsWindowSize<<" to "<< i << endl;
      }
      itsBufIndex = rightOverlap;
    }

    void AORFlagger::flagBaseline (uint leftOverlap, uint windowSize,
                                   uint rightOverlap, uint bl,
                                   FlagCounter& counter,
				   Strategy& strategy,
                                   StatisticsCollection& rfiStats)
    {
      NSTimer moveTimer, flagTimer, qualTimer;
      moveTimer.start();
      // Get the sizes of the axes.
      uint ntime  = leftOverlap + windowSize + rightOverlap;
      uint nchan  = itsBuf[0].getData().shape()[1];
      uint blsize = nchan * itsBuf[0].getData().shape()[0];
      // Fill the rficonsole buffers and flag.
      // Create the objects for the real and imaginary data of all corr.
      Image2DPtr realXX = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr imagXX = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr realXY = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr imagXY = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr realYX = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr imagYX = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr realYY = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Image2DPtr imagYY = Image2D::CreateUnsetImagePtr(ntime, nchan);
      Mask2DPtr origFlags = Mask2D::CreateUnsetMaskPtr(ntime, nchan);
      for (uint i=0; i<ntime; ++i) {
        const Complex* data = itsBuf[i].getData().data()  + bl*blsize;
        const bool*   flags = itsBuf[i].getFlags().data() + bl*blsize;
        for (uint j=0; j<nchan; ++j) {
          realXX->SetValue (i, j, data->real());
          imagXX->SetValue (i, j, data->imag());
          data++;
          realXY->SetValue (i, j, data->real());
          imagXY->SetValue (i, j, data->imag());
          data++;
          realYX->SetValue (i, j, data->real());
          imagYX->SetValue (i, j, data->imag());
          data++;
          realYY->SetValue (i, j, data->real());
          imagYY->SetValue (i, j, data->imag());
          data++;
          *(origFlags->ValuePtr(i, j)) = *flags;
          flags += 4;
	}
      }
      Mask2DCPtr falseMask = Mask2D::CreateSetMaskPtr<false> (ntime, nchan);
      Image2DCPtr zeroData = Image2D::CreateZeroImagePtr (ntime, nchan);
      // Create original data.
      TimeFrequencyData origData(realXX, imagXX, realXY, imagXY,
                                 realYX, imagYX, realYY, imagYY);
      origData.SetIndividualPolarisationMasks (falseMask, falseMask,
                                               falseMask, falseMask);
      // Create contaminated data.
      TimeFrequencyData contData(origData);
      // Create revised data.
      TimeFrequencyData revData(zeroData, zeroData, zeroData, zeroData,
                                zeroData, zeroData, zeroData, zeroData);
      revData.SetIndividualPolarisationMasks (falseMask, falseMask,
                                              falseMask, falseMask);
      ////      boost::mutex mutex;
      ////      ArtifactSet artifacts(&mutex);
      // Create and fill the artifact set. A mutex is not needed.
      ArtifactSet artifacts(0);
      artifacts.SetOriginalData (origData);
      artifacts.SetContaminatedData (contData);
      artifacts.SetRevisedData (revData);
      // Execute the strategy to do the flagging.
      moveTimer.stop();
      flagTimer.start();
      strategy.Perform (artifacts, itsProgressListener);
      flagTimer.stop();
      // Put back the true flags and count newly set flags.
      moveTimer.start();
      Mask2DCPtr maskXX = artifacts.ContaminatedData().GetMask (XXPolarisation);
      Mask2DCPtr maskXY = artifacts.ContaminatedData().GetMask (XYPolarisation);
      Mask2DCPtr maskYX = artifacts.ContaminatedData().GetMask (YXPolarisation);
      Mask2DCPtr maskYY = artifacts.ContaminatedData().GetMask (YYPolarisation);
      for (uint i=leftOverlap; i<windowSize+leftOverlap; ++i) {
        bool* flags = itsBuf[i].getFlags().data() + bl*blsize;
        for (uint j=0; j<nchan; ++j) {
          // Only set if not already set.
          // Note that if first corr flag is true, all are true.
          // If any corr is newly set, set all corr.
          if (! flags[0]) {
            bool setFlag = true;
            if (maskXX->Value(i,j)) {
              counter.incrCorrelation(0);
            } else if (maskXY->Value(i,j)) {
              counter.incrCorrelation(1);
            } else if (maskYX->Value(i,j)) {
              counter.incrCorrelation(2);
            } else if (maskYY->Value(i,j)) {
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
        addStats (rfiStats, realXX, imagXX, maskXX, origFlags, bl, 0);
        addStats (rfiStats, realXY, imagXY, maskXY, origFlags, bl, 1);
        addStats (rfiStats, realYX, imagYX, maskYX, origFlags, bl, 2);
        addStats (rfiStats, realYY, imagYY, maskYY, origFlags, bl, 3);
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

    void AORFlagger::addStats (StatisticsCollection& rfiStats,
                               const Image2DPtr& reals, const Image2DPtr& imags,
			       const Mask2DCPtr& mask, const Mask2DPtr& origFlags,
			       int bl, uint polarization)
    {
      uint nchan = reals->Height();
      uint imagestride = reals->Stride();
      uint maskstride = mask->Stride();
      for (uint i=0; i<itsWindowSize; ++i) {
        rfiStats.Add (getInfo().getAnt1()[bl], getInfo().getAnt2()[bl],
                      itsBuf[i].getTime(), 0, polarization,
                      reals->ValuePtr (i,0), imags->ValuePtr (i,0),
                      mask->ValuePtr (i,0), origFlags->ValuePtr (i,0),
                      nchan, imagestride, maskstride, maskstride);
      }
    }

    void AORFlagger::fillStrategy (boost::shared_ptr<Strategy>& pstrategy)
    {
      string fileName = itsStrategyName;
      if (! fileName.empty()) {
        if (! File(fileName).exists()) {
          fileName = "$LOFARROOT/share/rfistrategies/" + fileName;
          if (! File(fileName).exists()) {
            THROW (Exception, "Unknown rfistrategy file " << itsStrategyName);
          }
        }
        StrategyReader reader;
        pstrategy = boost::shared_ptr<Strategy>
          (reader.CreateStrategyFromFile(fileName));
        return;
      }
      pstrategy = boost::shared_ptr<Strategy> (new Strategy);
      Strategy& strategy = *pstrategy;
      strategy.Add(new SetFlaggingAction());
      ForEachPolarisationBlock* fepBlock = new ForEachPolarisationBlock();
      strategy.Add(fepBlock);
      ActionBlock* current = fepBlock;

      ForEachComplexComponentAction* focAction =
        new ForEachComplexComponentAction();
      focAction->SetOnAmplitude(true);
      focAction->SetOnImaginary(false);
      focAction->SetOnReal(false);
      focAction->SetOnPhase(false);
      focAction->SetRestoreFromAmplitude(false);
      current->Add(focAction);
      current = focAction;

      IterationBlock* iteration = new IterationBlock();
      iteration->SetIterationCount(2);
      iteration->SetSensitivityStart(4.0);
      current->Add(iteration);
      current = iteration;
		
      SumThresholdAction* t2 = new SumThresholdAction();
      t2->SetBaseSensitivity(1.0);
      if (itsPulsarMode) {
        t2->SetFrequencyDirectionFlagging(false);
      }
      current->Add(t2);

      CombineFlagResults* cfr2 = new CombineFlagResults();
      current->Add(cfr2);

      cfr2->Add(new FrequencySelectionAction());
      if (!itsPulsarMode) {
        cfr2->Add(new TimeSelectionAction());
      }
	
      current->Add(new SetImageAction());
      ChangeResolutionAction* changeResAction2 = new ChangeResolutionAction();
      if (itsPulsarMode) {
        changeResAction2->SetTimeDecreaseFactor(1);
      } else {
        changeResAction2->SetTimeDecreaseFactor(3);
      }
      changeResAction2->SetFrequencyDecreaseFactor(3);

      SlidingWindowFitAction* swfAction2 = new SlidingWindowFitAction();
      if (itsPulsarMode) {
        swfAction2->Parameters().timeDirectionWindowSize = 1;
      } else {
        swfAction2->Parameters().timeDirectionKernelSize = 2.5;
        swfAction2->Parameters().timeDirectionWindowSize = 10;
      }
      swfAction2->Parameters().frequencyDirectionKernelSize = 5.0;
      swfAction2->Parameters().frequencyDirectionWindowSize = 15;
      changeResAction2->Add(swfAction2);

      current->Add(changeResAction2);
			
			// This action causes iterations not to converge the thresholds towards the
			// noise, but rather keep using the whole image for threshold calculation.
			// The result is that strongly RFI contaminated sets are very weakly flagged.
			// Commented out on june 11, 2011.
			//current->Add(new SetFlaggingAction());

      current = focAction;
      SumThresholdAction* t3 = new SumThresholdAction();
      if (itsPulsarMode) {
        t3->SetFrequencyDirectionFlagging(false);
      }
      current->Add(t3);
		
      SetFlaggingAction* setFlagsInAllPolarizations = new SetFlaggingAction();
      setFlagsInAllPolarizations->SetNewFlagging
        (SetFlaggingAction::PolarisationsEqual);

      strategy.Add(setFlagsInAllPolarizations);
      strategy.Add(new StatisticalFlagAction());

      if (itsPedantic) {
        CombineFlagResults* cfr3 = new CombineFlagResults();
        strategy.Add(cfr3);
        cfr3->Add(new FrequencySelectionAction());
        if (!itsPulsarMode) {
          cfr3->Add(new TimeSelectionAction());
        }
      } else {
        if (!itsPulsarMode) {
          strategy.Add(new TimeSelectionAction());
        }
      }
    }

  } //# end namespace
}

/*
Hoi Ger,

Ik heb één en ander aan implementatie voor het bijhouden van de quality 
statistics gemaakt, en nu heb ik wat vragen...

Een samenvatting van wat ik gedaan heb: Ik heb twee programma's 
toegevoegd om een beetje met de statistieken te kunnen spelen/testen, in 
de LOFAR repository:

- aoquality
Kan gegeven een MS de statistieken verzamelen en wegschrijven (e.g. 
"aoquality collect SB000.MS" maakt de tabellen), en wat simpele queries 
uitvoeren op de quality tabellen.

- aoqplot
Grafisch programma waarmee je min-of-meer interactief door de 
statistieken kunt browsen, te gebruiken met e.g. "aoqplot SB000.MS".

Beide werk in uitvoering. Het zou mooi zijn als de implementatie 
verbonden kan worden aan NDPPP zodat NDPPP tijdens zijn run de 
statistieken wegschrijft, en de vraag is even hoe en wie dat doet. Het 
is niet zoveel werk meer denk ik. Hoewel ik in principe de statistieken 
in de flagger kan collecten, denk ik dat het netter is als dit buiten de 
flagger gebeurt, omdat er anders een interface moet komen vanuit de 
flagger om de statistieken achteraf weg te schrijven, en dat wordt wat 
clumsy.

Aangezien jij het beste weet waar wat zit in NDPPP, wil jij dit 
toevoegen aan NDPPP? Mijn suggestie is om het volgende te doen in NDPPP :

Benodigde headers:

----
#include <AOFlagger/quality/statisticscollection.h>
#include <AOFlagger/quality/qualitytablesformatter.h>
----

Tijdens de initialiatie moet een StatisticsCollection gemaakt worden. 
Deze moet weten hoeveel polarisaties de set heeft, en moet voor iedere 
band (neem aan dat er altijd maar een band in NDPPP is) een 
initialisatie krijgen:

----
StatisticsCollection collection(polarizationCount);
for(unsigned b=0;b<bandCount;++b)
{
   collection.InitializeBand(b, frequencies[b], bands[b].channelCount);
}
----

Parameters van InitializeBand() zijn:
unsigned subBandIndex;
double *frequencies; //array van de centrale frequenties van ieder 
kanaal in de subband
unsigned channelCount;

Dan moet ieder tijdstap aan de collection gegeven worden:

----
for(unsigned p = 0; p < polarizationCount; ++p)
{
   collection.Add(antenna1Index, antenna2Index, time, bandIndex, p, 
samples[p], isRFI[p]);
}
----

Dit komt dus neer op een Add per regel in de casatable per polarisatie.
samples[p] is een const std::vector<std::complex<float> > &, isRFI is 
een const bool *. Hiervoor zal dus de data uit de casa tables in deze 
structuren gekopieerd moeten worden. samples[p] en isRFI[p] moeten exact 
het aantal elementen hebben als de corresponderende channelCount gegeven 
in InitializeBand(). isRFI bevat de uitvoer van de flagger. Let er op 
dat Add() niet thread safe is, dus als dit multithreaded gebeurt moet de 
call gesynchroniseerd worden. Als dit te veel tijd blijkt te kosten 
kunnen er evt. meerdere "StatisticsCollection collection"'s gemaakt 
worden, zodat iedere thread zijn eigen heeft, en kunnen die aan het eind 
gecombineerd worden, maar denk dat dat niet nodig is.

Uiteindelijk moet dan het volgende stukje eenmaal worden uitgevoerd:

----
QualityTablesFormatter qualityData(msFilename);
collection.Save(qualityData);
----

Dit maakt dan de nieuwe tabellen en vult ze. De stukjes code komen bijna 
volledig uit AOFlagger/src/aoquality.cpp, dus daar zou je kunnen kijken 
voor extra info.

In overleg met verschillende andere mensen kwam ik erachter dat het 
misschien ook handig is om hier en daar nog wat extra statistieken toe 
te voegen op een later tijdstip, met name de statistieken van de 
CORRECTED_DATA kolom, en evt. de statistieken van een tweede flagging run.

Groeten,
André

-- 
..................................
. A.R. Offringa                  .
. Office: +31 50 363 4081        .
. Mobile: +31 6 2879 5672        .
.   Room: 181                    .
. Kapteyn Astronomical Institute .
..................................
*/
