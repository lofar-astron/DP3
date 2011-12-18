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
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>
#include <iostream>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    Demixer::Demixer (DPInput* input,
                      const ParSet& parset, const string& prefix)
      : itsInput     (input),
        itsName      (prefix),
        itsSources   (parset.getStringVector (prefix+"sources")),
        itsNChanAvg  (parset.getUint  (prefix+"freqstep", 1)),
        itsNTimeAvg  (parset.getUint  (prefix+"timestep", 1)),
        itsDMChanAvg (parset.getUint  (prefix+"demixfreqstep", itsNChanAvg)),
        itsDMTimeAvg (parset.getUint  (prefix+"demixtimestep", itsNTimeAvg)),
        itsTimeWindow(parset.getUint  (prefix+"timewindow", 1))
        ///itsJointSolve(parset.getBool  (prefix+"jointsolve", true)),
        ///itsSources   (parset.getStringVector(prefix+"demixsources")),
        ///itsExtra     (parset.getStringVector(prefix+"extrasources"))
    {
      if (itsSources.empty()) {
        // No sources means that nothing has to be demixed, only averaged.
        itsAverager = DPStep::ShPtr (new Averager(input, prefix,
                                                  itsNChanAvg, itsNTimeAvg));
        itsAvgResult = new ResultStep();
        itsAverager->setNextStep (DPStep::ShPtr(itsAvgResult));
        return;
      }
      // Create the steps for the sources to be removed.
      // Demixing consists of the following steps:
      // - phaseshift data to each demix source
      // - average data in each direction, also for original phasecenter.
      // - determine demix factors for all directions
      // - use BBS to predict and solve in each direction. It is possible to
      //   predict more directions than to solve (for strong sources in field).
      // - use BBS to subtract the solved sources using the demix factors.
      //   The averaging used here can be smaller than used when solving.
      itsFirstSteps.reserve     (itsSources.size());
      itsSecondSteps.reserve    (itsSources.size());
      itsDemixInputs.reserve    (itsSources.size());
      itsSubtractInputs.reserve (itsSources.size());
      for (uint i=0; i<itsSources.size(); ++i) {
        ///DPStep::ShPtr step1 (new PhaseShift(input, parset, prefix));
        DPStep::ShPtr step1 (new ResultStep());
        itsFirstSteps.push_back (step1);
        DPStep::ShPtr step2 (new Averager(input, parset, prefix));
        step1->setNextStep (step2);
        ResultStep* step3 = new ResultStep();
        step2->setNextStep (DPStep::ShPtr(step3));
        // There is a single demix step needing all sources.
        itsDemixInputs.push_back (step3);
        ///DPStep::ShPtr step4 (new BBSStep());
        DPStep::ShPtr step4 (new ResultStep());
        itsSecondSteps.push_back (step4);
        ///DPStep::ShPtr step5 (new ParmSmooth());
        /// Smoothing requires a time window; may not be necessary if
        /// BBS is handling outliers correctly.
        DPStep::ShPtr step5 (new ResultStep());
        step4->setNextStep (step5);
        ///DPStep::ShPtr step6 (new BBSStep());
        DPStep::ShPtr step6 (new ResultStep());
        step5->setNextStep (step6);
        ResultStep* step7 = new ResultStep();
        step6->setNextStep (DPStep::ShPtr(step7));
        // There is a single subtract step needing all sources.
        itsSubtractInputs.push_back (step7);
      }
    }

    Demixer::~Demixer()
    {}

    void Demixer::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
      // Let the internal steps update their data.
      // Use a copy of the DPInfo, otherwise it is updated multiple times.
      DPInfo infocp(info);
      if (itsSources.empty()) {
        itsAverager->updateInfo (infocp);
      }
      for (uint i=0; i<itsSources.size(); ++i) {
        infocp = info;
        DPStep::ShPtr step = itsFirstSteps[i];
        while (step) {
          step->updateInfo (infocp);
          step = step->getNextStep();
        }
        step = itsSecondSteps[i];
        while (step) {
          step->updateInfo (infocp);
          step = step->getNextStep();
        }
      }
      // Update the info of this object.
      info.update (itsNChanAvg, itsNTimeAvg);
    }

    void Demixer::show (std::ostream& os) const
    {
      os << "Demixer " << itsName << std::endl;
      os << "  freqstep:       " << itsNChanAvg << std::endl;
      os << "  timestep:       " << itsNTimeAvg << std::endl;
      os << "  demixfreqstep:  " << itsDMChanAvg << std::endl;
      os << "  demixtimestep:  " << itsDMTimeAvg << std::endl;
    }

    void Demixer::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " Demixer " << itsName << endl;
    }

    bool Demixer::process (const DPBuffer& buf)
    {
      itsTimer.start();
      if (itsSources.empty()) {
        return itsAverager->process (buf);
      }
      RefRows rowNrs(buf.getRowNrs());
      for (uint i=0; i<itsSources.size(); ++i) {
        itsFirstSteps[i]->process (buf);
      }
      // If the first steps have produced results, demix and process it.
      if (! itsDemixInputs[0]->get().getData().empty()) {
        DPBuffer buf1 = demix();
        for (uint i=0; i<itsSources.size(); ++i) {
          itsSecondSteps[i]->process (buf1);
          itsDemixInputs[i]->clear();
        }
        DPBuffer buf2 = subtract();
        getNextStep()->process (buf2);
      }
      return true;
    }

    void Demixer::finish()
    {
      // Process remaining entries.
      // Let the next steps finish.
      if (! itsDemixInputs[0]->get().getData().empty()) {
        DPBuffer buf1 = demix();
        for (uint i=0; i<itsSecondSteps.size(); ++i) {
          itsSecondSteps[i]->process (buf1);
        }
        DPBuffer buf2 = subtract();
        getNextStep()->process (buf2);
      }
      getNextStep()->finish();
    }

    DPBuffer Demixer::demix() const
    {
      DPBuffer buf;
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
      /*
      // Calculate matrix for field center.
   x = numpy.sin(ra)*numpy.cos(dec)
   y = numpy.cos(ra)*numpy.cos(dec)
   z = numpy.sin(dec)
   w = numpy.array([[x,y,z]]).T
   x = -numpy.sin(ra)*numpy.sin(dec)
   y = -numpy.cos(ra)*numpy.sin(dec)
   z = numpy.cos(dec)      return itsAvgResult->get();
   v = numpy.array([[x,y,z]]).T
   x = numpy.cos(ra)
   y = -numpy.sin(ra)
   z = 0
   u = numpy.array([[x,y,z]]).T
   T = numpy.concatenate([u,v,w], axis = -1 )
     // Calculate vector for each demix center.
      x1 = numpy.sin(ra1)*numpy.cos(dec1)
      y1 = numpy.cos(ra1)*numpy.cos(dec1)
      z1 = numpy.sin(dec1)
      w1 = numpy.array([[x1,y1,z1]]).T
     // Create demix array.
     Array<DComplex> a (ndir, ndir, npol, nchan, nbl);
     Array<DComplex> ainv (ndir, ndir, npol, nchan, nbl);
     // Read data from all mix directions.
     Array<Complex> data(ndir, npol, nchan, nbl);
     for (i=0; i<ndir; ++i) {
       // Get average weight per 
       weight = sums(buf.getWeight());
       */
      return buf;
    }

    DPBuffer Demixer::subtract() const
    {
      DPBuffer buf;
      return itsAvgResult->get();
      ///return buf;
    }


  } //# end namespace
}
