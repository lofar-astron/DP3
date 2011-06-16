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
        itsHalfWindow(parset.getUint  (prefix+"halfwindow", 1)),
        itsThreshold (parset.getDouble(prefix+"threshold", 1.))
    {
      ASSERTSTR (itsSources.size() > 0, "Demixer: no sources given ");
      // Create the step to average the data and to get the result.
      itsAverager = DPStep::ShPtr (new Averager(input, parset, prefix));
      itsAvgResult = new ResultStep();
      itsAverager->setNextStep (DPStep::ShPtr(itsAvgResult));
      // Create the steps for the sources to be removed.
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
      itsAverager->updateInfo (infocp);
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
      RefRows rowNrs(buf.getRowNrs());
      itsAverager->process (buf);
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
      return itsAvgResult->get();
      ///return buf;
    }

    DPBuffer Demixer::subtract() const
    {
      DPBuffer buf;
      return itsAvgResult->get();
      ///return buf;
    }


  } //# end namespace
}
