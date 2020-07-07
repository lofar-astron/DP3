// GainCal.cc: DPPP step class to CalculateSomething visibilities
// Copyright (C) 2013
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//
// @author Tammo Jan Dijkema

#include "CalculateSomething.h"

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    CalculateSomething::CalculateSomething (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix)
    : itsInput(input)
    {

    }

    CalculateSomething::~CalculateSomething()
    {}

    void CalculateSomething::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
    }

    void CalculateSomething::show (std::ostream& os) const
    {
      os << "CalculateSomething " << itsName << endl;
    }

    void CalculateSomething::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " CalculateSomething " << itsName << endl;
    }

    bool CalculateSomething::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      //CLA
      //itsBuffer.copy (bufin);
      //itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
      //itsInput->fetchWeights(bufin, itsBuffer, itsTimer);
      cout << "This is my first implementation of CalculateSomething" << endl;
      Array<Complex>::const_contiter indIter = bufin.getData().cbegin();
      int icount = 0;
      //while (outnIter != outnIterEnd) {
      while (icount < 10) {
	      cout << icount << "   " << *indIter << endl;
              icount++;
	      indIter++;
      }
      itsTimer.stop();
      //getNextStep()->process(itsBuffer);
      return false;
    }


    void CalculateSomething::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } // end namespace
}
