//# Counter.cc: DPPP step class to count flags
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
#include <DPPP/Counter.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    Counter::Counter (DPInput* input,
                      const ParSet& parset, const string& prefix)
      : itsName        (prefix),
        itsCount       (0),
        itsFlagCounter (input->msName(), parset, prefix)
    {
      itsFlagData = parset.getBool (prefix+"flagdata", false);
    }

    Counter::~Counter()
    {}

    void Counter::show (std::ostream& os) const
    {
      os << "Counter " << itsName << std::endl;
    }

    void Counter::showCounts (std::ostream& os) const
    {
      os << endl << "Cumulative flag counts in Counter " << itsName;
      os << endl << "=================================" << endl;
      itsFlagCounter.showBaseline (os, itsCount);
      itsFlagCounter.showChannel  (os, itsCount);
    }

    void Counter::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      // Visibility data must be read if needed, so NaNs are flagged.
      if (itsFlagData) {
        info().setNeedVisData();
      }
      // Initialize the flag counters.
      itsFlagCounter.init (infoIn);
    }

    bool Counter::process (const DPBuffer& buf)
    {
      const IPosition& shape = buf.getFlags().shape();
      uint nrcorr = shape[0];
      uint nrchan = shape[1];
      uint nrbl   = shape[2];
      const bool* flagPtr = buf.getFlags().data();
      for (uint i=0; i<nrbl; ++i) {
        for (uint j=0; j<nrchan; ++j) {
          if (*flagPtr) {
            itsFlagCounter.incrBaseline(i);
            itsFlagCounter.incrChannel(j);
          }
          flagPtr += nrcorr;    // only count 1st corr
        }
      }
      // Let the next step do its processing.
      getNextStep()->process (buf);
      itsCount++;
      return true;
    }

    void Counter::finish()
    {
      // Let the next step finish its processing.
      getNextStep()->finish();
    }

  } //# end namespace
}
