//# MSUpdater.cc: DPPP step updating an MS
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
#include <DPPP/MSUpdater.h>
#include <DPPP/MSReader.h>
#include <DPPP/MSWriter.h>
#include <DPPP/DPBuffer.h>
#include <Common/ParameterSet.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSUpdater::MSUpdater (MSReader* reader, const ParameterSet& parset,
                          const string& prefix)
      : itsReader      (reader),
        itsNrCorr      (reader->ncorr()),
        itsNrChan      (reader->nchan()),
        itsNrBl        (reader->nbaselines()),
        itsNrTimes     (0),
        itsFlagCounter ("MSUpdater")
    {
      itsCountFlags = parset.getBool (prefix+"countflag", false);
      MSWriter::writeHistory (reader->table(), parset);
      if (itsCountFlags) {
        itsFlagCounter.init (itsNrBl, itsNrChan, itsNrCorr);
      }
    }

    MSUpdater::~MSUpdater()
    {}

    bool MSUpdater::process (const DPBuffer& buf)
    {
      itsReader->putFlags (buf.getRowNrs(), buf.getFlags());
      // Count the flags if needed.
      if (itsCountFlags) {
        const bool* flagPtr = buf.getFlags().data();
        for (uint i=0; i<itsNrBl; ++i) {
          for (uint j=0; j<itsNrChan; ++j) {
            if (*flagPtr) {
              itsFlagCounter.incrBaseline(i);
              itsFlagCounter.incrChannel(j);
            }
            flagPtr += itsNrCorr;    // only count 1st corr
          }
        }
      }
      itsNrTimes++;
      return true;
    }

    void MSUpdater::finish()
    {}

    void MSUpdater::show (std::ostream& os) const
    {
      os << "MSUpdater" << std::endl;
      os << "  MS:             " << itsReader->msName() << std::endl;
    }

    void MSUpdater::showCounts (std::ostream& os) const
    {
      if (itsCountFlags) {
        os << endl << "Flag statistics of data updated";
        os << endl << "===============================" << endl;
        itsFlagCounter.showBaseline (os, itsReader->getAnt1(),
                                     itsReader->getAnt2(), itsNrTimes);
        itsFlagCounter.showChannel  (os, itsNrTimes);
      }
    }

  } //# end namespace
}
