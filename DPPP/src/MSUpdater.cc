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

#include <tables/Tables/Table.h>
#include <tables/Tables/ArrayColumn.h>
#include <iostream>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    MSUpdater::MSUpdater (MSReader* reader, const ParameterSet& parset,
                          const string& prefix, int needWrite)
      : itsReader      (reader),
        itsWriteData   ((needWrite & DPInfo::NeedWriteData) != 0),
        itsNrCorr      (reader->getInfo().ncorr()),
        itsNrChan      (reader->getInfo().nchan()),
        itsNrBl        (reader->getInfo().nbaselines()),
        itsNrDone      (0)
    {
      itsNrTimesFlush = parset.getUint (prefix+"flush", 0);
      NSTimer::StartStop sstime(itsTimer);
      MSWriter::writeHistory (reader->table(), parset);
    }

    MSUpdater::~MSUpdater()
    {}

    bool MSUpdater::process (const DPBuffer& buf)
    {
      NSTimer::StartStop sstime(itsTimer);
      putFlags (buf.getRowNrs(), buf.getFlags());
      if (itsWriteData) {
        putData (buf.getRowNrs(), buf.getData());
      }
      itsNrDone++;
      if (itsNrTimesFlush > 0  &&  itsNrDone%itsNrTimesFlush == 0) {
        itsReader->table().flush();
      }
      return true;
    }

    void MSUpdater::finish()
    {}

    void MSUpdater::show (std::ostream& os) const
    {
      os << "MSUpdater" << std::endl;
      os << "  MS:             " << itsReader->msName() << std::endl;
    }

    void MSUpdater::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " MSUpdater" << endl;
    }

    void MSUpdater::putFlags (const RefRows& rowNrs,
                              const Cube<bool>& flags)
    {
      if (! rowNrs.rowVector().empty()) {
        const Slicer& colSlicer = itsReader->colSlicer();
        Table& ms = itsReader->table();
        ms.reopenRW();
        ArrayColumn<bool> flagCol(ms, "FLAG");
        ScalarColumn<bool> flagRowCol(ms, "FLAG_ROW");
        // Loop over all rows of this subset.
	// (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<bool> flagIter (flags, 2);
        for (uint i=0; i<rows.size(); ++i) {
          flagCol.putSlice (rows[i], colSlicer, flagIter.array());
          // If a new flag in a row is clear, the ROW_FLAG should not be set.
          // If all new flags are set, we leave it because we might have a
          // subset of the channels, so other flags might still be clear.
          if (anyEQ (flagIter.array(), False)) {
            flagRowCol.put (rows[i], False);
          }
          flagIter.next();
	}
      }
    }

    void MSUpdater::putData (const RefRows& rowNrs,
                             const Cube<Complex>& data)
    {
      if (! rowNrs.rowVector().empty()) {
        const Slicer& colSlicer = itsReader->colSlicer();
        Table& ms = itsReader->table();
        ms.reopenRW();
        ArrayColumn<Complex> dataCol(ms, itsReader->dataColumnName());
        // Loop over all rows of this subset.
	// (it also avoids StandardStMan putCol with RefRows problem).
        Vector<uint> rows = rowNrs.convert();
        ReadOnlyArrayIterator<Complex> dataIter (data, 2);
        for (uint i=0; i<rows.size(); ++i) {
          dataCol.putSlice (rows[i], colSlicer, dataIter.array());
          dataIter.next();
	}
      }
    }

  } //# end namespace
}
