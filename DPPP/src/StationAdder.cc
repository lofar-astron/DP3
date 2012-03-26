//# StationAdder.cc: DPPP step class to add station to a superstation
//# Copyright (C) 2012
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
#include <DPPP/StationAdder.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/ParameterRecord.h>
#include <Common/LofarLogger.h>

#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasConvert.h>
#include <casa/Utilities/LinearSearch.h>
#include <iostream>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    StationAdder::StationAdder (DPInput* input,
                                const ParSet& parset, const string& prefix)
      : itsInput       (input),
        itsName        (prefix),
        itsStatRec     (parset.getRecord(prefix+"stations")),
        itsMinNStation (parset.getUint  (prefix+"minstations", 1))
    {
      // Check the superstation definition(s).
      const Vector<String>& antennaNames = itsInput->antennaNames();
      itsStations.resize (antennaNames.size());
      std::fill (itsStations.begin(), itsStations.end(), -1);
      for (ParameterRecord::const_iterator iter = itsStatRec.begin();
           iter != itsStatRec.end(); ++iter) {
        if (std::find(antennaNames.begin(), antennaNames.end(),
                      String(iter->first)) != antennaNames.end()  ||
            std::find(itsNewNames.begin(),
                      itsNewNames.end(), iter->first) != itsNewNames.end()) {
          THROW (Exception, "StationAdder: new station name " + iter->first +
                            " already exists");
        }
        vector<string> stations = iter->second.getStringVector();
        ASSERT (!stations.empty());
        MVPosition newPosition;
        for (uint i=0; i<stations.size(); ++i) {
          int inx = linearSearch1 (antennaNames, String(stations[i]));
          ASSERTSTR (inx>=0, "Station " + stations[i] + " does not exist");
          ASSERTSTR (itsStations[inx] < 0, "Station " + stations[i] +
                     " is used multiple times ");
          itsStations[inx] = itsNewNames.size();
          newPosition += MPosition::Convert (itsInput->antennaPos()[inx],
                                             MPosition::ITRF)().getValue();
        }
        itsNewNames.push_back (iter->first);
        newPosition *= 1./itsStations.size();
      }
    }

    StationAdder::~StationAdder()
    {}

    void StationAdder::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
    }

    void StationAdder::show (std::ostream& os) const
    {
      os << "StationAdder " << itsName << std::endl;
      os << "  stations:       " << itsStatRec << std::endl;
      os << "  minstations:    " << itsMinNStation << std::endl;
    }

    void StationAdder::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " StationAdder " << itsName << endl;
    }

    bool StationAdder::process (const DPBuffer& buf)
    {
      itsTimer.start();
      RefRows rowNrs(buf.getRowNrs());
      // Sum the data applying the weights.
      itsBuf.getUVW()     = itsInput->fetchUVW (buf, rowNrs, itsTimer);
      itsBuf.getWeights() = itsInput->fetchWeights (buf, rowNrs, itsTimer);
      Cube<bool> fullResFlags(itsInput->fetchFullResFlags (buf, rowNrs,
                                                           itsTimer));
      itsBuf.getData()    = buf.getData();
      IPosition shapeIn   = buf.getData().shape();
      itsBuf.setTime (buf.getTime());
      itsTimer.stop();
      getNextStep()->process (buf);
      return true;
    }

    void StationAdder::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

  } //# end namespace
}
