//# GainCal.cc: DPPP step class to H5ParmPredict visibilities
//# Copyright (C) 2013
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <lofar_config.h>
#include <DPPP/H5ParmPredict.h>

#include <iostream>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <Common/StreamUtil.h>

using namespace casacore;
using namespace LOFAR::BBS;

namespace LOFAR {
  namespace DPPP {

    H5ParmPredict::H5ParmPredict (DPInput* input,
                      const ParameterSet& parset,
                      const string& prefix):
                          itsInput(input),
                          itsH5ParmName(parset.getString(prefix+"h5parm")),
                          itsDirections(parset.getStringVector(
                              prefix+"directions", vector<string> ())),
                          itsSolTabName(parset.getString("correction"))
    {

    }

    H5ParmPredict::~H5ParmPredict()
    {}

    void H5ParmPredict::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();

      itsH5Parm = H5Parm(itsH5ParmName, false);


    }

    void H5ParmPredict::show (std::ostream& os) const
    {
      os << "H5ParmPredict " << itsName << endl;
      os << "  H5Parm:     " << itsH5ParmName;
      os << "  directions: " << itsDirections;
    }

    void H5ParmPredict::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " H5ParmPredict " << itsName << endl;
    }

    bool H5ParmPredict::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      itsInput->fetchUVW(bufin, itsBuffer, itsTimer);
      itsInput->fetchWeights(bufin, itsBuffer, itsTimer);

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }


    void H5ParmPredict::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
  } //# end namespace
}
