//# AverageInfo.cc: Info how the data are averaged in time and frequency
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
#include <DPPP/AverageInfo.h>
#include <Common/LofarLogger.h>

namespace LOFAR {
  namespace DPPP {

    AverageInfo::AverageInfo()
      : itsStartChan    (0),
        itsOrigNChan    (0),
        itsNChan        (0),
        itsChanAvg      (1),
        itsNTime        (0),
        itsTimeAvg      (1),
        itsTimeInterval (0)
    {}

    void AverageInfo::init (uint startChan, uint nchan, uint ntime,
                            double timeInterval)
    {
      itsStartChan    = startChan;
      itsOrigNChan    = nchan;
      itsNChan        = nchan;
      itsNTime        = ntime;
      itsTimeInterval = timeInterval;
    }

    void AverageInfo::update (uint chanAvg, uint timeAvg)
    {
      ASSERTSTR (itsNChan % chanAvg == 0,
                 "When averaging, nr of channels must divide integrally");
      itsChanAvg *= chanAvg;
      itsNChan = (itsNChan + chanAvg - 1) / chanAvg;
      itsTimeAvg *= timeAvg;
      itsNTime = (itsNTime + timeAvg - 1) / timeAvg;
      itsTimeInterval *= timeAvg;
    }

  } //# end namespace
}
