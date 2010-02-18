//# FlagCounter.cc: Class to keep counts of nr of flagged points
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

#include <DPPP/FlagCounter.h>

namespace LOFAR {
  namespace DPPP {

    void FlagCounter::init (uint nbaselines, uint nchan, uint ncorr)
    {
      itsBLCounts.resize (nbaselines);
      itsChanCounts.resize (nchan);
      itsCorrCounts.resize (ncorr);
      std::fill (itsBLCounts.begin(), itsBLCounts.end(), 0);
      std::fill (itsChanCounts.begin(),itsChanCounts.end(), 0);
      std::fill (itsCorrCounts.begin(),itsCorrCounts.end(), 0);
    }

  } //# end namespace
}
