//# DPStep.cc: Abstract base class for a DPPP step
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
#include <DPPP/DPStep.h>

namespace LOFAR {
  namespace DPPP {

    DPStep::~DPStep()
    {}

    void DPStep::updateAverageInfo (AverageInfo&)
    {}


    NullStep::~NullStep()
    {}

    bool NullStep::process (const DPBuffer&)
      { return true; }

    void NullStep::finish()
    {}

    void NullStep::show (std::ostream&)
    {}

  } //# end namespace
}
