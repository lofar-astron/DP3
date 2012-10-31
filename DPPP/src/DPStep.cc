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

    const DPInfo& DPStep::setInfo (const DPInfo& info)
    {
      // Update the info of this step using the given info.
      updateInfo (info);
      // If there is a next step, set its info using the info of this step.
      if (getNextStep()) {
        return getNextStep()->setInfo (getInfo());
      }
      return getInfo();
    }

    void DPStep::updateInfo (const DPInfo& infoIn)
      { info() = infoIn; }

    void DPStep::addToMS (const string&)
    {}

    void DPStep::showCounts (std::ostream&) const
    {}

    void DPStep::showTimings (std::ostream&, double) const
    {}


    NullStep::~NullStep()
    {}

    bool NullStep::process (const DPBuffer&)
      { return true; }

    void NullStep::finish()
    {}

    void NullStep::show (std::ostream&) const
    {}


    ResultStep::ResultStep()
    {
      setNextStep (DPStep::ShPtr (new NullStep()));
    }

    ResultStep::~ResultStep()
    {}

    bool ResultStep::process (const DPBuffer& buf)
    {
      itsBuffer = buf;
      getNextStep()->process (buf);
      return true;
    }

    void ResultStep::finish()
    {
      getNextStep()->finish();
    }

    void ResultStep::show (std::ostream&) const
    {}


    MultiResultStep::MultiResultStep (uint reserveSize)
    {
      setNextStep (DPStep::ShPtr (new NullStep()));
      itsBuffers.reserve (reserveSize);
    }

    MultiResultStep::~MultiResultStep()
    {}

    bool MultiResultStep::process (const DPBuffer& buf)
    {
      itsBuffers.push_back (buf);
      getNextStep()->process (buf);
      return true;
    }

    void MultiResultStep::finish()
    {
      getNextStep()->finish();
    }

    void MultiResultStep::show (std::ostream&) const
    {}

  } //# end namespace
}
