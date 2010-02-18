//# MSUpdater.h: DPPP step writing to an MS
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

#ifndef DPPP_MSUPDATER_H
#define DPPP_MSUPDATER_H

// @file
// @brief DPPP step writing to an MS

#include <DPPP/DPStep.h>

namespace LOFAR {
  class ParameterSet;

  namespace DPPP {

    //# Forward Declarations.
    class MSReader;

    // @ingroup DPPP

    class MSUpdater: public DPStep
    {
    public:
      MSUpdater (MSReader*, const ParameterSet& parset,
                 const std::string& prefix);

      virtual ~MSUpdater();

      // Process the next data chunk.
      // It returns false when at the end.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Show the step parameters.
      virtual void show (std::ostream&);

    private:
      MSReader* itsReader;
    };

  } //# end namespace
}

#endif
