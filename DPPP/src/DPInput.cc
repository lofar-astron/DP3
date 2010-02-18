//# DPInput.cc: Abstract base class for a DPStep generating input
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
#include <DPPP/DPInput.h>
#include <Common/Exception.h>
#include <casa/Utilities/Copy.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    DPInput::~DPInput()
    {}

    Cube<bool> DPInput::fetchPreAvgFlags (const DPBuffer& buf,
                                          const RefRows& rowNrs,
                                          bool merge)
    {
      // If already defined in the buffer, return those preAvg flags.
      if (! buf.getPreAvgFlags().empty()) {
        return buf.getPreAvgFlags();
      }
      // No preAvg flags in buffer, so get them from the input.
      Cube<bool> preAvgFlags (getPreAvgFlags(rowNrs));
      if (preAvgFlags.empty()) {
        // No preAvg flags in input; form them from the flags in the buffer.
        // Only use the XX flags; no averaging done, thus navgtime=1.
        IPosition shp(buf.getFlags().shape());
        IPosition ofShape(shp[1], 1, shp[2]);    // nchan,navgtime,nbl
        preAvgFlags.resize (ofShape);
        objcopy (preAvgFlags.data(), buf.getFlags().data(),
                 preAvgFlags.size(), shp[0], 1);    // only copy XX.
        return preAvgFlags;
      }
      // There are preAvg flags.
      // If needed, merge them with the buffer's flags.
      if (merge) {
        DPBuffer::mergePreAvgFlags (preAvgFlags, buf.getFlags());
      }
      return preAvgFlags;
    }

    Cube<float> DPInput::fetchWeights (const DPBuffer& buf,
                                       const RefRows& rowNrs)
    {
      // If already defined in the buffer, return those weights.
      if (! buf.getWeights().empty()) {
        return buf.getWeights();
      }
      // No weights in buffer, so get them from the input.
      return getWeights(rowNrs);
    }

    Matrix<double> DPInput::fetchUVW (const DPBuffer& buf,
                                      const RefRows& rowNrs)
    {
      // If already defined in the buffer, return those UVW.
      if (! buf.getUVW().empty()) {
        return buf.getUVW();
      }
      // No UVW in buffer, so get them from the input.
      return getUVW(rowNrs);
    }

    Matrix<double> DPInput::getUVW (const RefRows&)
      { throw Exception ("DPInput::getUVW not implemented"); }

    Cube<float> DPInput::getWeights (const RefRows&)
      { throw Exception ("DPInput::getWeights not implemented"); }

    Cube<bool> DPInput::getPreAvgFlags (const RefRows&)
      { throw Exception ("DPInput::getPreAvgFlags not implemented"); }

    Cube<Complex> DPInput::getData (const String&, const RefRows&)
      { throw Exception ("DPInput::getData not implemented"); }

  } //# end namespace
}
