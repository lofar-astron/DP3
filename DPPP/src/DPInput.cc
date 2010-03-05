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

    Vector<double> DPInput::chanFreqs (uint nchanAvg) const
    {
      uint nchan = itsChanFreqs.size() / nchanAvg;
      Vector<double> freqs(nchan); 
      for (uint i=0; i<nchan; ++i) {
        freqs[i] = 0.5 * (itsChanFreqs[i*nchanAvg] +
                          itsChanFreqs[(i+1)*nchanAvg - 1]);
      }
      return freqs;
    }

    Cube<bool> DPInput::fetchFullResFlags (const DPBuffer& buf,
                                           const RefRows& rowNrs,
                                           bool merge)
    {
      // If already defined in the buffer, return those fullRes flags.
      if (! buf.getFullResFlags().empty()) {
        return buf.getFullResFlags();
      }
      // No fullRes flags in buffer, so get them from the input.
      Cube<bool> fullResFlags (getFullResFlags(rowNrs));
      if (fullResFlags.empty()) {
        // No fullRes flags in input; form them from the flags in the buffer.
        // Only use the XX flags; no averaging done, thus navgtime=1.
        IPosition shp(buf.getFlags().shape());
        IPosition ofShape(3, shp[1], 1, shp[2]);    // nchan,navgtime,nbl
        fullResFlags.resize (ofShape);
        objcopy (fullResFlags.data(), buf.getFlags().data(),
                 fullResFlags.size(), 1, shp[0]);    // only copy XX.
        return fullResFlags;
      }
      // There are fullRes flags.
      // If needed, merge them with the buffer's flags.
      if (merge) {
        DPBuffer::mergeFullResFlags (fullResFlags, buf.getFlags());
      }
      return fullResFlags;
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

    Cube<bool> DPInput::getFullResFlags (const RefRows&)
      { throw Exception ("DPInput::getFullResFlags not implemented"); }

    Cube<Complex> DPInput::getData (const String&, const RefRows&)
      { throw Exception ("DPInput::getData not implemented"); }

  } //# end namespace
}
