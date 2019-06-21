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

#include "DPInput.h"
#include "Exceptions.h"

#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Utilities/Copy.h>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    DPInput::~DPInput()
    {}

    std::string DPInput::msName() const
    {
      return String();
    }

    const Cube<bool>& DPInput::fetchFullResFlags (const DPBuffer& bufin,
                                                  DPBuffer& bufout,
                                                  NSTimer& timer,
                                                  bool merge)
    {
      // If already defined in the buffer, return those fullRes flags.
      if (! bufin.getFullResFlags().empty()) {
        return bufin.getFullResFlags();
      }
      // No fullRes flags in buffer, so get them from the input.
      timer.stop();
      bool fnd = getFullResFlags (bufin.getRowNrs(), bufout);
      timer.start();
      Cube<bool>& fullResFlags = bufout.getFullResFlags();
      if (!fnd) {
        // No fullRes flags in input; form them from the flags in the buffer.
        // Only use the XX flags; no averaging done, thus navgtime=1.
        // (If any averaging was done, the flags would be in the buffer).
        IPosition shp(bufin.getFlags().shape());
        if (fullResFlags.shape()[0] != shp[1]  ||
            fullResFlags.shape()[1] != 1  ||
            fullResFlags.shape()[2] != shp[2])
          throw std::runtime_error("Invalid shape of full res flags");
        objcopy (fullResFlags.data(), bufin.getFlags().data(),
                 fullResFlags.size(), 1, shp[0]);    // only copy XX.
        return fullResFlags;
      }
      // There are fullRes flags.
      // If needed, merge them with the buffer's flags.
      if (merge) {
        DPBuffer::mergeFullResFlags (fullResFlags, bufin.getFlags());
      }
      return fullResFlags;
    }

    const Cube<float>& DPInput::fetchWeights (const DPBuffer& bufin,
                                              DPBuffer& bufout,
                                              NSTimer& timer)
    {
      // If already defined in the buffer, return those weights.
      if (! bufin.getWeights().empty()) {
        return bufin.getWeights();
      }
      // No weights in buffer, so get them from the input.
      // It might need the data and flags in the buffer.
      timer.stop();
      getWeights (bufin.getRowNrs(), bufout);
      timer.start();
      return bufout.getWeights();
    }

    const Matrix<double>& DPInput::fetchUVW (const DPBuffer& bufin,
                                             DPBuffer& bufout,
                                             NSTimer& timer)
    {
      // If already defined in the buffer, return those UVW.
      if (! bufin.getUVW().empty()) {
        return bufin.getUVW();
      }
      // No UVW in buffer, so get them from the input.
      timer.stop();
      getUVW (bufin.getRowNrs(), bufin.getTime(), bufout);
      timer.start();
      return bufout.getUVW();
    }

    void DPInput::getUVW (const RefRows&, double, DPBuffer&)
      { throw Exception ("DPInput::getUVW not implemented"); }

    void DPInput::getWeights (const RefRows&, DPBuffer&)
      { throw Exception ("DPInput::getWeights not implemented"); }

    bool DPInput::getFullResFlags (const RefRows&, DPBuffer&)
      { throw Exception ("DPInput::getFullResFlags not implemented"); }

    void DPInput::getModelData (const RefRows&, Cube<Complex>&)
      { throw Exception ("DPInput::getModelData not implemented"); }

#ifdef HAVE_LOFAR_BEAM
    void DPInput::fillBeamInfo (vector<LOFAR::StationResponse::Station::Ptr>&,
                                const Vector<casacore::String>&)
      { throw Exception ("DPInput::fillBeamInfo not implemented"); }
#endif

  } //# end namespace
}
