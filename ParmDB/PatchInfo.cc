//# PatchInfo.cc: Info about a patch
//#
//# Copyright (C) 2008
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
//# $Id: PatchInfo.h 25297 2013-06-12 11:39:35Z diepen $

// @file
// @brief Info about a patch
// @author Ger van Diepen (diepen AT astron nl)


//# Includes
#include "PatchInfo.h"

#include "../Blob/BlobIStream.h"
#include "../Blob/BlobOStream.h"

#include <casacore/casa/Quanta/MVAngle.h>

#include <cassert>

using namespace casacore;

namespace DP3 {
namespace BBS {

  void PatchSumInfo::add (double ra, double dec, double flux)
  {
    // Add the position in xyz coordinates.
    // Use the flux as the weight.
    double cosDec = cos(dec);
    itsSumX += cos(ra) * cosDec * flux;
    itsSumY += sin(ra) * cosDec * flux;
    itsSumZ += sin(dec) * flux;
    itsSumFlux += flux;
  }

  std::ostream& operator<< (std::ostream& os, const PatchInfo& info)
  {
    os << "patch=" << info.getName() << " cat=" << info.getCategory();
    os << " ra=";
    MVAngle(info.getRa()).print(os, MVAngle::Format(MVAngle::TIME, 9));
    os << " dec=";
    MVAngle(info.getDec()).print (os, MVAngle::Format(MVAngle::ANGLE, 9));
    os << " flux=" << info.apparentBrightness();
    return os;
  }

  // Write the contents of a PatchInfo object into a blob.
  BlobOStream operator<< (BlobOStream& bos, const PatchInfo& info)
  {
    int16_t cType = info.getCategory();
    bos.putStart ("patch", 1);
    bos << info.getName() << cType << info.apparentBrightness()
        << info.getRa() << info.getDec();
    bos.putEnd();
    return bos;
  }

  // Read the contents of a PatchInfo object from a blob.
  BlobIStream operator>> (BlobIStream& bis, PatchInfo& info)
  {
    string patchName;
    int16_t cType;
    double apparentBrightness, ra, dec;
    assert (bis.getStart ("patch") == 1);   // version must be 1
    bis >> patchName >> cType >> apparentBrightness >> ra >> dec;
    bis.getEnd();
    info = PatchInfo (patchName, ra, dec, cType, apparentBrightness);
    return bis;
  }


} // namespace BBS
} // namespace LOFAR
