//# PatchInfo.h: Info about a patch
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
//# $Id: PatchInfo.h 27648 2013-12-04 09:01:16Z diepen $

// @file
// @brief Info about a patch
// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_PATCHINFO_H
#define LOFAR_PARMDB_PATCHINFO_H

#include <cmath>
#include <string>

namespace DP3 {

//# Forward Declarations
class BlobOStream;
class BlobIStream;

namespace BBS {

  // @ingroup ParmDB
  // @{

  // @brief Info about a patch
  class PatchInfo
  {
  public:
    // Default constructor.
    PatchInfo()
    {}

    // Create from patch name, category, ra, dec, and apparent brightness (Jy).
    // Ra and dec must be in radians in J2000.
    PatchInfo (const std::string& name, double ra, double dec, int category,
               double apparentBrightness)
      : itsName          (name),
        itsRa            (ra),
        itsDec           (dec),
        itsCategory      (category),
        itsAppBrightness (apparentBrightness)
    {}

    // Get the patch name.
    const std::string& getName() const
      { return itsName; }

    // Get the right ascension in radians (J2000).
    double getRa() const
      { return itsRa; }

    // Get the declination in radians (J2000).
    double getDec() const
      { return itsDec; }

    // Get the category.
    int getCategory() const
      { return itsCategory; }

    // Get the apparent brightness of the patch (in Jy).
    double apparentBrightness() const
      { return itsAppBrightness; }

    // Set the right ascension in radians (J2000).
    void setRa (double ra)
      { itsRa = ra; }

     // Set the declination in radians (J2000).
    void setDec (double dec)
      { itsDec = dec; }

    // Set the apparent brightness of the patch (in Jy).
    void setApparentBrightness (double apparentBrightness)
      { itsAppBrightness = apparentBrightness; }

  private:
    std::string itsName;
    double itsRa;
    double itsDec;
    int    itsCategory;
    double itsAppBrightness;
  };


  // Show the contents of a PatchInfo object.
  std::ostream& operator<< (std::ostream& os, const PatchInfo& info);

  // Write the contents of a PatchInfo object into a blob.
  BlobOStream operator<< (BlobOStream& os, const PatchInfo& info);

  // Read the contents of a PatchInfo object from a blob.
  BlobIStream operator>> (BlobIStream& os, PatchInfo& info);



  // @brief Info about a patch direction
  class PatchSumInfo
  {
  public:
    // Create from patch name, category, ra, dec, and apparent brightness (Jy).
    // Ra and dec must be in radians in J2000.
    explicit PatchSumInfo (unsigned int patchId)
      : itsSumX    (0),
        itsSumY    (0),
        itsSumZ    (0),
        itsSumFlux (0),
        itsPatchId (patchId)
    {}

    // Add a source direction to determine the average patch direction.
    void add (double ra, double dec, double flux);

    // Get the total flux of the patch.
    double getFlux() const
      { return itsSumFlux; }

    // Get the patch direction (flux-weighted average direction of its sources).
    double getRa() const
      { return std::atan2(itsSumY/itsSumFlux, itsSumX/itsSumFlux); }
    double getDec() const
      { return std::asin(itsSumZ/itsSumFlux); }

    // Get the patchId.
    unsigned int getPatchId() const
      { return itsPatchId; }

  private:
    double itsSumX;
    double itsSumY;
    double itsSumZ;
    double itsSumFlux;
    unsigned int   itsPatchId;
  };


  // @}

} // namespace BBS
} // namespace LOFAR

#endif
