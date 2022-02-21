// PatchInfo.cc: Info about a patch
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// @file
// @brief Info about a patch
// @author Ger van Diepen (diepen AT astron nl)

#include "PatchInfo.h"

#include "../blob/BlobIStream.h"
#include "../blob/BlobOStream.h"

#include <casacore/casa/Quanta/MVAngle.h>

#include <cassert>

using namespace casacore;

namespace dp3 {
namespace parmdb {

void PatchSumInfo::add(double ra, double dec, double flux) {
  // Add the position in xyz coordinates.
  // Use the flux as the weight.
  double cosDec = cos(dec);
  itsSumX += cos(ra) * cosDec * flux;
  itsSumY += sin(ra) * cosDec * flux;
  itsSumZ += sin(dec) * flux;
  itsSumFlux += flux;
}

std::ostream& operator<<(std::ostream& os, const PatchInfo& info) {
  os << "patch=" << info.getName() << " cat=" << info.getCategory();
  os << " ra=";
  MVAngle(info.getRa()).print(os, MVAngle::Format(MVAngle::TIME, 9));
  os << " dec=";
  MVAngle(info.getDec()).print(os, MVAngle::Format(MVAngle::ANGLE, 9));
  os << " flux=" << info.apparentBrightness();
  return os;
}

// Write the contents of a PatchInfo object into a blob.
blob::BlobOStream operator<<(blob::BlobOStream& bos, const PatchInfo& info) {
  int16_t cType = info.getCategory();
  bos.putStart("patch", 1);
  bos << info.getName() << cType << info.apparentBrightness() << info.getRa()
      << info.getDec();
  bos.putEnd();
  return bos;
}

// Read the contents of a PatchInfo object from a blob.
blob::BlobIStream operator>>(blob::BlobIStream& bis, PatchInfo& info) {
  string patchName;
  int16_t cType;
  double apparentBrightness, ra, dec;
  if (bis.getStart("patch") != 1)  // version must be 1
    throw std::runtime_error("Version of patch must be 1");
  bis >> patchName >> cType >> apparentBrightness >> ra >> dec;
  bis.getEnd();
  info = PatchInfo(patchName, ra, dec, cType, apparentBrightness);
  return bis;
}

void toSkymodel(std::ostream& output, const PatchInfo& patch) {
  output << ", , " << patch.getName() << ", ";
  MVAngle(patch.getRa()).print(output, MVAngle::Format(MVAngle::TIME, 9));
  output << ", ";
  MVAngle(patch.getDec()).print(output, MVAngle::Format(MVAngle::ANGLE, 9));
  output << '\n';
}

}  // namespace parmdb
}  // namespace dp3
