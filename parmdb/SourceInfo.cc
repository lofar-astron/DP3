// SourceInfo.cc: Info about a source
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// @file
// @brief Info about a source
// @author Ger van Diepen (diepen AT astron nl)

#include "SourceInfo.h"

#include "../blob/BlobIStream.h"
#include "../blob/BlobOStream.h"
#include "../blob/BlobArray.h"

#include <casacore/casa/Arrays/Array.h>

namespace dp3 {
namespace parmdb {

SourceInfo::SourceInfo(const string& name, Type type, const string& refType,
                       bool useLogarithmicSI, unsigned int nSpectralTerms,
                       double spectralTermsRefFreqHz, bool useRotationMeasure)
    : itsName(name),
      itsType(type),
      itsRefType(refType),
      itsNSpTerms(nSpectralTerms),
      itsSpTermsRefFreq(spectralTermsRefFreqHz),
      itsHasLogarithmicSI(useLogarithmicSI),
      itsUseRotMeas(useRotationMeasure),
      itsShapeletScaleI(0),
      itsShapeletScaleQ(0),
      itsShapeletScaleU(0),
      itsShapeletScaleV(0) {}

SourceInfo::SourceInfo(const SourceInfo& that) { operator=(that); }

SourceInfo& SourceInfo::operator=(const SourceInfo& that) {
  if (this != &that) {
    itsName = that.itsName;
    itsType = that.itsType;
    itsRefType = that.itsRefType;
    itsNSpTerms = that.itsNSpTerms;
    itsSpTermsRefFreq = that.itsSpTermsRefFreq;
    itsHasLogarithmicSI = that.itsHasLogarithmicSI;
    itsUseRotMeas = that.itsUseRotMeas;
    itsShapeletScaleI = that.itsShapeletScaleI;
    itsShapeletScaleQ = that.itsShapeletScaleQ;
    itsShapeletScaleU = that.itsShapeletScaleU;
    itsShapeletScaleV = that.itsShapeletScaleV;
    itsShapeletCoeffI.assign(that.itsShapeletCoeffI);
    itsShapeletCoeffQ.assign(that.itsShapeletCoeffQ);
    itsShapeletCoeffU.assign(that.itsShapeletCoeffU);
    itsShapeletCoeffV.assign(that.itsShapeletCoeffV);
  }
  return *this;
}

void SourceInfo::setShapeletCoeff(const casacore::Array<double>& I,
                                  const casacore::Array<double>& Q,
                                  const casacore::Array<double>& U,
                                  const casacore::Array<double>& V) {
  itsShapeletCoeffI.assign(I);
  itsShapeletCoeffQ.assign(Q);
  itsShapeletCoeffU.assign(U);
  itsShapeletCoeffV.assign(V);
}

void SourceInfo::setShapeletScale(double scaleI, double scaleQ, double scaleU,
                                  double scaleV) {
  itsShapeletScaleI = scaleI;
  itsShapeletScaleQ = scaleQ;
  itsShapeletScaleU = scaleU;
  itsShapeletScaleV = scaleV;
}

void SourceInfo::write(blob::BlobOStream& bos) const {
  int16_t version = 2;
  // TODO itsHasLogarithmicSI has to be written in new blob version
  bos << version << itsName << int16_t(itsType) << itsRefType
      << itsHasLogarithmicSI << itsNSpTerms << itsSpTermsRefFreq
      << itsUseRotMeas;
  if (itsType == SHAPELET) {
    bos << itsShapeletScaleI << itsShapeletScaleQ << itsShapeletScaleU
        << itsShapeletScaleV << itsShapeletCoeffI << itsShapeletCoeffQ
        << itsShapeletCoeffU << itsShapeletCoeffV;
  }
}

// If ever version info is needed,
void SourceInfo::read(blob::BlobIStream& bis) {
  int16_t version, type;
  bis >> version >> itsName >> type >> itsRefType;
  if (version != 1 && version != 2)
    throw std::runtime_error("Version of sourcedb must be 1 or 2");
  if (version >= 2) {
    bis >> itsHasLogarithmicSI;
  } else {
    itsHasLogarithmicSI = true;
  }
  bis >> itsNSpTerms >> itsSpTermsRefFreq >> itsUseRotMeas;
  // Convert to enum.
  itsType = Type(type);
  if (itsType == SHAPELET) {
    bis >> itsShapeletScaleI >> itsShapeletScaleQ >> itsShapeletScaleU >>
        itsShapeletScaleV >> itsShapeletCoeffI >> itsShapeletCoeffQ >>
        itsShapeletCoeffU >> itsShapeletCoeffV;
  } else {
    itsShapeletScaleI = 0;
    itsShapeletScaleQ = 0;
    itsShapeletScaleU = 0;
    itsShapeletScaleV = 0;
    itsShapeletCoeffI.resize();
    itsShapeletCoeffQ.resize();
    itsShapeletCoeffU.resize();
    itsShapeletCoeffV.resize();
  }
}

void toSkymodel(std::ostream& output, const SourceInfo& source);
}  // namespace parmdb
}  // namespace dp3
