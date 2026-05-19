// SourceInfo.cc: Info about a source
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// @file
// @brief Info about a source
// @author Ger van Diepen (diepen AT astron nl)

#include "SourceInfo.h"

#include "blob/BlobIStream.h"
#include "blob/BlobOStream.h"
#include "blob/BlobArray.h"

#include <casacore/casa/Arrays/Array.h>

namespace dp3::sky_model {

SourceInfo::SourceInfo(const std::string& name, Type type,
                       const std::string& refType, bool useLogarithmicSI,
                       unsigned int nSpectralTerms,
                       double spectralTermsRefFreqHz, bool useRotationMeasure,
                       bool positionAngleIsAbsolute)
    : itsName(name),
      itsType(type),
      itsRefType(refType),
      itsNSpTerms(nSpectralTerms),
      itsSpTermsRefFreq(spectralTermsRefFreqHz),
      itsHasLogarithmicSI(useLogarithmicSI),
      itsUseRotMeas(useRotationMeasure),
      itsPositionAngleIsAbsolute(positionAngleIsAbsolute),
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
    itsPositionAngleIsAbsolute = that.itsPositionAngleIsAbsolute;
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

}  // namespace dp3::sky_model
