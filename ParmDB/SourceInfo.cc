//# SourceInfo.cc: Info about a source
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
//# $Id: SourceInfo.cc 37340 2017-05-11 12:39:06Z dijkema $

// @file
// @brief Info about a source
// @author Ger van Diepen (diepen AT astron nl)

//# Includes
#include "SourceInfo.h"

#include "../Blob/BlobIStream.h"
#include "../Blob/BlobOStream.h"
#include "../Blob/BlobArray.h"

#include <casacore/casa/Arrays/Array.h>


namespace DP3 {
namespace BBS {

  SourceInfo::SourceInfo (const string& name, Type type,
                          const string& refType,
                          bool useLogarithmicSI,
                          unsigned int nSpectralTerms,
                          double spectralTermsRefFreqHz,
                          bool useRotationMeasure)
    : itsName           (name),
      itsType           (type),
      itsRefType        (refType),
      itsNSpTerms    (nSpectralTerms),
      itsSpTermsRefFreq   (spectralTermsRefFreqHz),
      itsHasLogarithmicSI(useLogarithmicSI),
      itsUseRotMeas     (useRotationMeasure),
      itsShapeletScaleI (0),
      itsShapeletScaleQ (0),
      itsShapeletScaleU (0),
      itsShapeletScaleV (0)
  {}

  SourceInfo::SourceInfo (const SourceInfo& that)
  {
    operator= (that);
  }

  SourceInfo& SourceInfo::operator= (const SourceInfo& that)
  {
    if (this != &that) {
      itsName           = that.itsName;
      itsType           = that.itsType;
      itsRefType        = that.itsRefType;
      itsNSpTerms    = that.itsNSpTerms;
      itsSpTermsRefFreq   = that.itsSpTermsRefFreq;
      itsHasLogarithmicSI = that.itsHasLogarithmicSI;
      itsUseRotMeas     = that.itsUseRotMeas;
      itsShapeletScaleI = that.itsShapeletScaleI;
      itsShapeletScaleQ = that.itsShapeletScaleQ;
      itsShapeletScaleU = that.itsShapeletScaleU;
      itsShapeletScaleV = that.itsShapeletScaleV;
      itsShapeletCoeffI.assign (that.itsShapeletCoeffI);
      itsShapeletCoeffQ.assign (that.itsShapeletCoeffQ);
      itsShapeletCoeffU.assign (that.itsShapeletCoeffU);
      itsShapeletCoeffV.assign (that.itsShapeletCoeffV);
    }
    return *this;
  }

  void SourceInfo::setShapeletCoeff (const casacore::Array<double>& I,
                                     const casacore::Array<double>& Q,
                                     const casacore::Array<double>& U,
                                     const casacore::Array<double>& V)
  {
    itsShapeletCoeffI.assign (I);
    itsShapeletCoeffQ.assign (Q);
    itsShapeletCoeffU.assign (U);
    itsShapeletCoeffV.assign (V);
  }

  void SourceInfo::setShapeletScale (double scaleI, double scaleQ,
                                     double scaleU, double scaleV)
  {
    itsShapeletScaleI = scaleI;
    itsShapeletScaleQ = scaleQ;
    itsShapeletScaleU = scaleU;
    itsShapeletScaleV = scaleV;
  }

  void SourceInfo::write (BlobOStream& bos) const
  {
    int16_t version = 2;
      // TODO itsHasLogarithmicSI has to be written in new blob version
    bos << version << itsName << int16_t(itsType) << itsRefType
        << itsHasLogarithmicSI
        << itsNSpTerms << itsSpTermsRefFreq << itsUseRotMeas;
    if (itsType == SHAPELET) {
      bos << itsShapeletScaleI << itsShapeletScaleQ
          << itsShapeletScaleU << itsShapeletScaleV
          << itsShapeletCoeffI << itsShapeletCoeffQ
          << itsShapeletCoeffU << itsShapeletCoeffV;
    }
  }

  // If ever version info is needed, 
  void SourceInfo::read (BlobIStream& bis)
  {
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
      bis >> itsShapeletScaleI >> itsShapeletScaleQ
          >> itsShapeletScaleU >> itsShapeletScaleV
          >> itsShapeletCoeffI >> itsShapeletCoeffQ
          >> itsShapeletCoeffU >> itsShapeletCoeffV;
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

} // namespace BBS
} // namespace LOFAR
