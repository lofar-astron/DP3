//# SourceData.cc: Class for a Blob file holding sources and their parameters
//#
//# Copyright (C) 2012
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
//# $Id: SourceData.cc 37340 2017-05-11 12:39:06Z dijkema $

#include "SourceData.h"
#include "ParmMap.h"

#include "../Blob/BlobIStream.h"
#include "../Blob/BlobOStream.h"

#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Quanta/MVAngle.h>

#include <iostream>

using namespace casacore;
using namespace std;

namespace DP3 {
namespace BBS {


  SourceData::SourceData()
    : itsInfo (string(), SourceInfo::POINT)
  {}

  SourceData::SourceData (const SourceInfo& info,
                          const string& patchName,
                          double ra, double dec)
    : itsInfo      (info),
      itsPatchName (patchName),
      itsRa        (ra),
      itsDec       (dec)
  {}

  void SourceData::setParm (const ParmMap& parms, const string& name,
                            double defValue, double& value)
  {
    // Try to find the parameter with its name and suffixed with source name.
    ParmMap::const_iterator iter = parms.find(name);
    if (iter == parms.end()) {
      iter = parms.find(name + ':' + itsInfo.getName());
    }
    if (iter != parms.end()) {
      const casacore::Array<double>& arr =
        iter->second.getFirstParmValue().getValues();
      if (arr.size()!=1)
				throw std::runtime_error("Error: value " + name + " of source " +
                 itsInfo.getName() + " has multiple values");
      value = arr.data()[0];
    } else {
      value = defValue;
    }
  }

  void SourceData::setParms (const ParmMap& parms)
  {
    setParm (parms, "Ra", itsRa, itsRa);
    setParm (parms, "Dec", itsDec, itsDec);
    setParm (parms, "I", 0, itsI);
    setParm (parms, "Q", 0, itsQ);
    setParm (parms, "U", 0, itsU);
    setParm (parms, "V", 0, itsV);
    setParm (parms, "MajorAxis", 0, itsMajorAxis);
    setParm (parms, "MinorAxis", 0, itsMinorAxis);
    setParm (parms, "Orientation", 0, itsOrientation);
    setParm (parms, "PolarizationAngle", 0, itsPolAngle);
    setParm (parms, "PolarizedFraction", 0, itsPolFrac);
    setParm (parms, "RotationMeasure", 0, itsRM);
    itsSpTerms.resize (itsInfo.getNSpectralTerms());
    for (unsigned int i=0; i<itsSpTerms.size(); ++i) {
      ostringstream ostr;
      ostr << "SpectralIndex:" << i;
      setParm (parms, ostr.str(), 0, itsSpTerms[i]);
    }
  }

  void SourceData::makeParm (ParmMap& parms, const string& name,
                             double value, bool pertRel) const
  {
    ParmValueSet pvs(ParmValue(value), ParmValue::Scalar, 1e-6, pertRel);
    parms.define (name, pvs);
  }

  void SourceData::getParms (ParmMap& parms) const
  {
    makeParm (parms, "Ra", itsRa, true);
    makeParm (parms, "Dec", itsDec, true);
    makeParm (parms, "I", itsI);
    makeParm (parms, "Q", itsQ);
    makeParm (parms, "U", itsU);
    makeParm (parms, "V", itsV);
    if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
      makeParm (parms, "MajorAxis", itsMajorAxis);
      makeParm (parms, "MinorAxis", itsMinorAxis);
      makeParm (parms, "Orientation", itsOrientation);
    }
    if (itsInfo.getUseRotationMeasure()) {
      makeParm (parms, "PolarizationAngle", itsPolAngle);
      makeParm (parms, "PolarizedFraction", itsPolFrac);
      makeParm (parms, "RotationMeasure", itsRM);
    }
    for (unsigned int i=0; i<itsSpTerms.size(); ++i) {
      ostringstream ostr;
      ostr << "SpectralIndex:" << i;
      makeParm (parms, ostr.str(), itsSpTerms[i]);
    }
  }

  void SourceData::writeSource (BlobOStream& bos) const
  {
    bos.putStart ("source", 1);
    itsInfo.write (bos);
    bos << itsPatchName << itsRa << itsDec << itsI << itsQ << itsU << itsV;
    if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
      bos << itsMajorAxis << itsMinorAxis << itsOrientation;
    }
    if (itsInfo.getUseRotationMeasure()) {
      bos << itsPolAngle << itsPolFrac << itsRM;
    }
    if (itsInfo.getNSpectralTerms() > 0) {
      bos.put (itsSpTerms);
    }
    bos.putEnd();
  }

  void SourceData::readSource (BlobIStream& bis)
  {
    int version = bis.getStart ("source");
    assert (version == 1);
    itsInfo.read (bis);
    bis >> itsPatchName >> itsRa >> itsDec >> itsI >> itsQ >> itsU >> itsV;
    if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
      bis >> itsMajorAxis >> itsMinorAxis >> itsOrientation;
    } else {
      itsMajorAxis = itsMinorAxis = itsOrientation = 0;
    }
    if (itsInfo.getUseRotationMeasure()) {
      bis >> itsPolAngle >> itsPolFrac >> itsRM;
    } else {
      itsPolAngle = itsPolFrac = itsRM = 0;
    }
    if (itsInfo.getNSpectralTerms() > 0) {
      bis.get (itsSpTerms);
    } else {
      itsSpTerms.clear();
    }
    bis.getEnd();
  }

  void SourceData::print (ostream& os) const
  {
    os << "  ";
    MVAngle(itsRa).print (os, MVAngle::Format(MVAngle::TIME, 9));
    os << ' ';
    MVAngle(itsDec).print (os, MVAngle::Format(MVAngle::ANGLE, 9));
    os << ' ' << itsInfo.getRefType();
    os << "  " << itsInfo.getName() << ' ' << itsInfo.getType();
    os << "  iquv=(" << itsI << ',' << itsQ << ',' << itsU <<',' << itsV << ')';
    os << endl;
    if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
      os << "    major=" << itsMajorAxis << " arcsec  minor=" << itsMinorAxis
         << " arcsec  orientation=" << itsOrientation << " deg" << endl;
    }
    if (itsInfo.getNSpectralTerms() > 0) {
      os << "    nspinx=" << itsInfo.getNSpectralTerms()
         << " logSI=" << boolalpha << itsInfo.getHasLogarithmicSI()
         << " reffreq=" << itsInfo.getSpectralTermsRefFreq()/1e6 << " MHz"
         << endl;
    }
    if (itsInfo.getUseRotationMeasure()) {
      os << "    polangle=" << itsPolAngle << "  polfrac=" << itsPolFrac
         << "  rm=" << itsRM << endl;
    }
    if (itsInfo.getType() == SourceInfo::SHAPELET) {
      os << "    shapelet I " << itsInfo.getShapeletScaleI()
         << itsInfo.getShapeletCoeffI()
         << "             Q " << itsInfo.getShapeletScaleQ()
         << itsInfo.getShapeletCoeffQ()
         << "             U " << itsInfo.getShapeletScaleU()
         << itsInfo.getShapeletCoeffU()
         << "             V " << itsInfo.getShapeletScaleV()
         << itsInfo.getShapeletCoeffV();
    }
  }

} // namespace BBS
} // namespace LOFAR
