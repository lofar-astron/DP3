// SourceData.cc: Class for a Blob file holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SourceData.h"
#include "ParmMap.h"

#include "../blob/BlobIStream.h"
#include "../blob/BlobOStream.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Quanta/MVAngle.h>

#include <iostream>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace parmdb {

SourceData::SourceData() : itsInfo(string(), SourceInfo::POINT) {}

SourceData::SourceData(const SourceInfo& info, const string& patchName,
                       double ra, double dec)
    : itsInfo(info), itsPatchName(patchName), itsRa(ra), itsDec(dec) {}

void SourceData::setParm(const ParmMap& parms, const string& name,
                         double defValue, double& value) {
  // Try to find the parameter with its name and suffixed with source name.
  ParmMap::const_iterator iter = parms.find(name);
  if (iter == parms.end()) {
    iter = parms.find(name + ':' + itsInfo.getName());
  }
  if (iter != parms.end()) {
    const casacore::Array<double>& arr =
        iter->second.getFirstParmValue().getValues();
    if (arr.size() != 1)
      throw std::runtime_error("Error: value " + name + " of source " +
                               itsInfo.getName() + " has multiple values");
    value = arr.data()[0];
  } else {
    value = defValue;
  }
}

void SourceData::setParms(const ParmMap& parms) {
  setParm(parms, "Ra", itsRa, itsRa);
  setParm(parms, "Dec", itsDec, itsDec);
  setParm(parms, "I", 0, itsI);
  setParm(parms, "Q", 0, itsQ);
  setParm(parms, "U", 0, itsU);
  setParm(parms, "V", 0, itsV);
  setParm(parms, "MajorAxis", 0, itsMajorAxis);
  setParm(parms, "MinorAxis", 0, itsMinorAxis);
  setParm(parms, "Orientation", 0, itsOrientation);
  setParm(parms, "PolarizationAngle", 0, itsPolAngle);
  setParm(parms, "PolarizedFraction", 0, itsPolFrac);
  setParm(parms, "RotationMeasure", 0, itsRM);
  itsSpTerms.resize(itsInfo.getNSpectralTerms());
  for (unsigned int i = 0; i < itsSpTerms.size(); ++i) {
    ostringstream ostr;
    ostr << "SpectralIndex:" << i;
    setParm(parms, ostr.str(), 0, itsSpTerms[i]);
  }
}

void SourceData::makeParm(ParmMap& parms, const string& name, double value,
                          bool pertRel) const {
  ParmValueSet pvs(ParmValue(value), ParmValue::Scalar, 1e-6, pertRel);
  parms.define(name, pvs);
}

void SourceData::getParms(ParmMap& parms) const {
  makeParm(parms, "Ra", itsRa, true);
  makeParm(parms, "Dec", itsDec, true);
  makeParm(parms, "I", itsI);
  makeParm(parms, "Q", itsQ);
  makeParm(parms, "U", itsU);
  makeParm(parms, "V", itsV);
  if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
    makeParm(parms, "MajorAxis", itsMajorAxis);
    makeParm(parms, "MinorAxis", itsMinorAxis);
    makeParm(parms, "Orientation", itsOrientation);
  }
  if (itsInfo.getUseRotationMeasure()) {
    makeParm(parms, "PolarizationAngle", itsPolAngle);
    makeParm(parms, "PolarizedFraction", itsPolFrac);
    makeParm(parms, "RotationMeasure", itsRM);
  }
  for (unsigned int i = 0; i < itsSpTerms.size(); ++i) {
    ostringstream ostr;
    ostr << "SpectralIndex:" << i;
    makeParm(parms, ostr.str(), itsSpTerms[i]);
  }
}

void SourceData::writeSource(blob::BlobOStream& bos) const {
  bos.putStart("source", 1);
  itsInfo.write(bos);
  bos << itsPatchName << itsRa << itsDec << itsI << itsQ << itsU << itsV;
  if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
    bos << itsMajorAxis << itsMinorAxis << itsOrientation;
  }
  if (itsInfo.getUseRotationMeasure()) {
    bos << itsPolAngle << itsPolFrac << itsRM;
  }
  if (itsInfo.getNSpectralTerms() > 0) {
    bos.put(itsSpTerms);
  }
  bos.putEnd();
}

void SourceData::readSource(blob::BlobIStream& bis) {
  int version = bis.getStart("source");
  if (version != 1) {
    throw std::domain_error("Unsupported source version");
  }
  itsInfo.read(bis);
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
    bis.get(itsSpTerms);
  } else {
    itsSpTerms.clear();
  }
  bis.getEnd();
}

void SourceData::print(ostream& os) const {
  os << "  ";
  MVAngle(itsRa).print(os, MVAngle::Format(MVAngle::TIME, 9));
  os << ' ';
  MVAngle(itsDec).print(os, MVAngle::Format(MVAngle::ANGLE, 9));
  os << ' ' << itsInfo.getRefType();
  os << "  " << itsInfo.getName() << ' ' << itsInfo.getType();
  os << "  iquv=(" << itsI << ',' << itsQ << ',' << itsU << ',' << itsV << ')';
  os << endl;
  if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
    os << "    major=" << itsMajorAxis << " arcsec  minor=" << itsMinorAxis
       << " arcsec  orientation=" << itsOrientation << " deg" << endl;
  }
  if (itsInfo.getNSpectralTerms() > 0) {
    os << "    nspinx=" << itsInfo.getNSpectralTerms() << " logSI=" << boolalpha
       << itsInfo.getHasLogarithmicSI()
       << " reffreq=" << itsInfo.getSpectralTermsRefFreq() / 1e6 << " MHz"
       << endl;
  }
  if (itsInfo.getUseRotationMeasure()) {
    os << "    polangle=" << itsPolAngle << "  polfrac=" << itsPolFrac
       << "  rm=" << itsRM << endl;
  }
  if (itsInfo.getType() == SourceInfo::SHAPELET) {
    os << "    shapelet I " << itsInfo.getShapeletScaleI()
       << itsInfo.getShapeletCoeffI() << "             Q "
       << itsInfo.getShapeletScaleQ() << itsInfo.getShapeletCoeffQ()
       << "             U " << itsInfo.getShapeletScaleU()
       << itsInfo.getShapeletCoeffU() << "             V "
       << itsInfo.getShapeletScaleV() << itsInfo.getShapeletCoeffV();
  }
}

static ostream& operator<<(std::ostream& output, SourceInfo::Type type) {
  switch (type) {
    case SourceInfo::POINT:
      return output << "POINT";
    case SourceInfo::GAUSSIAN:
      return output << "GAUSSIAN";
    default:
      return output;
  }
}

void toSkymodel(std::ostream& output, const SourceData& source) {
  const SourceInfo& info = source.getInfo();
  output << info.getName() << ", " << info.getType() << ", "
         << source.getPatchName() << ", ";
  MVAngle(source.getRa()).print(output, MVAngle::Format(MVAngle::TIME, 9));
  output << ", ";
  MVAngle(source.getDec()).print(output, MVAngle::Format(MVAngle::ANGLE, 9));
  output << ", " << source.getI() << ", " << info.getSpectralTermsRefFreq()
         << ", [";
  const std::vector<double> spectral_terms = source.getSpectralTerms();
  for (size_t i = 0, e = spectral_terms.size(); i != e; ++i) {
    if (i) output << ", ";
    output << spectral_terms[i];
  }
  if (info.getType() == SourceInfo::GAUSSIAN)
    output << "], " << source.getMajorAxis() << ", " << source.getMinorAxis()
           << ", " << source.getOrientation() << '\n';
  else
    output << "]\n";
}

}  // namespace parmdb
}  // namespace dp3
