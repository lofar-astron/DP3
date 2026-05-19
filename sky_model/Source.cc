// SourceData.cc: Class for a Blob file holding sources and their parameters
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Source.h"

#include "../parmdb/ParmMap.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Quanta/MVAngle.h>

#include <iostream>

using casacore::MVAngle;
using dp3::parmdb::ParmMap;

namespace dp3::sky_model {

Source::Source() : itsInfo(string(), SourceInfo::POINT) {}

Source::Source(const SourceInfo& info, const std::string& patchName, double ra,
               double dec)
    : itsInfo(info), itsPatchName(patchName), itsRa(ra), itsDec(dec) {}

void Source::setParm(const std::map<std::string, parmdb::ParmValue>& parms,
                     const std::string& name, double defValue, double& value) {
  // Try to find the parameter with its name and suffixed with source name.
  auto iter = parms.find(name);
  if (iter == parms.end()) {
    iter = parms.find(name + ':' + itsInfo.getName());
  }
  if (iter != parms.end()) {
    const casacore::Array<double>& arr = iter->second.getValues();
    if (arr.size() != 1)
      throw std::runtime_error("Error: value " + name + " of source " +
                               itsInfo.getName() + " has multiple values");
    value = arr.data()[0];
  } else {
    value = defValue;
  }
}

void Source::setParms(const std::map<std::string, parmdb::ParmValue>& parms) {
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
    std::ostringstream ostr;
    ostr << "SpectralIndex:" << i;
    setParm(parms, ostr.str(), 0, itsSpTerms[i]);
  }
}

void Source::print(std::ostream& os) const {
  os << "  ";
  MVAngle(itsRa).print(os, MVAngle::Format(MVAngle::TIME, 9));
  os << ' ';
  MVAngle(itsDec).print(os, MVAngle::Format(MVAngle::ANGLE, 9));
  os << ' ' << itsInfo.getRefType();
  os << "  " << itsInfo.getName() << ' ' << itsInfo.getType();
  os << "  iquv=(" << itsI << ',' << itsQ << ',' << itsU << ',' << itsV << ')';
  os << '\n';
  if (itsInfo.getType() == SourceInfo::GAUSSIAN) {
    os << "    major=" << itsMajorAxis << " arcsec  minor=" << itsMinorAxis
       << " arcsec  orientation=" << itsOrientation << " deg";
    if (itsInfo.getPositionAngleIsAbsolute()) {
      os << " (absolute)" << '\n';
    } else {
      os << " (w.r.t. North at phase center)" << '\n';
    }
  }
  if (itsInfo.getNSpectralTerms() > 0) {
    os << "    nspinx=" << itsInfo.getNSpectralTerms()
       << " logSI=" << std::boolalpha << itsInfo.getHasLogarithmicSI()
       << " reffreq=" << itsInfo.getSpectralTermsRefFreq() / 1e6 << " MHz"
       << '\n';
  }
  if (itsInfo.getUseRotationMeasure()) {
    os << "    polangle=" << itsPolAngle << "  polfrac=" << itsPolFrac
       << "  rm=" << itsRM << '\n';
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

static std::ostream& operator<<(std::ostream& output, SourceInfo::Type type) {
  switch (type) {
    case SourceInfo::POINT:
      return output << "POINT";
    case SourceInfo::GAUSSIAN:
      return output << "GAUSSIAN";
    default:
      return output;
  }
}

void toSkyModel(std::ostream& output, const Source& source) {
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
  output << "], ";
  output << std::boolalpha << source.getInfo().getHasLogarithmicSI();
  if (info.getType() == SourceInfo::GAUSSIAN) {
    output << ", " << source.getMajorAxis() << ", " << source.getMinorAxis()
           << ", " << source.getOrientation() << ", ";
    output << std::boolalpha << source.getInfo().getPositionAngleIsAbsolute();
  }
  output << "\n";
}

}  // namespace dp3::sky_model
