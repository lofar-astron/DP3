// SourceData.h: Class holding all parameters of a source
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Base class for a table holding sources and their parameters
/// @author Ger van Diepen (diepen AT astron nl)

#ifndef LOFAR_PARMDB_SOURCEDATA_H
#define LOFAR_PARMDB_SOURCEDATA_H

#include "SourceInfo.h"

#include <iosfwd>

namespace dp3::parmdb {
class ParmValue;
}

namespace dp3::sky_model {

/// @brief Class holding data of a source.
class Source {
 public:
  Source();

  Source(const SourceInfo&, const std::string& patchName, double ra,
         double dec);

  /// Get the various source parameters.
  const SourceInfo& getInfo() const { return itsInfo; }
  const std::string& getPatchName() const { return itsPatchName; }
  /// Get right ascension in radians.
  double getRa() const { return itsRa; }
  /// Get declination in radians.
  double getDec() const { return itsDec; }
  double getI() const { return itsI; }
  double getQ() const { return itsQ; }
  double getU() const { return itsU; }
  double getV() const { return itsV; }
  /// Get major axis in arcsec.
  double getMajorAxis() const { return itsMajorAxis; }
  /// Get minor axis in arcsec.
  double getMinorAxis() const { return itsMinorAxis; }
  /// Get orientation in degrees.
  double getOrientation() const { return itsOrientation; }
  double getPolarizationAngle() const { return itsPolAngle; }
  double getPolarizedFraction() const { return itsPolFrac; }
  double getRotationMeasure() const { return itsRM; }
  const std::vector<double>& getSpectralTerms() const { return itsSpTerms; }

  /// Set the various source parameters.
  ///@{
  void setInfo(const SourceInfo& info) { itsInfo = info; }
  void setPatchName(const std::string& patchName) { itsPatchName = patchName; }
  void setRa(double ra) { itsRa = ra; }
  void setDec(double dec) { itsDec = dec; }
  void setI(double i) { itsI = i; }
  void setQ(double q) { itsQ = q; }
  void setU(double u) { itsU = u; }
  void setV(double v) { itsV = v; }
  void setMajorAxis(double majorAxis) { itsMajorAxis = majorAxis; }
  void setMinorAxis(double minorAxis) { itsMinorAxis = minorAxis; }
  void setOrientation(double orientation) { itsOrientation = orientation; }
  void setPolarizationAngle(double polarizationAngle) {
    itsPolAngle = polarizationAngle;
  }
  void setPolarizedFraction(double polarizedFraction) {
    itsPolFrac = polarizedFraction;
  }
  void setRotationMeasure(double potationMeasure) { itsRM = potationMeasure; }
  void setSpectralTerms(const std::vector<double>& spectralTerms) {
    itsSpTerms = spectralTerms;
  }
  ///@}

  void setParms(
      const std::map<std::string, parmdb::ParmValue>& defaultParameters);

  /// Print the source data.
  void print(std::ostream&) const;

 private:
  /// Set a parameter.
  /// If defined, its value is taken from the map.
  /// Otherwise the default value is used.
  void setParm(const std::map<std::string, parmdb::ParmValue>& parms,
               const std::string& name, double defValue, double& value);

  SourceInfo itsInfo;
  string itsPatchName;
  double itsRa;   ///< radians
  double itsDec;  ///< radians
  double itsI;
  double itsQ;
  double itsU;
  double itsV;
  double itsMajorAxis;    ///< arcsec
  double itsMinorAxis;    ///< arcsec
  double itsOrientation;  ///< degrees
  double itsPolAngle;
  double itsPolFrac;
  double itsRM;
  std::vector<double> itsSpTerms;
};

/// Output a source to a sky_model text file.
///
/// The output format is used for @code showsourcedb mode=sky_model @endcode
void toSkyModel(std::ostream& output, const Source& source);

}  // namespace dp3::sky_model

#endif
