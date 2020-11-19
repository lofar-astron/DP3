// ScaleData.h: DPPP step class for freq-dependent scaling of the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class for freq-dependent scaling of the data
/// @author Ger van Diepen

#ifndef DPPP_SCALEDATA_H
#define DPPP_SCALEDATA_H

#include "DPInput.h"
#include "DPBuffer.h"
#include "BDABuffer.h"

#include <casacore/casa/Arrays/Cube.h>

namespace DP3 {

class ParameterSet;

namespace DPPP {
/// @brief DPPP step class for freq-dependent scaling of the data

/// This class is a DPStep class scaling the data using a polynomial
/// in frequency (in MHz) for LBA and HBA. The coefficients can be given
/// as ParSet parameters having a default determined by Adam Deller.
///
/// The polynomial coefficients can depend on station by giving them
/// per station name regular expression. The default coefficients are
/// used for the station not matching any regular expression.
///
/// The data are multiplied with a factor sqrt(scale[ant1] * scale[ant2]).
/// An extra scale factor can be applied to compensate for the different
/// number of dipoles or tiles or for missing ones. By default that
/// extra scale factor is only applied to stations using the default
/// coefficients, because it is assumed that coefficients are scaled well
/// when specifying them explicitly for stations.

class ScaleData : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ScaleData(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~ScaleData();

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer&);

  /// Process the DBA data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<BDABuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Return which datatype this step outputs.
  MSType outputs() const override;

  /// Boolean if this step can process this type of data.
  bool accepts(MSType dt) const override;

 private:
  /// Fill the scale factors for stations having different nr of tiles.
  void fillSizeScaleFactors(unsigned int nNominal, std::vector<double>& fact);

  string itsName;
  bool itsScaleSizeGiven;
  bool itsScaleSize;
  std::vector<string> itsStationExp;  ///< station regex strings
  std::vector<string> itsCoeffStr;    ///< coeff per station regex
  std::vector<std::vector<double> >
      itsStationFactors;              ///< scale factor per station,freq
  casacore::Cube<double> itsFactors;  ///< scale factor per baseline,freq,pol
  NSTimer itsTimer;
};

}  // namespace DPPP
}  // namespace DP3

#endif
