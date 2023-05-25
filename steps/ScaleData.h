// ScaleData.h: DPPP step class for freq-dependent scaling of the data
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step class for freq-dependent scaling of the data
/// @author Ger van Diepen

#ifndef DPPP_SCALEDATA_H
#define DPPP_SCALEDATA_H

#include "InputStep.h"

#include <xtensor/xtensor.hpp>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/BDABuffer.h>

#include <casacore/casa/Arrays/Cube.h>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DPPP step class for freq-dependent scaling of the data

/// This class is a Step class scaling the data using a polynomial
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

class ScaleData : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  ScaleData(const common::ParameterSet&, const string& prefix,
            MsType input_type);

  ~ScaleData() override;

  common::Fields getRequiredFields() const override { return kDataField; }

  common::Fields getProvidedFields() const override { return kDataField; }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer>) override;

  /// Process the DBA data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::BDABuffer>) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  bool accepts(MsType dt) const override { return dt == itsInputType; }

  MsType outputs() const override { return itsInputType; }

 private:
  /// Fill the scale factors for stations having different nr of tiles.
  void fillSizeScaleFactors(unsigned int nNominal, std::vector<float>& fact);

  std::string itsName;
  const MsType itsInputType;
  bool itsScaleSizeGiven;
  bool itsScaleSize;
  std::vector<std::string> itsStationExp;  ///< station regex strings
  std::vector<std::string> itsCoeffStr;    ///< coeff per station regex
  std::vector<std::vector<float>>
      itsStationFactors;             ///< scale factor per station,freq
  xt::xtensor<float, 3> itsFactors;  ///< scale factor per baseline,freq,pol
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
