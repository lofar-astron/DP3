// AddNoise.h: DPPP step class to add HBA/LBA random noise to data
// Copyright (C) 2013
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief DPPP step class to add HBA/LBA random noise to data
/// @author Tammo Jan Dijkema
//

#ifndef DPPP_AddNoise_H
#define DPPP_AddNoise_H

#include "DPBuffer.h"
#include "DPInput.h"

#include <utility>

#define POL_DEGREE 5

namespace DP3 {

class ParameterSet;

namespace DPPP {
/// @brief DPPP step class to AddNoise visibilities from a source model

class AddNoise : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  AddNoise(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~AddNoise();

  /// Process the data calculating or adding to DATA the estimated random
  /// noise (van Haarlem et al. 2013). When processed, it invokes the process
  /// function of the next step.
  virtual bool process(const DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

 private:
  DPInput* itsInput;
  string itsName;
  DPBuffer itsBuffer;

  NSTimer itsTimer;
  /// The mode parameter select the kind of output from the function
  /// SET mode = 0 data are modified with the noise > data = data + noise
  /// ADD mode = 1 a new array is created > newdata = data + noise
  int mode;
  /// The antennaSet parameter selects between: LBA=LBA_ALL, LBA_INNER
  ///                                        LBA_OUTER and HBA (read from MS)
  string antennaSet;
  /// Coefficients for polynomial interpolation (from constant -first- to highet
  /// order -last-)
  double coeffs_all_LBA[POL_DEGREE + 1] = {
      4.194759691372669e+05,  -0.040624470842582,    1.657038099833047e-09,
      -3.318591685264548e-17, 3.220530883859981e-25, -1.199723767939448e-33};
  double coeffs_inner_LBA[POL_DEGREE + 1] = {
      6.468156199342838e+05,  -0.065296541139271,    2.752773309538937e-09,
      -5.659747549065881e-17, 5.612743945799180e-25, -2.133417264334753e-33};
  double coeffs_outer_LBA[POL_DEGREE + 1] = {
      4.716452746313004e+05,  -0.042301679603444,    1.632924288392071e-09,
      -3.129891572910005e-17, 2.925471926740372e-25, -1.048279929083482e-33};
  double coeffs_cs_HBA[POL_DEGREE + 1] = {
      1.403173283860732e+06,  -0.044309811184301,    5.564568767538230e-10,
      -3.461658816150442e-18, 1.064770242682354e-26, -1.292413137320037e-35};
  double coeffs_rs_HBA[POL_DEGREE + 1] = {
      2.643667639489899e+05,  -0.007554769559170,    8.554267734878638e-11,
      -4.722267097362212e-19, 1.251625317661630e-27, -1.236074872829482e-36};
  /// system efficiency: roughly 1.0
  double eta = 0.95;
};

}  // namespace DPPP
}  // namespace DP3

#endif
