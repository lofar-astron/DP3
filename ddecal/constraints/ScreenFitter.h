// screenfitter.h: Class to perform screen fitting
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file ScreenFitter.h Implements TEC model screen filter @ref ScreenFitter.
 * @author Maaijke Mevius
 * @date 2017-02-01
 */

#ifndef DP3_DDECAL_SCREEN_FITTER_H_
#define DP3_DDECAL_SCREEN_FITTER_H_

#include <armadillo>
#include <vector>

namespace dp3 {
namespace ddecal {

/// \brief Class to perform screen fitting
class ScreenFitter {
 public:
  ScreenFitter();
  double* PhaseData() { return _phases.data(); }

 private:
  std::vector<double> _phases, _frequencies, _weights;
  arma::mat _corrmatrix;  ///< correlation matrix
};

}  // namespace ddecal
}  // namespace dp3

#endif
