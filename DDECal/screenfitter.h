// screenfitter.h: Class to perform screen fitting
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file screenfitter.h Implements TEC model screen filter @ref ScreenFitter.
 * @author Maaijke Mevius
 * @date 2017-02-01
 */

#ifndef SCREEN_FITTER_H
#define SCREEN_FITTER_H
#include <armadillo>
#include <vector>
using namespace arma;

/// \brief Class to perform screen fitting
class ScreenFitter {
 public:
  ScreenFitter();
  double* PhaseData() { return _phases.data(); }

 private:
  std::vector<double> _phases, _frequencies, _weights;
  mat _corrmatrix;  ///< correlation matrix
};
#endif
