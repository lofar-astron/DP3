// Copyright (C) 2020
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

#ifndef KLFITTER_H
#define KLFITTER_H

#include <vector>

#include <armadillo>

#include "PiercePoint.h"

namespace DP3 {

/// \brief Creates KH base and fits screens from collection of PiercePoints.
class KLFitter {
 public:
  KLFitter(double r0 = 1000., double beta = 5. / 3., int order = 3);
  void calculateCorrMatrix(const std::vector<PiercePoint> pp);
  void calculateCorrMatrix(const std::vector<PiercePoint*> pp);
  void doFit();
  size_t getOrder() const { return itsOrder; }
  double* PhaseData() { return _phases.memptr(); }
  double* WData() { return _weights.memptr(); }
  double* ParData() { return itsPar.memptr(); }
  double* TECFitWhiteData() { return itsTECFitWhite.memptr(); }
  double* PPData() { return itsPiercePoints.memptr(); }
  void setR0(double r0) { itsR0 = r0; }
  void setBeta(double beta) { itsBeta = beta; }
  void setOrder(double order) { itsOrder = order; }
  size_t getNumberofPP() { return itsPiercePoints.n_rows; }

 private:
  size_t itsOrder;
  double itsR0, itsBeta;
  arma::Mat<double> itsPiercePoints;
  arma::Col<double> _phases;
  arma::Mat<double> _weights;       ///< Weights of the data points.
  arma::Mat<double> itsCorrMatrix;  ///< Correlation Matrix for KL fit.
  arma::Mat<double> itsinvC;        ///< For quick interpolation.
  arma::Mat<double> itsU;
  arma::Mat<double> itsinvU;
  arma::Mat<double> itsTECFitWhite;
  arma::Col<double> itsPar;
};
}  // namespace DP3
#endif
