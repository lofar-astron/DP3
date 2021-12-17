// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_KLFITTER_H_
#define DP3_DDECAL_KLFITTER_H_

#include <vector>

#include <armadillo>

#include "PiercePoint.h"

namespace dp3 {
namespace ddecal {

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

}  // namespace ddecal
}  // namespace dp3

#endif
