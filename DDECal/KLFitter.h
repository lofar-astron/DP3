#ifndef KLFITTER_H
#define KLFITTER_H

#include <vector>

#include <armadillo>

#include "PiercePoint.h"

namespace DP3{
class KLFitter
{//creates KH base and fits screens from collection of PiercePoints
public:
  KLFitter(double r0=1000.,double beta=5./3.,int order=3);
  void calculateCorrMatrix(const std::vector<PiercePoint> pp);
  void calculateCorrMatrix(const std::vector<PiercePoint*> pp);
  void doFit();
  size_t getOrder() const {return itsOrder;}
  double* PhaseData() { return _phases.memptr(); }
  double* WData() { return _weights.memptr(); }
  double* ParData() { return itsPar.memptr(); }
  double* TECFitWhiteData() { return itsTECFitWhite.memptr(); }
  double* PPData() { return itsPiercePoints.memptr(); }
  void setR0(double r0) {itsR0=r0;}
  void setBeta(double beta) {itsBeta=beta;}
  void setOrder(double order) {itsOrder=order;}
  size_t getNumberofPP() {return itsPiercePoints.n_rows;}

private:
  size_t                  itsOrder;
  double                  itsR0,itsBeta;
  arma::Mat<double>             itsPiercePoints;
  arma::Col<double>             _phases;
  arma::Mat<double>             _weights; //weights of the data points
  arma::Mat<double>             itsCorrMatrix; //Correlation Matrix for KL fit
  arma::Mat<double>             itsinvC;  //for quick interpolation
  arma::Mat<double>             itsU;
  arma::Mat<double>             itsinvU;
  arma::Mat<double>             itsTECFitWhite;
  arma::Col<double>             itsPar;
};
}
#endif
