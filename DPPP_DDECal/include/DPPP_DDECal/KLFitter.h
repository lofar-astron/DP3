#ifndef KLFITTER_H
#define KLFITTER_H
#include<vector>
#include <armadillo>
#include <DPPP_DDECal/PiercePoint.h>

using namespace arma;
namespace LOFAR{
class KLFitter
{//creates KH base and fits screens from collection of PiercePoints
public:
  KLFitter(double r0=1000.,double beta=5./3.,int order=3);
  void calculateCorrMatrix(const vector<PiercePoint> pp);
  void calculateCorrMatrix(const vector<PiercePoint*> pp);
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
  Mat<double>             itsPiercePoints;
  Col<double>             _phases;
  Mat<double>             _weights; //weights of the data points
  Mat<double>             itsCorrMatrix; //Correlation Matrix for KL fit
  Mat<double>             itsinvC;  //for quick interpolation
  Mat<double>             itsU;
  Mat<double>             itsinvU;
  Mat<double>             itsTECFitWhite;
  Col<double>             itsPar;
};
}
#endif
