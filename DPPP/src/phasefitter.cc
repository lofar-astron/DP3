#include <lofar_config.h>
#include "DPPP/phasefitter.h"

#include <limits>

double PhaseFitter::TECModelCost(double alpha, double beta) const
{
  double costVal = 0.0;
  for(int i=0; i!=Size(); ++i) {
    double estphase = TECModelFuncWrapped(_frequencies[i], alpha, beta);
    double dCost = fmod(std::fabs(estphase - _phases[i]), 2.0*M_PI);
    if(dCost > M_PI) dCost = 2.0*M_PI - dCost;
		dCost *= _weights[i];
    costVal += dCost;
  }
  return costVal;
}

double PhaseFitter::fitTECModelBeta(double alpha, double betaEstimate) const {
	double weightSum = 0.0;
  for(size_t iter=0; iter!=3; ++iter) {
    double sum = 0.0;
    for(size_t i=0; i!=Size(); ++i) {
      double p = _phases[i], e = TECModelFunc(_frequencies[i], alpha, betaEstimate);
      double dist = fmod(p - e, 2.0*M_PI);
      if(dist < -M_PI)
				dist += 2.0*M_PI;
      else if(dist > M_PI)
				dist -= 2.0*M_PI;
      sum += dist * _weights[i];
			weightSum += _weights[i];
    }
    betaEstimate = betaEstimate + sum / weightSum;
  }
  return fmod(betaEstimate, 2.0*M_PI);
}

void PhaseFitter::bruteForceSearchTECModel(double& lowerAlpha, double& upperAlpha, double& beta) const
{
  double minCost = std::numeric_limits<double>::max();
  double alphaOversampling = 128;
  //size_t betaOversampling = 16;
  double dAlpha = upperAlpha - lowerAlpha;
  const double bandw = _frequencies.back() - _frequencies.front();
  int alphaIndex = 0;
  for(int i=0; i!=alphaOversampling; ++i) {
    // make r between [0, 1]
    double r = double(i)/alphaOversampling;
    double alpha = lowerAlpha + r*dAlpha;
    double curBeta = fitTECModelBeta(alpha, beta);
    double costVal = TECModelCost(alpha, curBeta);
    if(costVal < minCost) {
      beta = curBeta;
      minCost = costVal;
      alphaIndex = i;
    }
  }
  double newLowerAlpha = double(alphaIndex-1)/alphaOversampling*dAlpha + lowerAlpha;
  upperAlpha = double(alphaIndex+1)/alphaOversampling*dAlpha + lowerAlpha;
  lowerAlpha = newLowerAlpha;
}

double PhaseFitter::ternarySearchTECModelAlpha(double startAlpha, double endAlpha, double& beta) const
{
  size_t iter = 0;
  double dCost, lAlpha, rAlpha;
  do {
    lAlpha = startAlpha + (endAlpha - startAlpha) * (1.0/3.0);
    rAlpha = startAlpha + (endAlpha - startAlpha) * (2.0/3.0);
    double lBeta = fitTECModelBeta(lAlpha, beta);
    double rBeta = fitTECModelBeta(rAlpha, beta);
    double lCost = TECModelCost(lAlpha, lBeta);
    double rCost = TECModelCost(rAlpha, rBeta);
    if(lCost < rCost) {
      endAlpha = rAlpha;
      beta = lBeta;
    } else {
      startAlpha = lAlpha;
      beta = rBeta;
    }
    dCost = std::fabs(lCost - rCost);
    ++iter;
  } while(dCost > _fittingAccuracy && iter < 100);
  double finalAlpha = (lAlpha + rAlpha) * 0.5;
  beta = fitTECModelBeta(finalAlpha, beta);
  return finalAlpha;
}

void PhaseFitter::fillDataWithTECModel(double alpha, double beta)
{
  for(size_t ch=0; ch!=Size(); ++ch)
    _phases[ch] = TECModelFunc(_frequencies[ch], alpha, beta);
}

void PhaseFitter::FitTECModelParameters(double& alpha, double& beta) const
{
  double lowerAlpha = 0.0, upperAlpha = 40000.0e6;
  bruteForceSearchTECModel(lowerAlpha, upperAlpha, beta);
  alpha = (lowerAlpha + upperAlpha) * 0.5;
  //beta = fitBeta(alpha, beta);
  alpha = ternarySearchTECModelAlpha(lowerAlpha, upperAlpha, beta);
}

void PhaseFitter::FitDataToTECModel(double& alpha, double& beta)
{
  FitTECModelParameters(alpha, beta);
  fillDataWithTECModel(alpha, beta);
}

