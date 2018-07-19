//# phasefitter.cc: Class to perform phase fitting (TEC), allowing phase wraps
//# Copyright (C) 2016
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: phasefitter.cc 21598 2012-07-16 08:07:34Z offringa $
//#
//# @author Andre Offringa

#ifdef AOPROJECT
#include "PhaseFitter.h"
#else
#include <lofar_config.h>
#include <DPPP/PhaseFitter.h>
#endif

#include <limits>

double PhaseFitter::TEC2ModelCost(double alpha, double beta) const
{
  double costVal = 0.0, weightSum = 0.0;
  for(size_t i=0; i!=Size(); ++i) {
    double estphase = TEC2ModelFuncWrapped(_frequencies[i], alpha, beta);
    double dCost = fmod(std::fabs(estphase - _phases[i]), 2.0*M_PI);
    if(dCost > M_PI)
      dCost = 2.0*M_PI - dCost;
    dCost *= _weights[i];
    costVal += dCost;
    weightSum += _weights[i];
  }
  if(weightSum == 0.0)
    return 0.0;
  else
    return costVal / weightSum;
}

double PhaseFitter::fitTEC2ModelBeta(double alpha, double betaEstimate) const {
  double weightSum = 0.0;
  for(size_t iter=0; iter!=3; ++iter) {
    double sum = 0.0;
    for(size_t i=0; i!=Size(); ++i) {
      double p = _phases[i], e = TEC2ModelFunc(_frequencies[i], alpha, betaEstimate);
      double dist = fmod(p - e, 2.0*M_PI);
      if(dist < -M_PI)
        dist += 2.0*M_PI;
      else if(dist > M_PI)
        dist -= 2.0*M_PI;
      sum += dist * _weights[i];
      weightSum += _weights[i];
    }
    if(weightSum != 0.0)
      betaEstimate = betaEstimate + sum / weightSum;
  }
  return fmod(betaEstimate, 2.0*M_PI);
}

void PhaseFitter::bruteForceSearchTEC2Model(double& lowerAlpha, double& upperAlpha, double& beta) const
{
  double minCost = std::numeric_limits<double>::max();
  double alphaOversampling = 256;
  //size_t betaOversampling = 16;
  double dAlpha = upperAlpha - lowerAlpha;
  int alphaIndex = 0;
  for(int i=0; i!=alphaOversampling; ++i) {
    // make r between [0, 1]
    double r = double(i)/alphaOversampling;
    double alpha = lowerAlpha + r*dAlpha;
    double curBeta = fitTEC2ModelBeta(alpha, beta);
    double costVal = TEC2ModelCost(alpha, curBeta);
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

double PhaseFitter::ternarySearchTEC2ModelAlpha(double startAlpha, double endAlpha, double& beta) const
{
  size_t iter = 0;
  double dCost, lAlpha, rAlpha;
  do {
    lAlpha = startAlpha + (endAlpha - startAlpha) * (1.0/3.0);
    rAlpha = startAlpha + (endAlpha - startAlpha) * (2.0/3.0);
    double lBeta = fitTEC2ModelBeta(lAlpha, beta);
    double rBeta = fitTEC2ModelBeta(rAlpha, beta);
    double lCost = TEC2ModelCost(lAlpha, lBeta);
    double rCost = TEC2ModelCost(rAlpha, rBeta);
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
  beta = fitTEC2ModelBeta(finalAlpha, beta);
  return finalAlpha;
}

void PhaseFitter::fillDataWithTEC2Model(double alpha, double beta)
{
  for(size_t ch=0; ch!=Size(); ++ch)
    _phases[ch] = TEC2ModelFunc(_frequencies[ch], alpha, beta);
}

void PhaseFitter::fillDataWithTEC1Model(double alpha)
{
  for(size_t ch=0; ch!=Size(); ++ch)
    _phases[ch] = TEC1ModelFuncWrapped(_frequencies[ch], alpha);
}

void PhaseFitter::FitTEC2ModelParameters(double& alpha, double& beta) const
{
  double lowerAlpha = -40000.0e6, upperAlpha = 40000.0e6;
  bruteForceSearchTEC2Model(lowerAlpha, upperAlpha, beta);
  alpha = (lowerAlpha + upperAlpha) * 0.5;
  //beta = fitBeta(alpha, beta);
  alpha = ternarySearchTEC2ModelAlpha(lowerAlpha, upperAlpha, beta);
}

double PhaseFitter::FitDataToTEC2Model(double& alpha, double& beta)
{
  FitTEC2ModelParameters(alpha, beta);
  double cost = TEC2ModelCost(alpha, beta);
  fillDataWithTEC2Model(alpha, beta);
  return cost;
}

double PhaseFitter::FitDataToTEC1Model(double& alpha)
{
  FitTEC1ModelParameters(alpha);
  double cost = TEC1ModelCost(alpha);
  fillDataWithTEC1Model(alpha);
  return cost;
}

void PhaseFitter::FitTEC1ModelParameters(double& alpha) const
{
  double lowerAlpha = -40000.0e6, upperAlpha = 40000.0e6;
  bruteForceSearchTEC1Model(lowerAlpha, upperAlpha);
  alpha = ternarySearchTEC1ModelAlpha(lowerAlpha, upperAlpha);
}

#include <iostream>
void PhaseFitter::bruteForceSearchTEC1Model(double& lowerAlpha, double& upperAlpha) const
{
  double minCost = std::numeric_limits<double>::max();
  double alphaOversampling = 256;
  double dAlpha = upperAlpha - lowerAlpha;
  int alphaIndex = 0;
  for(int i=0; i!=alphaOversampling; ++i) {
    // make r between [0, 1]
    double r = double(i)/alphaOversampling;
    double alpha = lowerAlpha + r*dAlpha;
    // We have to have some freedom in the fit to make sure
    // we do rule out an area with an unwripping that is correct
    // Hence we use the two-parameter model and allow beta to be fitted.
    // The ternary search will fix alpha to accomodate a zero beta.
    double curBeta = fitTEC2ModelBeta(alpha, 0.0);
    double costVal = TEC2ModelCost(alpha, curBeta);
    if(costVal < minCost) {
      minCost = costVal;
      alphaIndex = i;
    }
  }
  double newLowerAlpha = double(alphaIndex-1)/alphaOversampling*dAlpha + lowerAlpha;
  upperAlpha = double(alphaIndex+1)/alphaOversampling*dAlpha + lowerAlpha;
  lowerAlpha = newLowerAlpha;
  //std::cout << "alpha in " << lowerAlpha << "-" << upperAlpha << '\n';
}

double PhaseFitter::TEC1ModelCost(double alpha) const
{
  double costVal = 0.0, weightSum = 0.0;
  for(size_t i=0; i!=Size(); ++i) {
    double estphase = TEC1ModelFuncWrapped(_frequencies[i], alpha);
    double dCost = fmod(std::fabs(estphase - _phases[i]), 2.0*M_PI);
    if(dCost > M_PI)
      dCost = 2.0*M_PI - dCost;
    dCost *= _weights[i];
    costVal += dCost;
    weightSum += _weights[i];
  }
  if(weightSum == 0.0)
    return 0.0;
  else
    return costVal / weightSum;
}

double PhaseFitter::ternarySearchTEC1ModelAlpha(double startAlpha, double endAlpha) const
{
  size_t iter = 0;
  double dCost, lAlpha, rAlpha;
  do {
    lAlpha = startAlpha + (endAlpha - startAlpha) * (1.0/3.0);
    rAlpha = startAlpha + (endAlpha - startAlpha) * (2.0/3.0);
    double lCost = TEC1ModelCost(lAlpha);
    double rCost = TEC1ModelCost(rAlpha);
    if(lCost < rCost) {
      endAlpha = rAlpha;
    } else {
      startAlpha = lAlpha;
    }
    dCost = std::fabs(lCost - rCost);
    ++iter;
    //std::cout << iter << '\t' << startAlpha << '\t' << endAlpha << '\n';
  } while(dCost > _fittingAccuracy && iter < 100);
  double finalAlpha = (lAlpha + rAlpha) * 0.5;
  return finalAlpha;
}

