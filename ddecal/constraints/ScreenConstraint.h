// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCREEN_CONSTRAINT_H
#define SCREEN_CONSTRAINT_H

#include "../../base/PhaseFitter.h"

#include "Constraint.h"
#include "PiercePoint.h"
#include "KLFitter.h"

#include "../../common/ParameterSet.h"

#include <cmath>
#include <complex>
#include <memory>
#include <ostream>
#include <vector>

namespace dp3 {
namespace common {
class ParameterSet;
}

class ScreenConstraint : public Constraint {
  using DComplex = std::complex<double>;
  static const double phtoTEC;  //=1./8.4479745e9;
  static const double TECtoph;  //=8.4479745e9;
  static const size_t maxIter;  // number of iterations to store in debug mode

 public:
  ScreenConstraint(const common::ParameterSet& parset, const string& prefix);

  /** Initialize metadata with frequencies, resize some members.
   * Should be called after InitializeDimensions.
   */
  void initialize(const double* frequencies);
  virtual std::vector<Constraint::Result> Apply(
      std::vector<std::vector<DComplex> >& solutions, double time,
      std::ostream* statStream);
  virtual void CalculatePiercepoints();

  void setAntennaPositions(
      const std::vector<std::array<double, 3> > antenna_pos);
  void setDirections(const std::vector<std::pair<double, double> > source_pos);
  void setTime(double time);
  void initPiercePoints();
  void setCoreAntennas(const std::vector<size_t>& coreAntennas) {
    _coreAntennas = coreAntennas;
    if (itsMode == "csfull")
      _screenFitters.resize(_nAntennas - _coreAntennas.size() + 1);
  }
  void setOtherAntennas(const std::vector<size_t>& otherAntennas) {
    _otherAntennas = otherAntennas;
  }
  void getPPValue(std::vector<std::vector<std::complex<double> > >&, size_t,
                  size_t, double&, double&) const;

 private:
  std::vector<std::array<double, 3> > itsAntennaPos;
  std::vector<std::vector<double> > itsSourcePos;
  std::vector<double> itsFrequencies;
  std::vector<double> itsprevsol;
  std::vector<double> _iterphases;
  /// antenna positions
  /// source positions
  /// measures instance ofzo
  std::vector<std::vector<PiercePoint> >
      itsPiercePoints;  // temporary hold calculated piercepoints per antenna
  std::vector<KLFitter> _screenFitters;
  std::vector<size_t> _coreAntennas;
  std::vector<size_t> _otherAntennas;  // has to be a vector for openmp looping
  double itsCurrentTime;
  double itsBeta;
  double itsHeight;
  double itsOrder;
  double itsRdiff;
  string itsMode;
  string itsAVGMode;
  int itsDebugMode;
  size_t itsIter;
};
}  // namespace dp3
#endif
