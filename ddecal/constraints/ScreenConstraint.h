// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCREEN_CONSTRAINT_H
#define SCREEN_CONSTRAINT_H

#include "../../base/PhaseFitter.h"
#include "../../base/Direction.h"

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

namespace ddecal {

class ScreenConstraint : public Constraint {
  using DComplex = std::complex<double>;
  static const double phtoTEC;  //=1./8.4479745e9;
  static const double TECtoph;  //=8.4479745e9;
  static const size_t maxIter;  // number of iterations to store in debug mode

 public:
  ScreenConstraint(const common::ParameterSet& parset, const string& prefix);

  void Initialize(size_t nAntennas, size_t nDirections,
                  const std::vector<double>& frequencies) override;
  std::vector<Constraint::Result> Apply(
      std::vector<std::vector<DComplex>>& solutions, double time,
      std::ostream* statStream) override;

  void SetCoreAntennas(const std::set<size_t>& core_antennas);
  void InitPiercePoints(const std::vector<std::array<double, 3>>& antenna_pos,
                        const std::vector<base::Direction>& source_directions);

  const std::vector<size_t>& GetCoreAntennas() const { return _coreAntennas; }
  const std::vector<std::vector<PiercePoint>>& GetPiercePoints() const {
    return itsPiercePoints;
  }

 private:
  void SetTime(double time);
  void CalculatePiercepoints();
  void GetPpValue(const std::vector<std::vector<std::complex<double>>>&,
                  size_t solution_index, size_t direction_index,
                  double& avg_tec, double& error) const;

  std::vector<double> itsFrequencies;
  std::vector<double> itsprevsol;
  std::vector<double> _iterphases;
  /// antenna positions
  /// source positions
  /// measures instance ofzo
  std::vector<std::vector<PiercePoint>>
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

}  // namespace ddecal
}  // namespace dp3

#endif
