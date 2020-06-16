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

#ifndef SCREEN_CONSTRAINT_H
#define SCREEN_CONSTRAINT_H

#include "../DPPP/PhaseFitter.h"

#include "MultiDirSolver.h"
#include "PiercePoint.h"
#include "KLFitter.h"

#include "../Common/ParameterSet.h"

#include <cmath>
#include <vector>
#include <memory>
#include <ostream>

namespace DP3 {
class ParameterSet;
class ScreenConstraint : public Constraint
{ 
  static const  double phtoTEC;//=1./8.4479745e9;
  static const  double TECtoph;//=8.4479745e9;
  static const  size_t maxIter;//number of iterations to store in debug mode

public:
  ScreenConstraint(const ParameterSet& parset,
                   const string& prefix);

  /** Initialize metadata with frequencies, resize some members.
   * Should be called after InitializeDimensions.
   */
  void initialize(const double* frequencies);
  virtual std::vector<Constraint::Result> Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions,
                   double time, std::ostream* statStream);
  virtual void CalculatePiercepoints();

  void setAntennaPositions(const std::vector<std::vector<double> > antenna_pos);
  void setDirections(const std::vector<std::pair<double, double> > source_pos);
  void setTime(double time);
  void initPiercePoints();
  void setCoreAntennas(const std::vector<size_t>& coreAntennas)
  {
    _coreAntennas = coreAntennas;
    if (itsMode == "csfull")
       _screenFitters.resize(_nAntennas-_coreAntennas.size()+1);
  }
 void setOtherAntennas(const std::vector<size_t>& otherAntennas)
  {
    _otherAntennas = otherAntennas;
  }
 void getPPValue(std::vector<std::vector<std::complex<double> > >&, size_t, size_t, double&,double&) const;
private:
  std::vector<std::vector<double> > itsAntennaPos;
  std::vector<std::vector<double> > itsSourcePos;
  std::vector<double>               itsFrequencies;
  std::vector<double>               itsprevsol;
  std::vector<double>               _iterphases;
  /// antenna positions
  /// source positions
  /// measures instance ofzo               
  std::vector<std::vector<PiercePoint> > itsPiercePoints;  //temporary hold calculated piercepoints per antenna
  std::vector<KLFitter> _screenFitters;
  std::vector<size_t> _coreAntennas;
  std::vector<size_t> _otherAntennas; //has to be a vector for openmp looping
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
}
#endif
