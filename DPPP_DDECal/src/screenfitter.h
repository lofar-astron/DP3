//# screenfitter.h: Class to perform screen fitting 
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
//#
//# @author Maaijke Mevius

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

class ScreenFitter{
 public:
  ScreenFitter();
  double* PhaseData() { return _phases.data(); }

 private:
  std::vector<double> _phases,_frequencies, _weights;
  mat _corrmatrix;  //correlation matrix
   
};
#endif
