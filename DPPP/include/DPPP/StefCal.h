//# StefCal.h: Perform StefCal algorithm for gain calibration
//# Copyright (C) 2013
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
//# $Id: StefCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_STEFCAL_H
#define DPPP_STEFCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>

namespace LOFAR {

  namespace DPPP {
    // @ingroup NDPPP

    class StefCal
    {
    public:
      enum Status {CONVERGED=1, NOTCONVERGED=2, STALLED=3};

      // mode can be "diagonal", "fulljones", "phaseonly", "scalarphase"
      StefCal(uint solInt, uint nChan, const string& mode, double tolerance,
              uint maxAntennas, bool detectStalling, uint debugLevel);

      // Sets visibility matrices to zero, resizes them to nSt stations
      void resetVis(uint nSt);

      // Initializes a new run of stefcal, resizes all internal vectors
      void init();

      // Perform sone iteration of stefcal. Returns CONVERGED, NOTCONVERGED
      // or STALLED
      Status doStep(uint iter);

      // Returns the solution. The return matrix has a length of maxAntennas,
      // which is zero for antennas for which no solution was computed.
      // The mapping is stored in the antenna map
      casa::Matrix<casa::DComplex> getSolution();

      // Returns a reference to the visibility matrix
      casa::Array<casa::DComplex>& getVis() {
        return _vis;
      }

      // Returns a reference to the model visibility matrix
      casa::Array<casa::DComplex>& getMVis() {
        return _mvis;
      }

      // Returns a reference to the antenna map. This map has length
      // maxAntennas, its values are the indices in the stefcal internal
      // numbering, or zero when an antenna was not solved for
      std::vector<int>& getAntMap() {
        return _antMap;
      }

      // Number of correlations in the solution (1,2 or 4)
      uint numCorrelations() {
        return _savedNCr;
      }

    private:
      // Number of unknowns (nSt or 2*nSt, depending on _mode)
      uint nUn();

      // Perform relaxation
      Status relax(uint iter);

      void doStep_polarized();
      void doStep_unpolarized();

      uint _savedNCr;
      std::vector<int> _antMap; // Length antennaNames, contains size(antennaNames)-nSt times the value -1
      casa::Array<casa::DComplex> _vis; // Visibility matrix
      casa::Array<casa::DComplex> _mvis; // Model visibility matrix
      casa::Matrix<casa::DComplex> _g; // Solution, indexed by station, correlation
      casa::Matrix<casa::DComplex> _gx; // Previous solution
      casa::Matrix<casa::DComplex> _gxx; // Solution before previous solution
      casa::Matrix<casa::DComplex> _gold; // Previous solution
      casa::Matrix<casa::DComplex> _h; // Hermitian transpose of previous solution
      casa::Matrix<casa::DComplex> _z; // Internal stefcal vector

      uint _nSt; // number of stations in the current solution
      uint _nUn; // number of unknowns
      uint _nCr; // number of correlations (1 or 4)
      uint _nSp; // number that is two for scalarphase, one else
      uint _solInt; // solution interval
      uint _nChan;  // number of channels
      string _mode; // diagonal, scalarphase, fulljones or phaseonly
      double _tolerance;
      bool _detectStalling;
      uint _debugLevel;

      double _dg, _dgx; // previous convergence
      std::vector<double> _dgs; // convergence history
    };

  } //# end namespace
}

#endif
