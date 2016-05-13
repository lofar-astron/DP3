//# GainCal.h: DPPP step class to calibrate (direction independent) gains
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
//# $Id: GainCal.h 21598 2012-07-16 08:07:34Z diepen $
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
      StefCal(uint solInt, uint nChan, string mode, uint maxAntennas);
      void doStep_polarized();
      void doStep_unpolarized(bool phaseOnly);
      Status relax(uint iter);
      void resetVis(uint nSt);
      void init();
      casa::Matrix<casa::DComplex> getSolution();

      bool converged;
      std::vector<int> _antMap; // Length antennaNames, contains size(antennaNames)-nSt times the value -1
                                // Values are indices in the stefcal internal numbering
      casa::Matrix<casa::DComplex> g; // Station, correlation
      casa::Array<casa::DComplex> vis;
      casa::Array<casa::DComplex> mvis;
      uint savedNCr; // number of correlations stored (1,2 or 4)

    private:
      casa::Matrix<casa::DComplex> gx;
      casa::Matrix<casa::DComplex> gxx;
      casa::Matrix<casa::DComplex> gold;
      casa::Matrix<casa::DComplex> h; // Station, correlation
      casa::Matrix<casa::DComplex> z;

      uint nSt; // number of stations
      uint nCr; // number of correlations (1 or 4)
      uint nUn; // number of unknowns
      uint nSp; // number that is two for scalarphase, one else
      uint _solInt;
      uint _nChan;
      string _mode;

      double dg, dgx;
      std::vector<double> dgs;
    };

  } //# end namespace
}

#endif
