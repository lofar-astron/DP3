//# PhaseShift.cc: DPPP step class to shift the data to another phase center
//# Copyright (C) 2010
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
//# $Id$
//#
//# @author Ger van Diepen

#include "PhaseShift.h"
#include "DPBuffer.h"
#include "DPInfo.h"
#include "Exceptions.h"

#include "../Common/ParallelFor.h"
#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/BasicSL/Constants.h>

#include <iostream>
#include <iomanip>

#include <boost/algorithm/string/case_conv.hpp>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    PhaseShift::PhaseShift (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix)
      : itsInput   (input),
        itsName    (prefix),
        itsCenter  (parset.getStringVector(prefix+"phasecenter"))
    {}

    PhaseShift::PhaseShift (DPInput* input,
                            const ParameterSet& parset,
                            const string& prefix,
                            const vector<string>& defVal)
      : itsInput   (input),
        itsName    (prefix),
        itsCenter  (parset.getStringVector(prefix+"phasecenter", defVal))
    {}

    PhaseShift::~PhaseShift()
    {}

    void PhaseShift::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
      info().setMetaChanged();
      // Default phase center is the original one.
      MDirection newDir(itsInput->getInfo().phaseCenter());
      ////      bool original = true;
      bool original = false;
      if (! itsCenter.empty()) {
        newDir = handleCenter();
        original = false;
      }
      double newRa  = newDir.getValue().get()[0];
      double newDec = newDir.getValue().get()[1];
      double oldRa  = infoIn.phaseCenter().getValue().get()[0];
      double oldDec = infoIn.phaseCenter().getValue().get()[1];
      Matrix<double> oldUVW(3,3);
      Matrix<double> newUVW(3,3);
      fillTransMatrix (oldUVW, oldRa, oldDec);
      fillTransMatrix (newUVW, newRa, newDec);

      itsMat1.reference (product(transpose(newUVW), oldUVW));
      Matrix<double> wold(oldUVW(IPosition(2,0,2),IPosition(2,2,2)));
      Matrix<double> wnew(newUVW(IPosition(2,0,2),IPosition(2,2,2)));
      Matrix<double> tt= product(transpose(Matrix<double>(wold-wnew)), oldUVW);
      itsXYZ[0] = tt(0,0);
      itsXYZ[1] = tt(0,1);
      itsXYZ[2] = tt(0,2);
      ///      cout << itsXYZ[0]<<' '<<itsXYZ[1]<<' '<<itsXYZ[2]<<" ps"<<endl;

      info().setPhaseCenter (newDir, original);
      // Calculate 2*pi*freq/C to get correct phase term (in wavelengths).
      const Vector<double>& freq = infoIn.chanFreqs();
      itsFreqC.reserve (freq.size());
      for (uint i=0; i<freq.size(); ++i) {
        itsFreqC.push_back (2. * C::pi * freq[i] / C::c);
      }
      itsPhasors.resize (infoIn.nchan(), infoIn.nbaselines());
    }

    void PhaseShift::show (std::ostream& os) const
    {
      os << "PhaseShift " << itsName << std::endl;
      os << "  phasecenter:    " << itsCenter << std::endl;
    }

    void PhaseShift::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " PhaseShift " << itsName << endl;
    }

    bool PhaseShift::process (const DPBuffer& buf)
    {
      itsTimer.start();
      ///itsBuf.referenceFilled (buf);
      itsBuf.copy (buf);
      itsInput->fetchUVW (buf, itsBuf, itsTimer);
      int ncorr  = itsBuf.getData().shape()[0];
      int nchan  = itsBuf.getData().shape()[1];
      int nbl    = itsBuf.getData().shape()[2];
      assert (itsPhasors.nrow() == uint(nchan)  &&
                 itsPhasors.ncolumn() == uint(nbl));
      const double* mat1 = itsMat1.data();
      //# If ever in the future a time dependent phase center is used,
      //# the machine must be reset for each new time, thus each new call
      //# to process.
      ParallelFor<size_t> loop(getInfo().nThreads());
      loop.Run(0, nbl, [&](size_t bl, size_t /*thread*/) {
        Complex*  __restrict__ data    = itsBuf.getData().data() + bl*nchan*ncorr;
        double*   __restrict__ uvw     = itsBuf.getUVW().data() + bl*3;
        DComplex* __restrict__ phasors = itsPhasors.data() + bl*nchan;
        double u = uvw[0]*mat1[0] + uvw[1]*mat1[3] + uvw[2]*mat1[6];
        double v = uvw[0]*mat1[1] + uvw[1]*mat1[4] + uvw[2]*mat1[7];
        double w = uvw[0]*mat1[2] + uvw[1]*mat1[5] + uvw[2]*mat1[8];
        double phase = itsXYZ[0]*uvw[0] + itsXYZ[1]*uvw[1] + itsXYZ[2]*uvw[2];
        for (int j=0; j<nchan; ++j) {
          // Shift the phase of the data of this baseline.
          // Converting the phase term to wavelengths (and applying 2*pi)
          //      u_wvl = u_m / wvl = u_m * freq / c
                // has been done once in the beginning (in updateInfo).
          double phasewvl = phase * itsFreqC[j];
          DComplex phasor(cos(phasewvl), sin(phasewvl));
          *phasors++ = phasor;
          for (int k=0; k<ncorr; ++k) {
            *data = DComplex(*data) * phasor;
            data++;
          }
        }
        uvw[0] = u;
        uvw[1] = v;
        uvw[2] = w;
        uvw += 3;
      });
      itsTimer.stop();
      getNextStep()->process (itsBuf);
      return true;
    }

    void PhaseShift::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }

    MDirection PhaseShift::handleCenter()
    {
      // A case-insensitive name can be given for a moving source (e.g. SUN)
      // or a known source (e.g. CygA).
      if (itsCenter.size() == 1) {
        return MDirection::makeMDirection (itsCenter[0]);
      }
      // The phase center must be given in J2000 as two values (ra,dec).
      // In the future time dependent frames can be done as in UVWFlagger.
      if (itsCenter.size() != 2)
        throw Exception(
                 "2 values must be given in PhaseShift phasecenter");
      ///ASSERTSTR (itsCenter.size() < 4,
      ///"Up to 3 values can be given in PhaseShift phasecenter");
      casacore::MDirection phaseCenter;
      if (itsCenter.size() == 1) {
        string str = boost::to_upper_copy(itsCenter[0]);
        MDirection::Types tp;
        if (!MDirection::getType(tp, str))
          throw Exception(str + " is an invalid source type"
                   " in PhaseShift phasecenter");
        return MDirection(tp);
      }
      Quantity q0, q1;
      if (!MVAngle::read (q0, itsCenter[0]))
        throw Exception(itsCenter[0] + " is an invalid RA or longitude"
                 " in PhaseShift phasecenter");
      if (!MVAngle::read (q1, itsCenter[1]))
        throw Exception(itsCenter[1] + " is an invalid DEC or latitude"
                 " in PhaseShift phasecenter");
      MDirection::Types type = MDirection::J2000;
      if (itsCenter.size() > 2) {
        string str = boost::to_upper_copy(itsCenter[2]);
        MDirection::Types tp;
        if (!MDirection::getType(tp, str))
          throw Exception(str + " is an invalid direction type"
                   " in PhaseShift phasecenter");
      }
      return MDirection(q0, q1, type);
    }

    void PhaseShift::fillTransMatrix (Matrix<double>& mat,
                                      double ra, double dec)
    {
      assert (mat.nrow()==3 && mat.ncolumn()==3);
      double sinra  = sin(ra);
      double cosra  = cos(ra);
      double sindec = sin(dec);
      double cosdec = cos(dec);
      mat(0,0) = cosra;
      mat(1,0) = -sinra;
      mat(2,0) = 0;
      mat(0,1) = -sinra*sindec;
      mat(1,1) = -cosra*sindec;
      mat(2,1) = cosdec;
      mat(0,2) = sinra*cosdec;
      mat(1,2) = cosra*cosdec;
      mat(2,2) = sindec;
    }

  } //# end namespace
}
