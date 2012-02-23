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

#include <lofar_config.h>
#include <DPPP/PhaseShift.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <DPPP/ParSet.h>
#include <Common/LofarLogger.h>
#include <Common/StreamUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/BasicSL/Constants.h>
#include <iostream>
#include <iomanip>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    PhaseShift::PhaseShift (DPInput* input,
                            const ParSet& parset, const string& prefix)
      : itsInput   (input),
        itsName    (prefix),
        itsCenter  (parset.getStringVector(prefix+"phasecenter"))
    {}

    PhaseShift::PhaseShift (DPInput* input,
                            const ParSet& parset, const string& prefix,
                            const string& defVal)
      : itsInput   (input),
        itsName    (prefix),
        itsCenter  (parset.getStringVector(prefix+"phasecenter",
                                           vector<string>(1, defVal)))
    {}

    PhaseShift::~PhaseShift()
    {}

    void PhaseShift::updateInfo (DPInfo& info)
    {
      info.setNeedVisData();
      info.setNeedWrite();
      // Default phase center is the original one.
      MDirection newDir(itsInput->phaseCenter());
      ////      bool original = true;
      bool original = false;
      if (! itsCenter.empty()) {
        newDir = handleCenter();
        original = false;
      }
      double newRa  = newDir.getValue().get()[0];
      double newDec = newDir.getValue().get()[1];
      double oldRa  = info.phaseCenter().getValue().get()[0];
      double oldDec = info.phaseCenter().getValue().get()[1];
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

      info.setPhaseCenter (newDir, original);
      // Calculate 2*pi*freq/C to get correct phase term (in wavelengths).
      const Vector<double>& freq = itsInput->chanFreqs(info.nchanAvg());
      itsFreqC.reserve (freq.size());
      for (uint i=0; i<freq.size(); ++i) {
        itsFreqC.push_back (2. * C::pi * freq[i] / C::c);
      }
      itsPhasors.resize (info.nchan(), info.nbaselines());
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
      DPBuffer newBuf(buf);
      RefRows rowNrs(newBuf.getRowNrs());
      if (newBuf.getUVW().empty()) {
        newBuf.getUVW().reference (itsInput->fetchUVW (newBuf, rowNrs,
                                                       itsTimer));
      }
      newBuf.getData().unique();
      newBuf.getUVW().unique();
      uint ncorr  = newBuf.getData().shape()[0];
      uint nchan  = newBuf.getData().shape()[1];
      uint nbl    = newBuf.getData().shape()[2];
      DBGASSERT (itsPhasors.nrow() == nchan  &&  itsPhasors.ncolumn() == nbl);
      Complex* data     = newBuf.getData().data();
      double* uvw       = newBuf.getUVW().data();
      DComplex* phasors = itsPhasors.data();
      //# If ever in the future a time dependent phase center is used,
      //# the machine must be reset for each new time, thus each new call
      //# to process.
      for (uint i=0; i<nbl; ++i) {
        const double* mat1 = itsMat1.data();
        double u = uvw[0]*mat1[0] + uvw[1]*mat1[3] + uvw[2]*mat1[6];
        double v = uvw[0]*mat1[1] + uvw[1]*mat1[4] + uvw[2]*mat1[7];
        double w = uvw[0]*mat1[2] + uvw[1]*mat1[5] + uvw[2]*mat1[8];
        double phase = itsXYZ[0]*uvw[0] + itsXYZ[1]*uvw[1] + itsXYZ[2]*uvw[2];
        for (uint j=0; j<nchan; ++j) {
          // Shift the phase of the data of this baseline.
          // Convert the phase term to wavelengths (and apply 2*pi).
          // u_wvl = u_m / wvl = u_m * freq / c
          double phasewvl = phase * itsFreqC[j];
          DComplex phasor(cos(phasewvl), sin(phasewvl));
          *phasors++ = phasor;
          for (uint k=0; k<ncorr; ++k) {
            *data = DComplex(*data) * phasor;
            data++;
          }
        }
        uvw[0] = u;
        uvw[1] = v;
        uvw[2] = w;
        uvw += 3;
      }
      itsTimer.stop();
      getNextStep()->process (newBuf);
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
      ASSERTSTR (itsCenter.size() == 2,
                 "2 values must be given in PhaseShift phasecenter");
      ///ASSERTSTR (itsCenter.size() < 4,
      ///"Up to 3 values can be given in UVWFlagger phasecenter");
      MDirection phaseCenter;
      if (itsCenter.size() == 1) {
        string str = toUpper(itsCenter[0]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid source type"
                   " in UVWFlagger phasecenter");
        return MDirection(tp);
      }
      Quantity q0, q1;
      ASSERTSTR (MVAngle::read (q0, itsCenter[0]),
                 itsCenter[0] << " is an invalid RA or longitude"
                 " in UVWFlagger phasecenter");
      ASSERTSTR (MVAngle::read (q1, itsCenter[1]),
                 itsCenter[1] << " is an invalid DEC or latitude"
                 " in UVWFlagger phasecenter");
      MDirection::Types type = MDirection::J2000;
      if (itsCenter.size() > 2) {
        string str = toUpper(itsCenter[2]);
        MDirection::Types tp;
        ASSERTSTR (MDirection::getType(tp, str),
                   str << " is an invalid direction type in UVWFlagger"
                   " in UVWFlagger phasecenter");
      }
      return MDirection(q0, q1, type);
    }

    void PhaseShift::fillTransMatrix (Matrix<double>& mat,
                                      double ra, double dec)
    {
      DBGASSERT (mat.nrow()==3 && mat.ncolumn()==3);
      double sinra  = sin(ra);
      double sindec = sin(dec);
      double cosra  = cos(ra);
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
