//# ApplyCal.cc: DPPP step class to apply a calibration correction to the data
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
//# $Id: ApplyCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/ApplyCal.h>
#include <DPPP/DPBuffer.h>
#include <DPPP/DPInfo.h>
#include <Common/ParameterSet.h>
#include <Common/StringUtil.h>
#include <Common/LofarLogger.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/OS/File.h>
#include <iostream>
#include <iomanip>

using namespace casa;
using namespace LOFAR::BBS;

/// Look at BBSKernel MeasurementExprLOFARUtil.cc and Apply.cc

namespace LOFAR {
  namespace DPPP {

    ApplyCal::ApplyCal (DPInput* input,
                        const ParameterSet& parset,
                        const string& prefix)
      : itsInput       (input),
        itsName        (prefix),
        itsParmDBName  (parset.getString (prefix + "parmdb")),
        itsCorrectType (parset.getString (prefix + "correction")),
        itsSigma       (parset.getDouble (prefix + "sigma", 0.)),
        itsTimeInterval (-1),
        itsLastTime    (-1),
        itsUseAP       (False),
        itsNChan       (0),
        itsNPol        (0)
    {
      ASSERT (!itsParmDBName.empty());
      // Possible corrections one (or more?) of:
      //   Gain (real/imag or ampl/phase), RM, TEC, Clock, Bandpass
    }

    ApplyCal::~ApplyCal()
    {}

    void ApplyCal::updateInfo (const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setNeedWrite();
      itsTimeInterval = infoIn.timeInterval();

      itsParmDB.reset(new BBS::ParmFacade(itsParmDBName));

      /*
      vector<string> parNames;
      parNames=itsParmDB->getNames("*");
      for (int i=0;i<parNames.size();++i) {
        cout<<parNames[i]<<" ";
      }
      cout<<endl;
      */

      // Form the frequency axis for this time slot.
      int numFreqs = infoIn.chanFreqs().size();

      itsFreqInterval = infoIn.chanWidths()[0];
      itsMinFreq = infoIn.chanFreqs()[0]-0.5*itsFreqInterval;
      itsMaxFreq = infoIn.chanFreqs()[numFreqs-1]+0.5*itsFreqInterval;


      // Handle the correction type.

      string corrType = toLower(itsCorrectType);

      if (corrType == "gain") {
        itsUseAP = itsParmDB->getNames("Gain:*:Real:*").empty();
        itsHasCrossGain = !itsParmDB->getNames("Gain:0:1:").empty();
        itsParmExprs.push_back("Gain:0:0:");
        itsParmExprs.push_back("Gain:1:1:");
      } else if (corrType =="rm" || corrType == "rotationmeasure") {
        itsParmExprs.push_back("RotationMeasure:");
      } else if (corrType == "tec") {
        itsParmExprs.push_back("TEC:");
      } else if (corrType == "bandpass") { /*Bandpass:0:0 and Bandpass:1:1*/
        itsParmExprs.push_back("Bandpass:");
      } else {
        THROW (Exception, "Correction type " + itsCorrectType +
                         " is unknown");
      }
    }

    void ApplyCal::show (std::ostream& os) const
    {
      os << "ApplyCal " << itsName << std::endl;
      os << "  parmdb:         " << itsParmDBName << endl;
      os << "  correction:     " << itsCorrectType << endl;
      os << "  sigma:          " << itsSigma << endl;
    }

    void ApplyCal::showTimings (std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1 (os, itsTimer.getElapsed(), duration);
      os << " ApplyCal " << itsName << endl;
    }

    bool ApplyCal::process (const DPBuffer& bufin)
    {
      itsTimer.start();
      DPBuffer buf(bufin);
      buf.getData().unique();
      RefRows rowNrs(buf.getRowNrs());

      int numAnts = info().antennaNames().size();

      // If needed, cache parm values for the next 10 time slots.
      // getTime returns the center of the interval
      const int numparmbufsteps(10);
      double bufStartTime = buf.getTime() - 0.5*itsTimeInterval;

      if (buf.getTime() > itsLastTime) {
        itsLastTime = bufStartTime + numparmbufsteps * itsTimeInterval;

        map<string, vector<double> > parmMap;
        vector<vector<double> > oneParm;
        oneParm.reserve(info().antennaNames().size());

        for (int parmNum =0; parmNum<itsParmExprs.size();++parmNum) {
          parmMap = itsParmDB->getValuesMap( itsParmExprs[parmNum] + "*",
          itsMinFreq, itsMaxFreq, itsFreqInterval,
          bufStartTime, itsLastTime, itsTimeInterval, true);

          for (int ant = 0; ant < numAnts; ++ant) {
            // TODO: checken dat entry er is
            oneParm[ant]=parmMap.find("Gain:0:0:Real:"+info().antennaNames()[ant])->second;
          }

          itsParms[parmNum] = oneParm;
        }
      }

      // Loop through all baselines in the buffer.
      int nbl = bufin.getData().shape()[2];

      Complex* data = buf.getData().data();
      int npol  = buf.getData().shape()[0];
      int nchan = buf.getData().shape()[1];

      if (toLower(itsCorrectType)=="gain") {
        vector<vector<double> > gains00a;
        vector<vector<double> > gains00b;
        vector<vector<double> > gains11a;
        vector<vector<double> > gains11b;

        gains00a.reserve(info().antennaNames().size());
        gains00b.reserve(info().antennaNames().size());
        gains11a.reserve(info().antennaNames().size());
        gains11b.reserve(info().antennaNames().size());


        //#pragma omp parallel for
        for (int bl=0; bl<nbl; ++bl) {
          for (int i=bl*npol*nchan;i<(bl+1)*npol*nchan;i++) {
             // do nothing yet
          }
        }
      } // End of 'gain'


      itsTimer.stop();
      getNextStep()->process(buf);
      return false;
    }

    void ApplyCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


    // Corrections can be constant or can vary in freq.
    // TODO: this is a scalar effect, so should not be implemented as matrix
    void ApplyCal::applyTEC (Complex* vis, const DComplex& tec)
    {
      ///Matrix phase = (tec * -8.44797245e9) / freq;
      ///Matrix shift = tocomplex(cos(phase), sin(phase));
    }

    void ApplyCal::applyClock (Complex* vis, const DComplex* lhs,
                               const DComplex* rhs)
    {
      ///Matrix phase = freq * (delay() * casa::C::_2pi);
      ///Matrix shift = tocomplex(cos(phase), sin(phase));
      for (uint i=0; i<itsNChan; ++i) {
        DComplex factor = 1. / (lhs[i] * conj(rhs[i]));
        for (uint j=0; j<itsNPol; ++j) {
          *vis++ *= factor;
        }
      }
    }

    void ApplyCal::applyBandpass (Complex* vis, const DComplex* lhs,
                                  const DComplex* rhs)
    {
    }

    void ApplyCal::applyRM (Complex* vis, const DComplex* lhs,
                            const DComplex* rhs)
    {
      // Precompute lambda squared for the current frequency point.
      /*
      const double lambda = C::c / grid[FREQ]->center(f);
      const double lambda2 = lambda * lambda;

      double *sample = origin + f;
      for (unsigned int t = 0; t < nTime; ++t) {
        *sample = lambda2;
        sample += nFreq;
      }

      Matrix chi = rm() * lambdaSqr;
      Matrix cosChi = cos(chi);
      Matrix sinChi = sin(chi);

      JonesMatrix::View result;
      result.assign(0, 0, cosChi);
      result.assign(0, 1, -sinChi);
      result.assign(1, 0, sinChi);
      result.assign(1, 1, cosChi);
      */
    }

    void ApplyCal::applyJones (Complex* vis, const DComplex* lhs,
                               const DComplex* rhs)
    {
      for (uint i=0; i<itsNChan; ++i) {
        // Compute the Mueller matrix.

        DComplex mueller[4][4];
        mueller[0][0] = lhs[0] * conj(rhs[0]);
        mueller[0][1] = lhs[0] * conj(rhs[1]);
        mueller[1][0] = lhs[0] * conj(rhs[2]);
        mueller[1][1] = lhs[0] * conj(rhs[3]);

        mueller[0][2] = lhs[1] * conj(rhs[0]);
        mueller[0][3] = lhs[1] * conj(rhs[1]);
        mueller[1][2] = lhs[1] * conj(rhs[2]);
        mueller[1][3] = lhs[1] * conj(rhs[3]);

        mueller[2][0] = lhs[2] * conj(rhs[0]);
        mueller[2][1] = lhs[2] * conj(rhs[1]);
        mueller[3][0] = lhs[2] * conj(rhs[2]);
        mueller[3][1] = lhs[2] * conj(rhs[3]);

        mueller[2][2] = lhs[3] * conj(rhs[0]);
        mueller[2][3] = lhs[3] * conj(rhs[1]);
        mueller[3][2] = lhs[3] * conj(rhs[2]);
        mueller[3][3] = lhs[3] * conj(rhs[3]);

        // Apply Mueller matrix to visibilities.
        DComplex xx, xy, yx, yy;
        xx = (mueller[0][0] * DComplex(vis[0]) +
              mueller[0][1] * DComplex(vis[1]) +
              mueller[0][2] * DComplex(vis[2]) +
              mueller[0][3] * DComplex(vis[3]));

        xy = (mueller[1][0] * DComplex(vis[0]) +
              mueller[1][1] * DComplex(vis[1]) +
              mueller[1][2] * DComplex(vis[2]) +
              mueller[1][3] * DComplex(vis[3]));

        yx = (mueller[2][0] * DComplex(vis[0]) +
              mueller[2][1] * DComplex(vis[1]) +
              mueller[2][2] * DComplex(vis[2]) +
              mueller[2][3] * DComplex(vis[3]));

        yy = (mueller[3][0] * DComplex(vis[0]) +
              mueller[3][1] * DComplex(vis[1]) +
              mueller[3][2] * DComplex(vis[2]) +
              mueller[3][3] * DComplex(vis[3]));

        vis[0] = xx;
        vis[1] = xy;
        vis[2] = yx;
        vis[3] = yy;
        vis += 4;

      }
    }

    // Inverts complex input matrix (in place??)
    // TODO: what does this sigma term do? It is added, should it be added back?
    void ApplyCal::invert (DComplex* v)
    {
      // Add the variance of the nuisance term to the elements on the diagonal.
      const double variance = 0;//itsSigma * itsSigma;
      DComplex v0 = v[0] + variance;
      DComplex v3 = v[3] + variance;
      // Compute inverse in the usual way.
      DComplex invDet(1.0 / (v0 * v3 - v[1] * v[2]));
      v[0] = v3 * invDet;
      v[2] = v[2] * -invDet;
      v[1] = v[1] * -invDet;
      v[3] = v0 * invDet;
    }

  } //# end namespace
}
