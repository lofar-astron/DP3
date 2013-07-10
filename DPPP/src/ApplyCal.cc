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
        itsBufStep     (0),
        itsSigma       (parset.getDouble (prefix + "sigma", 0.)),
        itsTimeInterval (-1),
        itsLastTime    (-1),
        itsUseAP       (false),
        itsHasCrossGain (false)
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

      // Handle the correction type.

      string corrType = toLower(itsCorrectType);

      if (corrType == "gain") {
        itsUseAP = itsParmDB->getNames("Gain:*:Real:*").empty();
        itsHasCrossGain = !itsParmDB->getNames("Gain:0:1:").empty();
        if (itsUseAP) {
          itsParmExprs.push_back("Gain:0:0:Ampl:");
          itsParmExprs.push_back("Gain:0:0:Phase:");
          itsParmExprs.push_back("Gain:1:1:Ampl:");
          itsParmExprs.push_back("Gain:1:1:Phase:");
          if (itsHasCrossGain) {
            itsParmExprs.push_back("Gain:0:1:Ampl:");
            itsParmExprs.push_back("Gain:0:1:Phase:");
            itsParmExprs.push_back("Gain:1:0:Ampl:");
            itsParmExprs.push_back("Gain:1:0:Phase:");
          }
        } else {
          itsParmExprs.push_back("Gain:0:0:Real:");
          itsParmExprs.push_back("Gain:0:0:Imag:");
          itsParmExprs.push_back("Gain:1:1:Real:");
          itsParmExprs.push_back("Gain:1:1:Imag:");
          if (itsHasCrossGain) {
            itsParmExprs.push_back("Gain:0:1:Real:");
            itsParmExprs.push_back("Gain:0:1:Imag:");
            itsParmExprs.push_back("Gain:1:0:Real:");
            itsParmExprs.push_back("Gain:1:0:Imag:");
          }
        }
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

      itsParms.resize(itsParmExprs.size());
      for (uint i=0;i<itsParms.size();++i) {
        itsParms[i].resize(info().antennaNames().size());
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

      double bufStartTime = buf.getTime() - 0.5*itsTimeInterval;

      if (buf.getTime() > itsLastTime) {
        updateParms(bufStartTime);
        itsBufStep=0;
      }
      else {
        itsBufStep++;
      }

      // Loop through all baselines in the buffer.
      int nbl = bufin.getData().shape()[2];

      Complex* data = buf.getData().data();
      int npol  = buf.getData().shape()[0];
      ASSERT(npol==4);
      int nchan = buf.getData().shape()[1];

      //#pragma omp parallel for
      for (int bl=0; bl<nbl; ++bl) {
        for (int chan=0;chan<nchan;chan++) {
            applyGain( &data[bl * npol * nchan + chan * npol ],
              info().getAnt1()[bl], info().getAnt2()[bl], chan, itsBufStep);
        }
      }

      itsTimer.stop();
      getNextStep()->process(buf);
      return false;
    }

    void ApplyCal::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }


    void ApplyCal::updateParms (const double bufStartTime)
    {
      int numAnts = info().antennaNames().size();

      // If needed, cache parm values for the next 10 time slots.
      const int numparmbufsteps(10);

      int numFreqs         (info().chanFreqs().size());
      double freqInterval  (info().chanWidths()[0]);
      double minFreq       (info().chanFreqs()[0]-0.5*freqInterval);
      double maxFreq (info().chanFreqs()[numFreqs-1]+0.5*freqInterval);

      itsLastTime = bufStartTime + numparmbufsteps * itsTimeInterval;

      map<string, vector<double> > parmMap;
      map<string, vector<double> >::iterator parmIt;

      for (uint parmNum =0; parmNum<itsParmExprs.size();++parmNum) {
        parmMap = itsParmDB->getValuesMap( itsParmExprs[parmNum] + "*",
        minFreq, maxFreq, freqInterval,
        bufStartTime, itsLastTime, itsTimeInterval, true);

        for (int ant = 0; ant < numAnts; ++ant) {
          parmIt = parmMap.find(
                    itsParmExprs[parmNum]+info().antennaNames()[ant]);
          ASSERT( parmIt != parmMap.end() );

          itsParms[parmNum][ant].swap(parmIt->second);
        }
      }
    }


    void ApplyCal::applyGain (Complex* vis, const int ant1, const int ant2,
        int chan, int time) {
      vector<vector<Complex> > gainA(2);
      vector<vector<Complex> > gainB(2);

      int timeFreqOffset=(time*info().nchan())+chan;

      for (uint i=0;i<2;i++) {
        gainA[i].resize(2);
      }

      if (itsUseAP) {
        gainA[0][0] = polar(itsParms[0][ant1][timeFreqOffset],
                        itsParms[1][ant1][timeFreqOffset]);
        gainA[1][1] = polar(itsParms[2][ant1][timeFreqOffset],
                        itsParms[3][ant1][timeFreqOffset]);
        gainB[0][0] = polar(itsParms[0][ant2][timeFreqOffset],
                        itsParms[1][ant2][timeFreqOffset]);
        gainB[1][1] = polar(itsParms[2][ant2][timeFreqOffset],
                        itsParms[3][ant2][timeFreqOffset]);
        if (itsHasCrossGain) {
          gainA[0][1] = polar(itsParms[4][ant1][timeFreqOffset],
                          itsParms[5][ant1][timeFreqOffset]);
          gainA[1][0] = polar(itsParms[6][ant1][timeFreqOffset],
                          itsParms[7][ant1][timeFreqOffset]);
          gainB[0][1] = polar(itsParms[4][ant2][timeFreqOffset],
                          itsParms[5][ant2][timeFreqOffset]);
          gainB[1][0] = polar(itsParms[6][ant2][timeFreqOffset],
                          itsParms[7][ant2][timeFreqOffset]);
        }
      } else {
        gainA[0][0] = Complex(itsParms[0][ant1][timeFreqOffset],
                          itsParms[1][ant1][timeFreqOffset]);
        gainA[1][1] = Complex(itsParms[2][ant1][timeFreqOffset],
                          itsParms[3][ant1][timeFreqOffset]);
        gainB[0][0] = Complex(itsParms[0][ant2][timeFreqOffset],
                          itsParms[1][ant2][timeFreqOffset]);
        gainB[1][1] = Complex(itsParms[2][ant2][timeFreqOffset],
                          itsParms[3][ant2][timeFreqOffset]);
        if (itsHasCrossGain) {
          gainA[0][1] = Complex(itsParms[4][ant1][timeFreqOffset],
                            itsParms[5][ant1][timeFreqOffset]);
          gainA[1][0] = Complex(itsParms[6][ant1][timeFreqOffset],
                            itsParms[7][ant1][timeFreqOffset]);
          gainB[0][1] = Complex(itsParms[4][ant2][timeFreqOffset],
                            itsParms[5][ant2][timeFreqOffset]);
          gainB[1][0] = Complex(itsParms[6][ant2][timeFreqOffset],
                            itsParms[7][ant2][timeFreqOffset]);
        }
      }

      if (itsHasCrossGain) {
        THROW ( Exception, "cross gain not implemented");
        vector<vector<Complex> > gainAxvis(2);

        for (uint i=0;i<2;++i) {
          gainAxvis.resize(2);
          for (uint j=0;j<2;++j) {
            gainAxvis[i][j] = 0;
            for (uint k=0;k<2;++k) {
              gainAxvis[i][j] += gainA[i][j+k] * vis[2 * (i+k) + j];
            }
          }
        }

        for (uint i=0;i<2;++i) {
          for (uint j=0;j<2;++j) {
            vis[2 * i + j] = 0;
            for (uint k=0;k<2;++k) {
              vis[2 * i + j] += gainAxvis[i][j+k] * conj(gainB[i+k][j]);
            }
          }
        }
      }
      else {
        vis[0] *= gainA[0][0] * conj(gainB[0][0]);
        vis[1] *= gainA[0][0] * conj(gainB[1][1]);
        vis[2] *= gainA[1][1] * conj(gainB[0][0]);
        vis[3] *= gainA[1][1] * conj(gainB[1][1]);
      }
    }

    // Corrections can be constant or can vary in freq.
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
      /*
        DComplex factor = 1. / (lhs[i] * conj(rhs[i]));
        for (uint j=0; j<itsNPol; ++j) {
          *vis++ *= factor;
        }
      */
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
