//# ApplyBeam.tcc: DPPP step class to ApplyBeam visibilities from a source model
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
//# $Id:
//#
//# @author Tammo Jan Dijkema

#ifndef DPPP_ApplyBeam_TCC
#define DPPP_ApplyBeam_TCC

#include "ApplyBeam.h"

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"
#include "../Common/OpenMP.h"

#include "../ParmDB/ParmDBMeta.h"
#include "../ParmDB/PatchInfo.h"
#include "../ParmDB/SourceDB.h"

#include "DPInfo.h"
#include "FlagCounter.h"
#include "Position.h"
#include "Simulator.h"
#include "Simulate.h"
#include "ApplyCal.h"

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/AntennaField.h>
#endif

#include "Stokes.h"
#include "PointSource.h"
#include "GaussianSource.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/tables/Tables/RefRows.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

namespace DP3 {
  namespace DPPP {


// applyBeam is templated on the type of the data, could be double or float
template<typename T>
void ApplyBeam::applyBeam(
        const DPInfo& info, double time, T* data0, float* weight0,
        const LOFAR::StationResponse::vector3r_t& srcdir,
        const LOFAR::StationResponse::vector3r_t& refdir,
        const LOFAR::StationResponse::vector3r_t& tiledir,
        const vector<LOFAR::StationResponse::Station::Ptr>& antBeamInfo,
        vector<LOFAR::StationResponse::matrix22c_t>& beamValues, bool useChannelFreq,
        bool invert, int mode, bool doUpdateWeights)
    { 
      // Get the beam values for each station.
      uint nCh = info.chanFreqs().size();
      uint nSt = beamValues.size() / nCh;
      uint nBl = info.nbaselines();

      // Store array factor in diagonal matrix (in other modes this variable
      // is not used).
      LOFAR::StationResponse::diag22c_t af_tmp;

      double reffreq=info.refFreq();

      // Apply the beam values of both stations to the ApplyBeamed data.
      dcomplex tmp[4];
      for (size_t ch = 0; ch < nCh; ++ch) {
        if (useChannelFreq) {
          reffreq=info.chanFreqs()[ch];
        }

        switch (mode) {
        case DEFAULT:
          // Fill beamValues for channel ch
          for (size_t st = 0; st < nSt; ++st) {
            beamValues[nCh * st + ch] = antBeamInfo[st]->response(time,
                                          info.chanFreqs()[ch], srcdir,
                                          reffreq, refdir, tiledir);
            if (invert) {
              ApplyCal::invert((dcomplex*)(&(beamValues[nCh * st + ch])));
            }
          }
          break;
        case ARRAY_FACTOR:
          // Fill beamValues for channel ch
          for (size_t st = 0; st < nSt; ++st) {
            af_tmp = antBeamInfo[st]->arrayFactor(time,
                                          info.chanFreqs()[ch], srcdir,
                                          reffreq, refdir, tiledir);
            beamValues[nCh * st + ch][0][1]=0.;
            beamValues[nCh * st + ch][1][0]=0.;

            if (invert) {
              beamValues[nCh * st + ch][0][0]=1./af_tmp[0];
              beamValues[nCh * st + ch][1][1]=1./af_tmp[1];
            } else {
              beamValues[nCh * st + ch][0][0]=af_tmp[0];
              beamValues[nCh * st + ch][1][1]=af_tmp[1];
            }
          }
          break;
        case ELEMENT:
          // Fill beamValues for channel ch
          for (size_t st = 0; st < nSt; ++st) {
            LOFAR::StationResponse::AntennaField::ConstPtr field =
                *(antBeamInfo[st]->beginFields());

            beamValues[nCh * st + ch] = field->elementResponse(time,
                                                 info.chanFreqs()[ch],
                                                 srcdir);
            if (invert) {
              ApplyCal::invert((dcomplex*)(&(beamValues[nCh * st + ch])));
            }
          }
          break;
        }

        // Apply beam for channel ch on all baselines
        // For mode=ARRAY_FACTOR, too much work is done here because we know
        // that r and l are diagonal
        for (size_t bl = 0; bl < nBl; ++bl) {
          T* data = data0 + bl * 4 * nCh + ch * 4;
          LOFAR::StationResponse::matrix22c_t *left = &(beamValues[nCh
              * info.getAnt1()[bl]]);
          LOFAR::StationResponse::matrix22c_t *right = &(beamValues[nCh
              * info.getAnt2()[bl]]);
          dcomplex l[] = { left[ch][0][0], left[ch][0][1],
                           left[ch][1][0], left[ch][1][1] };
          // Form transposed conjugate of right.
          dcomplex r[] = { conj(right[ch][0][0]), conj(right[ch][1][0]),
                           conj(right[ch][0][1]), conj(right[ch][1][1]) };
          // left*data
          tmp[0] = l[0] * dcomplex(data[0]) + l[1] * dcomplex(data[2]);
          tmp[1] = l[0] * dcomplex(data[1]) + l[1] * dcomplex(data[3]);
          tmp[2] = l[2] * dcomplex(data[0]) + l[3] * dcomplex(data[2]);
          tmp[3] = l[2] * dcomplex(data[1]) + l[3] * dcomplex(data[3]);
          // data*conj(right)
          data[0] = tmp[0] * r[0] + tmp[1] * r[2];
          data[1] = tmp[0] * r[1] + tmp[1] * r[3];
          data[2] = tmp[2] * r[0] + tmp[3] * r[2];
          data[3] = tmp[2] * r[1] + tmp[3] * r[3];

          if (doUpdateWeights) {
            ApplyCal::applyWeights(l, r, weight0 + bl * 4 * nCh + ch * 4);
          }
        }
      }
    }
  }
}

#endif
