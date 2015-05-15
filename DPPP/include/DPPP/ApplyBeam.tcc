#ifndef DPPP_ApplyBeam_TCC
#define DPPP_ApplyBeam_TCC

#include <lofar_config.h>
#include <DPPP/ApplyBeam.h>

#include <iostream>
//#include <iomanip>
#include <Common/ParameterSet.h>
#include <Common/Timer.h>
#include <Common/OpenMP.h>
#include <ParmDB/ParmDBMeta.h>
#include <ParmDB/PatchInfo.h>
#include <DPPP/DPInfo.h>
#include <DPPP/FlagCounter.h>
#include <DPPP/Position.h>
#include <DPPP/Simulator.h>
#include <DPPP/Simulate.h>
#include <DPPP/ApplyCal.h>

#include <DPPP/Stokes.h>
#include <DPPP/PointSource.h>
#include <DPPP/GaussianSource.h>
#include <ParmDB/SourceDB.h>

#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <tables/Tables/RefRows.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

namespace LOFAR {
  namespace DPPP {
template<typename T>
void ApplyBeam::applyBeam(
        const DPInfo& info, double time, T* data0,
        const StationResponse::vector3r_t& srcdir,
        const StationResponse::vector3r_t& refdir,
        const StationResponse::vector3r_t& tiledir,
        const vector<StationResponse::Station::Ptr>& antBeamInfo,
        vector<StationResponse::matrix22c_t>& beamValues, bool useChannelFreq,
        bool invert)
    { 
      // Get the beam values for each station.
      uint nCh = info.chanFreqs().size();
      uint nSt = beamValues.size() / nCh;
      uint nBl = info.nbaselines();

      double reffreq=info.refFreq();

      // Apply the beam values of both stations to the ApplyBeamed data.
      dcomplex tmp[4];
      for (size_t ch = 0; ch < nCh; ++ch) {
        if (useChannelFreq) {
          reffreq=info.chanFreqs()[ch];
        }

        // Fill beamValues for channel ch
        for (size_t st = 0; st < nSt; ++st) {
          beamValues[nCh * st + ch] = antBeamInfo[st]->response(time,
                                        info.chanFreqs()[ch], srcdir,
                                        reffreq, refdir, tiledir);
          if (invert) {
            ApplyCal::invert((dcomplex*)(&(beamValues[nCh * st])));
          }
        }

        // Apply beam for channel ch on all baselines
        for (size_t bl = 0; bl < nBl; ++bl) {
          T* data = data0 + bl * 4 * nCh + ch * 4;
          StationResponse::matrix22c_t *left = &(beamValues[nCh
              * info.getAnt1()[bl]]);
          StationResponse::matrix22c_t *right = &(beamValues[nCh
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
        }
      }
    }
  }
}

#endif
