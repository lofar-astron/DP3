//# ApplyBeam.cc: DPPP step class to ApplyBeam visibilities
//# Copyright (C) 2015
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
//# $Id: GainCal.cc 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Tammo Jan Dijkema

#include <iostream>

#include "../Common/ParameterSet.h"
#include "../Common/Timer.h"
#include "../Common/StreamUtil.h"
#include "../Common/StringUtil.h"

#include "ApplyBeam.h"
#include "ApplyCal.h" // for matrix inversion
#include "DPInfo.h"
#include "Exceptions.h"
#include "FlagCounter.h"
#include "Position.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>

#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>

#include <casa/Quanta/MVAngle.h>

using namespace casacore;

namespace DP3 {
  namespace DPPP {

    ApplyBeam::ApplyBeam(DPInput* input, const ParameterSet& parset,
                         const string& prefix, bool substep)
        :
          itsInput(input),
          itsName(prefix),
          itsUpdateWeights(parset.getBool(prefix + "updateweights", false)),
          itsDirectionStr(parset.getStringVector(prefix+"direction", std::vector<std::string>())),
          itsUseChannelFreq(parset.getBool(prefix + "usechannelfreq", true)),
          itsDebugLevel(parset.getInt(prefix + "debuglevel", 0))
    {
      // only read 'invert' parset key if it is a separate step
      // if applybeam is called from gaincal/predict, the invert key should always be false
      if (substep) {
        itsInvert=false;
      } else {
        itsInvert=parset.getBool(prefix + "invert", true);
      }
      string mode=boost::to_lower_copy(parset.getString(prefix + "beammode","default"));
      assert (mode=="default" || mode=="array_factor" || mode=="element");
      if (mode=="default") {
        itsMode=FullBeamCorrection;
      } else if (mode=="array_factor") {
        itsMode=ArrayFactorBeamCorrection;
      } else if (mode=="element") {
        itsMode=ElementBeamCorrection;
      } else {
        throw Exception("Beammode should be DEFAULT, ARRAY_FACTOR or ELEMENT");
      }
    }

    ApplyBeam::ApplyBeam()
    {
    }

    ApplyBeam::~ApplyBeam()
    {
    }

    void ApplyBeam::updateInfo(const DPInfo& infoIn)
    {
      info() = infoIn;
      info().setNeedVisData();
      info().setWriteData();
      if (itsUpdateWeights) {
        info().setWriteWeights();
      }
      
      // Parse direction parset value
      if (itsDirectionStr.empty())
        itsDirection = info().phaseCenter();
      else {
        if (itsDirectionStr.size() != 2)
          throw std::runtime_error("2 values must be given in ApplyBeam");
        casacore::MDirection phaseCenter;
        Quantity q0, q1;
        if (!MVAngle::read (q0, itsDirectionStr[0]))
          throw Exception(itsDirectionStr[0] + " is an invalid RA or longitude in ApplyBeam direction");
        if (!MVAngle::read (q1, itsDirectionStr[1]))
          throw Exception(itsDirectionStr[1] + " is an invalid DEC or latitude in ApplyBeam direction");
        MDirection::Types type = MDirection::J2000;
        itsDirection = MDirection(q0, q1, type);
      }
      
      if(info().beamCorrectionMode() != NoBeamCorrection)
        throw std::runtime_error("In applying the beam: the metadata of this observation indicate that the beam has already been applied");
      info().setBeamCorrectionMode(itsMode);
      info().setBeamCorrectionDir(itsDirection);

      const size_t nSt = info().nantenna();
      const size_t nCh = info().nchan();

      const size_t nThreads = getInfo().nThreads();
      itsBeamValues.resize(nThreads);

      // Create the Measure ITRF conversion info given the array position.
      // The time and direction are filled in later.
      itsMeasConverters.resize(nThreads);
      itsMeasFrames.resize(nThreads);
      itsAntBeamInfo.resize(nThreads);

      for (size_t thread = 0; thread < nThreads; ++thread) {
        itsBeamValues[thread].resize(nSt * nCh);
        itsMeasFrames[thread].set(info().arrayPosCopy());
        itsMeasFrames[thread].set(
            MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
        itsMeasConverters[thread].set(
            MDirection::J2000,
            MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
        itsInput->fillBeamInfo(itsAntBeamInfo[thread], info().antennaNames());
      }
    }

    void ApplyBeam::show(std::ostream& os) const
    {
      os << "ApplyBeam " << itsName << '\n';
      os << "  mode:              ";
      if (itsMode==FullBeamCorrection)
        os<<"default";
      else if (itsMode==ArrayFactorBeamCorrection)
        os<<"array_factor";
      else os<<"element";
      os << '\n';
      os << "  use channelfreq:   " << boolalpha << itsUseChannelFreq << '\n';
      os << "  direction:         " << itsDirectionStr << '\n';
      os << "  invert:            " << boolalpha << itsInvert << '\n';
      os << "  update weights:    " << boolalpha << itsUpdateWeights << '\n';
    }

    void ApplyBeam::showTimings(std::ostream& os, double duration) const
    {
      os << "  ";
      FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
      os << " ApplyBeam " << itsName << '\n';
    }

    bool ApplyBeam::processMultithreaded(const DPBuffer& bufin, size_t thread)
    {
      itsTimer.start();
      itsBuffer.copy (bufin);
      Complex* data=itsBuffer.getData().data();

      if (itsUpdateWeights) {
        itsInput->fetchWeights (bufin, itsBuffer, itsTimer);
      }
      float* weight = itsBuffer.getWeights().data();

      double time = itsBuffer.getTime();

      //Set up directions for beam evaluation
      LOFAR::StationResponse::vector3r_t refdir, tiledir, srcdir;

      /**
       * I'm not sure this is correct the way it is. These loops
       * seem to initialize variables that are never used in a
       * multi-threaded way, and if they were used from multiple
       * threads, it would imply process() is called multiple times,
       * and hence this initialization is already subject to a race
       * condition... ???
       * itsMeasFrames seems not to be actually used.
       * André, 2018-10-07
       */
      for (size_t threadIter = 0; threadIter < getInfo().nThreads(); ++threadIter) {
        itsMeasFrames[threadIter].resetEpoch(
            MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
        //Do a conversion on all threads, because converters are not
        //thread safe and apparently need to be used at least once
        refdir = dir2Itrf(info().delayCenter(), itsMeasConverters[threadIter]);
        tiledir = dir2Itrf(info().tileBeamDir(), itsMeasConverters[threadIter]);
        srcdir = dir2Itrf(itsDirection, itsMeasConverters[threadIter]);
      }

      applyBeam(info(), time, data, weight, srcdir, refdir, tiledir,
                itsAntBeamInfo[thread], itsBeamValues[thread],
                itsUseChannelFreq, itsInvert, itsMode, itsUpdateWeights);

      itsTimer.stop();
      getNextStep()->process(itsBuffer);
      return false;
    }

    LOFAR::StationResponse::vector3r_t ApplyBeam::dir2Itrf(
        const MDirection& dir, MDirection::Convert& measConverter)
    {
      const MDirection& itrfDir = measConverter(dir);
      const Vector<Double>& itrf = itrfDir.getValue().getValue();
      LOFAR::StationResponse::vector3r_t vec;
      vec[0] = itrf[0];
      vec[1] = itrf[1];
      vec[2] = itrf[2];
      return vec;
    }

    void ApplyBeam::finish()
    {
      // Let the next steps finish.
      getNextStep()->finish();
    }
    
// applyBeam is templated on the type of the data, could be complex<double> or complex<float>
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
  using dcomplex = std::complex<double>;
  // Get the beam values for each station.
  uint nCh = info.chanFreqs().size();
  uint nSt = beamValues.size() / nCh;
  uint nBl = info.nbaselines();

  // Store array factor in diagonal matrix (in other modes this variable
  // is not used).
  LOFAR::StationResponse::diag22c_t af_tmp;

  double reffreq=info.refFreq();

  // Apply the beam values of both stations to the ApplyBeamed data.
  std::complex<double> tmp[4];
  for (size_t ch = 0; ch < nCh; ++ch) {
    if (useChannelFreq) {
      reffreq=info.chanFreqs()[ch];
    }

    switch (mode) {
    case FullBeamCorrection:
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
    case ArrayFactorBeamCorrection:
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
    case ElementBeamCorrection:
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

template
void ApplyBeam::applyBeam(const DPInfo& info, double time, std::complex<double>* data0, float* weight0,
  const LOFAR::StationResponse::vector3r_t& srcdir,
  const LOFAR::StationResponse::vector3r_t& refdir,
  const LOFAR::StationResponse::vector3r_t& tiledir,
  const vector<LOFAR::StationResponse::Station::Ptr>& antBeamInfo,
  vector<LOFAR::StationResponse::matrix22c_t>& beamValues, bool useChannelFreq,
  bool invert, int mode, bool doUpdateWeights);
    
}} //# end namespaces

