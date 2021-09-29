// PhaseShift.cc: DPPP step class to shift the data to another phase center
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "PhaseShift.h"

#include "../base/Direction.h"
#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <aocommon/parallelfor.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/BasicSL/Constants.h>

#include <iostream>
#include <iomanip>

#include <boost/algorithm/string/case_conv.hpp>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::common::operator<<;

using casacore::IPosition;
using casacore::Matrix;
using casacore::MDirection;

namespace dp3 {
namespace steps {

PhaseShift::PhaseShift(InputStep* input, const common::ParameterSet& parset,
                       const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsCenter(parset.getStringVector(prefix + "phasecenter")) {}

PhaseShift::PhaseShift(InputStep* input, const common::ParameterSet& parset,
                       const string& prefix,
                       const std::vector<std::string>& defVal)
    : itsInput(input),
      itsName(prefix),
      itsCenter(parset.getStringVector(prefix + "phasecenter", defVal)) {}

PhaseShift::~PhaseShift() {}

void PhaseShift::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  info().setMetaChanged();
  // Default phase center is the original one.
  MDirection newDir(itsInput->getInfo().phaseCenter());
  ////      bool original = true;
  bool original = false;
  if (!itsCenter.empty()) {
    newDir = handleCenter();
    original = false;
  }
  const base::Direction new_direction(newDir.getValue().get()[0],
                                      newDir.getValue().get()[1]);
  const base::Direction old_direction(infoIn.phaseCenter().getValue().get()[0],
                                      infoIn.phaseCenter().getValue().get()[1]);
  Matrix<double> oldUVW(3, 3);
  Matrix<double> newUVW(3, 3);
  fillTransMatrix(oldUVW, old_direction);
  fillTransMatrix(newUVW, new_direction);

  itsMat1.reference(product(transpose(newUVW), oldUVW));
  Matrix<double> wold(oldUVW(IPosition(2, 0, 2), IPosition(2, 2, 2)));
  Matrix<double> wnew(newUVW(IPosition(2, 0, 2), IPosition(2, 2, 2)));
  Matrix<double> tt = product(transpose(Matrix<double>(wold - wnew)), oldUVW);
  itsXYZ[0] = tt(0, 0);
  itsXYZ[1] = tt(0, 1);
  itsXYZ[2] = tt(0, 2);
  ///      cout << itsXYZ[0]<<' '<<itsXYZ[1]<<' '<<itsXYZ[2]<<" ps"<<'\n';

  info().setPhaseCenter(newDir, original);
  // Calculate 2*pi*freq/C to get correct phase term (in wavelengths).
  const casacore::Vector<double>& freq = infoIn.chanFreqs();
  itsFreqC.reserve(freq.size());
  for (unsigned int i = 0; i < freq.size(); ++i) {
    itsFreqC.push_back(2. * casacore::C::pi * freq[i] / casacore::C::c);
  }
  itsPhasors.resize(infoIn.nchan(), infoIn.nbaselines());
}

void PhaseShift::show(std::ostream& os) const {
  os << "PhaseShift " << itsName << '\n';
  os << "  phasecenter:    " << itsCenter << '\n';
}

void PhaseShift::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " PhaseShift " << itsName << '\n';
}

bool PhaseShift::process(const DPBuffer& buf) {
  itsTimer.start();
  /// itsBuf.referenceFilled (buf);
  itsBuf.copy(buf);
  itsInput->fetchUVW(buf, itsBuf, itsTimer);
  int ncorr = itsBuf.getData().shape()[0];
  int nchan = itsBuf.getData().shape()[1];
  int nbl = itsBuf.getData().shape()[2];
  const double* mat1 = itsMat1.data();
  // If ever in the future a time dependent phase center is used,
  // the machine must be reset for each new time, thus each new call
  // to process.
  aocommon::ParallelFor<size_t> loop(getInfo().nThreads());
  loop.Run(0, nbl, [&](size_t bl, size_t /*thread*/) {
    casacore::Complex* __restrict__ data =
        itsBuf.getData().data() + bl * nchan * ncorr;
    double* __restrict__ uvw = itsBuf.getUVW().data() + bl * 3;
    casacore::DComplex* __restrict__ phasors = itsPhasors.data() + bl * nchan;
    double u = uvw[0] * mat1[0] + uvw[1] * mat1[3] + uvw[2] * mat1[6];
    double v = uvw[0] * mat1[1] + uvw[1] * mat1[4] + uvw[2] * mat1[7];
    double w = uvw[0] * mat1[2] + uvw[1] * mat1[5] + uvw[2] * mat1[8];
    double phase = itsXYZ[0] * uvw[0] + itsXYZ[1] * uvw[1] + itsXYZ[2] * uvw[2];
    for (int j = 0; j < nchan; ++j) {
      // Shift the phase of the data of this baseline.
      // Converting the phase term to wavelengths (and applying 2*pi)
      //      u_wvl = u_m / wvl = u_m * freq / c
      // has been done once in the beginning (in updateInfo).
      double phasewvl = phase * itsFreqC[j];
      casacore::DComplex phasor(cos(phasewvl), sin(phasewvl));
      *phasors++ = phasor;
      for (int k = 0; k < ncorr; ++k) {
        *data = casacore::DComplex(*data) * phasor;
        data++;
      }
    }
    uvw[0] = u;
    uvw[1] = v;
    uvw[2] = w;
    uvw += 3;
  });
  itsTimer.stop();
  getNextStep()->process(itsBuf);
  return true;
}

void PhaseShift::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

MDirection PhaseShift::handleCenter() {
  // A case-insensitive name can be given for a moving source (e.g. SUN)
  // or a known source (e.g. CygA).
  if (itsCenter.size() == 1) {
    return MDirection::makeMDirection(itsCenter[0]);
  }
  // The phase center must be given in J2000 as two values (ra,dec).
  // In the future time dependent frames can be done as in UVWFlagger.
  if (itsCenter.size() != 2)
    throw Exception("2 values must be given in PhaseShift phasecenter");
  /// ASSERTSTR (itsCenter.size() < 4,
  ///"Up to 3 values can be given in PhaseShift phasecenter");
  casacore::MDirection phaseCenter;
  if (itsCenter.size() == 1) {
    string str = boost::to_upper_copy(itsCenter[0]);
    MDirection::Types tp;
    if (!MDirection::getType(tp, str))
      throw Exception(str +
                      " is an invalid source type"
                      " in PhaseShift phasecenter");
    return MDirection(tp);
  }
  casacore::Quantity q0, q1;
  if (!casacore::MVAngle::read(q0, itsCenter[0]))
    throw Exception(itsCenter[0] +
                    " is an invalid RA or longitude"
                    " in PhaseShift phasecenter");
  if (!casacore::MVAngle::read(q1, itsCenter[1]))
    throw Exception(itsCenter[1] +
                    " is an invalid DEC or latitude"
                    " in PhaseShift phasecenter");
  MDirection::Types type = MDirection::J2000;
  if (itsCenter.size() > 2) {
    string str = boost::to_upper_copy(itsCenter[2]);
    MDirection::Types tp;
    if (!MDirection::getType(tp, str))
      throw Exception(str +
                      " is an invalid direction type"
                      " in PhaseShift phasecenter");
  }
  return MDirection(q0, q1, type);
}

void PhaseShift::fillTransMatrix(Matrix<double>& mat,
                                 const base::Direction& direction) {
  double sinra = sin(direction.ra);
  double cosra = cos(direction.ra);
  double sindec = sin(direction.dec);
  double cosdec = cos(direction.dec);
  mat(0, 0) = cosra;
  mat(1, 0) = -sinra;
  mat(2, 0) = 0;
  mat(0, 1) = -sinra * sindec;
  mat(1, 1) = -cosra * sindec;
  mat(2, 1) = cosdec;
  mat(0, 2) = sinra * cosdec;
  mat(1, 2) = cosra * cosdec;
  mat(2, 2) = sindec;
}

}  // namespace steps
}  // namespace dp3
