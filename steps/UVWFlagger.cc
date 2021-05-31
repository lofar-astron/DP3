// UVWFlagger.cc: DPPP step class to flag data on UVW coordinates
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "UVWFlagger.h"

#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/Exceptions.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Utilities/GenSort.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/make_unique.hpp>

#include <algorithm>
#include <iostream>

using casacore::Cube;
using casacore::IPosition;
using casacore::Matrix;
using casacore::MDirection;
using casacore::MVAngle;
using casacore::Quantity;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

UVWFlagger::UVWFlagger(InputStep* input, const common::ParameterSet& parset,
                       const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsNTimes(0),
      itsRecWavel(),

      itsRangeUVm(fillUVW(parset, prefix, "uvm", true)),
      itsRangeUm(fillUVW(parset, prefix, "um", false)),
      itsRangeVm(fillUVW(parset, prefix, "vm", false)),
      itsRangeWm(fillUVW(parset, prefix, "wm", false)),
      itsRangeUVl(fillUVW(parset, prefix, "uvlambda", false)),
      itsRangeUl(fillUVW(parset, prefix, "ulambda", false)),
      itsRangeVl(fillUVW(parset, prefix, "vlambda", false)),
      itsRangeWl(fillUVW(parset, prefix, "wlambda", false)),
      itsIsDegenerate(itsRangeUVm.size() + itsRangeUVl.size() +
                          itsRangeUm.size() + itsRangeVm.size() +
                          itsRangeWm.size() + itsRangeUl.size() +
                          itsRangeVl.size() + itsRangeWl.size() ==
                      0),
      itsUVWCalc(),
      itsCenter(parset.getStringVector(prefix + "phasecenter",
                                       std::vector<string>())),
      itsTimer(),
      itsUVWTimer(),
      itsFlagCounter(input->msName(), parset, prefix + "count.") {}

UVWFlagger::~UVWFlagger() {}

void UVWFlagger::show(std::ostream& os) const {
  if (itsIsDegenerate) {
    return;
  }

  os << "UVWFlagger " << itsName << '\n';
  std::vector<double> uvm(itsRangeUVm);
  for (unsigned int i = 0; i < uvm.size(); ++i) {
    if (uvm[i] > 0) {
      uvm[i] = sqrt(uvm[i]);
    }
  }
  os << "  uvm:            " << uvm << '\n';
  os << "  um:             " << itsRangeUm << '\n';
  os << "  vm:             " << itsRangeVm << '\n';
  os << "  wm:             " << itsRangeWm << '\n';
  os << "  uvlambda:       " << itsRangeUVl << '\n';
  os << "  ulambda:        " << itsRangeUl << '\n';
  os << "  vlambda:        " << itsRangeVl << '\n';
  os << "  wlambda:        " << itsRangeWl << '\n';
  os << "  phasecenter:    " << itsCenter << '\n';
}

void UVWFlagger::showCounts(std::ostream& os) const {
  if (itsIsDegenerate) {
    return;
  }
  os << '\n' << "Flags set by UVWFlagger " << itsName;
  os << '\n' << "=======================" << '\n';
  itsFlagCounter.showBaseline(os, itsNTimes);
  itsFlagCounter.showChannel(os, itsNTimes);
}

void UVWFlagger::showTimings(std::ostream& os, double duration) const {
  double flagDur = itsTimer.getElapsed();
  os << "  ";
  base::FlagCounter::showPerc1(os, flagDur, duration);
  os << " UVWFlagger " << itsName << '\n';
  if (!itsCenter.empty()) {
    os << "          ";
    base::FlagCounter::showPerc1(os, itsUVWTimer.getElapsed(), flagDur);
    os << " of it spent in calculating UVW coordinates" << '\n';
  }
}

void UVWFlagger::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setWriteFlags();
  // Convert the given frequencies to possibly averaged frequencies.
  // Divide it by speed of light to get reciproke of wavelengths.
  itsRecWavel = infoIn.chanFreqs();
  const double inv_c = 1.0 / casacore::C::c;
  for (double& wl : itsRecWavel) {
    wl *= inv_c;
  }
  // Handle the phase center (if given).
  if (!itsCenter.empty()) {
    handleCenter();
  }
  // Initialize the flag counters.
  itsFlagCounter.init(getInfo());
}

bool UVWFlagger::process(const DPBuffer& buf) {
  if (itsIsDegenerate) {
    getNextStep()->process(buf);
    return true;
  }

  itsTimer.start();
  // Because no buffers are kept, we can reference the filled arrays
  // in the input buffer instead of copying them.
  itsBuffer.referenceFilled(buf);
  Cube<bool>& flags = itsBuffer.getFlags();
  // Loop over the baselines and flag as needed.
  const IPosition& shape = flags.shape();
  unsigned int nrcorr = shape[0];
  unsigned int nrchan = shape[1];
  unsigned int nrbl = shape[2];
  unsigned int nr = nrcorr * nrchan;
  assert(nrchan == itsRecWavel.size());
  // Input uvw coordinates are only needed if no new phase center is used.
  Matrix<double> uvws;
  if (itsCenter.empty()) {
    uvws.reference(itsInput->fetchUVW(buf, itsBuffer, itsTimer));
  }
  const double* uvwPtr = uvws.data();
  bool* flagPtr = flags.data();
  const bool* origPtr = buf.getFlags().data();
  for (unsigned int i = 0; i < nrbl; ++i) {
    std::array<double, 3> uvw{0.0, 0.0, 0.0};
    if (!itsCenter.empty()) {
      // A different phase center is given, so calculate UVW for it.
      common::NSTimer::StartStop ssuvwtimer(itsUVWTimer);
      uvw = itsUVWCalc->getUVW(getInfo().getAnt1()[i], getInfo().getAnt2()[i],
                               buf.getTime());
      uvwPtr = uvw.data();
      /// cout << "uvw = " << uvw << '\n';
    }
    double uvdist = uvwPtr[0] * uvwPtr[0] + uvwPtr[1] * uvwPtr[1];
    bool flagBL = false;
    if (!itsRangeUVm.empty()) {
      // UV-distance is sqrt(u^2 + v^2).
      // The sqrt is not needed because itsRangeUVm is squared.
      flagBL = testUVWm(uvdist, itsRangeUVm);
    }
    if (!(flagBL || itsRangeUm.empty())) {
      flagBL = testUVWm(uvwPtr[0], itsRangeUm);
    }
    if (!(flagBL || itsRangeVm.empty())) {
      flagBL = testUVWm(uvwPtr[1], itsRangeVm);
    }
    if (!(flagBL || itsRangeWm.empty())) {
      flagBL = testUVWm(uvwPtr[2], itsRangeWm);
    }
    if (flagBL) {
      // Flag entire baseline.
      std::fill(flagPtr, flagPtr + nr, true);
    } else {
      if (!itsRangeUVl.empty()) {
        // UV-distance is sqrt(u^2 + v^2).
        testUVWl(sqrt(uvdist), itsRangeUVl, flagPtr, nrcorr);
      }
      if (!itsRangeUl.empty()) {
        testUVWl(uvwPtr[0], itsRangeUl, flagPtr, nrcorr);
      }
      if (!itsRangeVl.empty()) {
        testUVWl(uvwPtr[1], itsRangeVl, flagPtr, nrcorr);
      }
      if (!itsRangeWl.empty()) {
        testUVWl(uvwPtr[2], itsRangeWl, flagPtr, nrcorr);
      }
    }
    // Count the flags set newly.
    for (unsigned int j = 0; j < nrchan; ++j) {
      if (*flagPtr && !*origPtr) {
        itsFlagCounter.incrBaseline(i);
        itsFlagCounter.incrChannel(j);
      }
      flagPtr += nrcorr;
      origPtr += nrcorr;
    }
    uvwPtr += 3;
  }
  // Let the next step do its processing.
  itsTimer.stop();
  itsNTimes++;
  getNextStep()->process(itsBuffer);
  return true;
}

void UVWFlagger::finish() {
  // Let the next step finish its processing.
  getNextStep()->finish();
}

bool UVWFlagger::testUVWm(double uvw, const std::vector<double>& ranges) {
  for (size_t i = 0; i < ranges.size(); i += 2) {
    if (uvw > ranges[i] && uvw < ranges[i + 1]) {
      return true;
    }
  }
  return false;
}

void UVWFlagger::testUVWl(double uvw, const std::vector<double>& ranges,
                          bool* flagPtr, unsigned int nrcorr) {
  // This loop could be made more efficient if it is guaranteed that
  // itsRecWavel is in strict ascending or descending order.
  // It is expected that the nr of ranges is so small that it is not
  // worth the trouble, but it could be done if ever needed.
  for (unsigned int j = 0; j < itsRecWavel.size(); ++j) {
    double uvwl = uvw * itsRecWavel[j];
    for (size_t i = 0; i < ranges.size(); i += 2) {
      if (uvwl > ranges[i] && uvwl < ranges[i + 1]) {
        std::fill(flagPtr, flagPtr + nrcorr, true);
        break;
      }
    }
    flagPtr += nrcorr;
  }
}

std::vector<double> UVWFlagger::fillUVW(const common::ParameterSet& parset,
                                        const string& prefix,
                                        const string& name, bool square) {
  // Get possible range, minimum, and maximum.
  std::vector<string> uvs =
      parset.getStringVector(prefix + name + "range", std::vector<string>());
  double minuv = parset.getDouble(prefix + name + "min", 0.);
  double maxuv = parset.getDouble(prefix + name + "max", 0.);
  // Process the ranges.
  std::vector<double> vals;
  vals.reserve(2 * uvs.size());
  for (std::vector<string>::const_iterator str = uvs.begin(); str != uvs.end();
       ++str) {
    // Each range can be given as st..end or val+-halfwidth.
    // Find the .. or +- token.
    bool usepm = false;
    string::size_type pos;
    pos = str->find("..");
    if (pos == string::npos) {
      usepm = true;
      pos = str->find("+-");
      if (pos == string::npos)
        throw Exception("UVWFlagger " + name + "range '" + *str +
                        "' should be range using .. or +-");
    }
    string str1 = str->substr(0, pos);
    string str2 = str->substr(pos + 2);
    double v1 = common::strToDouble(str1);
    double v2 = common::strToDouble(str2);
    if (usepm) {
      double hw = v2;
      v2 = v1 + hw;
      v1 -= hw;
    }
    vals.push_back(v1);
    vals.push_back(v2);
  }
  // If minimum or maximum is given, add them as a range as well.
  if (minuv > 0) {
    vals.push_back(-1e15);
    vals.push_back(minuv);
  }
  if (maxuv > 0) {
    vals.push_back(maxuv);
    vals.push_back(1e15);
  }
  if (square) {
    for (std::vector<double>::iterator iter = vals.begin(); iter != vals.end();
         ++iter) {
      if (*iter != -1e15) {
        *iter = *iter * *iter;
      }
    }
  }
  return vals;
}

void UVWFlagger::handleCenter() {
  // The phase center can be given as one, two, or three values.
  // I.e., as source name, ra,dec or ra,dec,frame.
  if (itsCenter.size() >= 4)
    throw Exception("Up to 3 values can be given in UVWFlagger phasecenter");
  MDirection phaseCenter;
  if (itsCenter.size() == 1) {
    string str = boost::to_upper_copy(itsCenter[0]);
    MDirection::Types tp;
    if (!MDirection::getType(tp, str))
      throw Exception(str +
                      " is an invalid source type"
                      " in UVWFlagger phasecenter");
    phaseCenter = MDirection(tp);
  } else {
    Quantity q0, q1;
    if (!MVAngle::read(q0, itsCenter[0]))
      throw Exception(itsCenter[0] +
                      " is an invalid RA or longitude"
                      " in UVWFlagger phasecenter");
    if (!MVAngle::read(q1, itsCenter[1]))
      throw Exception(itsCenter[1] +
                      " is an invalid DEC or latitude"
                      " in UVWFlagger phasecenter");
    MDirection::Types type = MDirection::J2000;
    if (itsCenter.size() > 2) {
      string str = boost::to_upper_copy(itsCenter[2]);
      MDirection::Types tp;
      if (!MDirection::getType(tp, str))
        throw Exception(str +
                        " is an invalid direction type in UVWFlagger"
                        " in UVWFlagger phasecenter");
    }
    phaseCenter = MDirection(q0, q1, type);
  }
  // Create the UVW calculator.
  itsUVWCalc = boost::make_unique<base::UVWCalculator>(
      phaseCenter, getInfo().arrayPos(), getInfo().antennaPos());
}

}  // namespace steps
}  // namespace dp3
