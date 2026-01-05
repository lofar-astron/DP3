// UVWFlagger.cc: DP3 step class to flag data on UVW coordinates
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "UVWFlagger.h"

#include "base/DPInfo.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Utilities/GenSort.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <algorithm>
#include <array>
#include <iostream>

using dp3::base::BdaBuffer;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

UVWFlagger::UVWFlagger(const common::ParameterSet& parset,
                       const std::string& prefix, MsType inputType)
    : itsInputType(inputType),
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
                                       std::vector<std::string>())),
      itsTimer(),
      itsUVWTimer(),
      itsFlagCounter(parset, prefix + "count.") {}

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
  if (itsIsDegenerate) {
    return;
  }

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
  Step::updateInfo(infoIn);
  // Convert the given frequencies to possibly averaged frequencies.
  // Divide it by speed of light to get reciprocal of wavelengths.
  itsRecWavel = infoIn.BdaChanFreqs();
  const double inv_c = 1.0 / casacore::C::c;

  for (std::vector<double>& baseline_channel_frequencies : itsRecWavel) {
    for (double& wavelength : baseline_channel_frequencies) {
      wavelength *= inv_c;
    }
  }
  // Handle the phase center (if given).
  if (!itsCenter.empty()) {
    handleCenter();
  }
  // Initialize the flag counters.
  itsFlagCounter.init(getInfoOut());
}

bool UVWFlagger::process(std::unique_ptr<base::DPBuffer> buffer) {
  if (itsIsDegenerate) {
    getNextStep()->process(std::move(buffer));
    return true;
  }

  itsTimer.start();
  // Loop over the baselines and flag as needed.
  unsigned int n_baselines = buffer->GetFlags().shape(0);
  unsigned int n_channels = buffer->GetFlags().shape(1);
  unsigned int n_correlations = buffer->GetFlags().shape(2);
  assert(n_channels == itsRecWavel[0].size());
  // Input uvw coordinates are only needed if no new phase center is used.
  const double* uvwPtr = nullptr;
  if (itsCenter.empty()) {
    assert(buffer->GetUvw().size() != 0);
    uvwPtr = buffer->GetUvw().data();
  }
  bool* flagPtr = buffer->GetFlags().data();

  // A std::vector<bool> doesn't store a contiguous array of bools, instead it
  // stores 8 bools per byte. Since the code wants to use an array of bools use
  // a std::unique_ptr instead.
  const size_t stride = n_correlations * n_channels;
  std::unique_ptr<bool[]> F{new bool[stride]};
  for (unsigned int i = 0; i < n_baselines; ++i) {
    std::array<double, 3> uvw;
    if (itsCenter.empty()) {
      std::copy_n(&uvwPtr[3 * i], 3, uvw.data());
    } else {
      // A different phase center is given, so calculate UVW for it.
      common::NSTimer::StartStop ssuvwtimer(itsUVWTimer);
      uvw = itsUVWCalc->getUVW(getInfoOut().getAnt1()[i],
                               getInfoOut().getAnt2()[i], buffer->GetTime());
    }

    // Copy the original flags so it's possible to count the number of newly
    // set flags.
    bool* origPtr = F.get();
    std::copy_n(flagPtr, stride, origPtr);
    doFlag(uvw, flagPtr, n_correlations, n_channels);

    // Count the flags set newly.
    for (unsigned int j = 0; j < n_channels; ++j) {
      if (*flagPtr && !*origPtr) {
        itsFlagCounter.incrBaseline(i);
        itsFlagCounter.incrChannel(j);
      }
      flagPtr += n_correlations;
      origPtr += n_correlations;
    }
  }
  // Let the next step do its processing.
  itsTimer.stop();
  itsNTimes++;
  getNextStep()->process(std::move(buffer));
  return true;
}

bool UVWFlagger::process(std::unique_ptr<BdaBuffer> buffer) {
  if (itsIsDegenerate) {
    getNextStep()->process(std::move(buffer));
    return true;
  }
  itsTimer.start();

  std::vector<BdaBuffer::Row> rows = buffer->GetRows();
  for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
    BdaBuffer::Row& row = rows[row_index];

    // Loop over the baselines and flag as needed.
    assert(row.n_channels == itsRecWavel[row.baseline_nr].size());

    // Input uvw coordinates are only needed if no new phase center is used.
    std::array<double, 3> uvw;
    if (itsCenter.empty()) {
      std::copy_n(row.uvw, 3, uvw.data());
    } else {
      // A different phase center is given, so calculate UVW for it.
      common::NSTimer::StartStop ssuvwtimer(itsUVWTimer);
      uvw =
          itsUVWCalc->getUVW(getInfoOut().getAnt1()[row.baseline_nr],
                             getInfoOut().getAnt2()[row.baseline_nr], row.time);
    }

    bool* row_flags_pointer = buffer->GetFlags(row_index);

    // Copy original flags to calculate flagged visibilities
    // Use aocommon's UVector since std::vector<bool>::data() is inacessible.
    const aocommon::UVector<bool> original_flags(
        row_flags_pointer, row_flags_pointer + row.GetDataSize());
    const bool* original_flags_pointer = original_flags.data();

    doFlag(uvw, row_flags_pointer, row.n_correlations, row.n_channels,
           row.baseline_nr);
    // Count the flags set newly.
    for (std::size_t channel = 0; channel < row.n_channels; ++channel) {
      if (*row_flags_pointer && !*original_flags_pointer) {
        itsFlagCounter.incrBaseline(row.baseline_nr);
        itsFlagCounter.incrChannel(channel);
      }
      row_flags_pointer += row.n_correlations;
      original_flags_pointer += row.n_correlations;
    }
  }

  itsTimer.stop();
  itsNTimes++;
  getNextStep()->process(std::move(buffer));
  return true;
}

void UVWFlagger::doFlag(const std::array<double, 3>& uvw, bool* flagPtr,
                        unsigned int n_correlations, unsigned int n_channels,
                        unsigned int baseline_id) {
  double uvdist = uvw[0] * uvw[0] + uvw[1] * uvw[1];
  bool flagBL = false;
  if (!itsRangeUVm.empty()) {
    // UV-distance is sqrt(u^2 + v^2).
    // The sqrt is not needed because itsRangeUVm is squared.
    flagBL = testUVWm(uvdist, itsRangeUVm);
  }
  if (!(flagBL || itsRangeUm.empty())) {
    flagBL = testUVWm(uvw[0], itsRangeUm);
  }
  if (!(flagBL || itsRangeVm.empty())) {
    flagBL = testUVWm(uvw[1], itsRangeVm);
  }
  if (!(flagBL || itsRangeWm.empty())) {
    flagBL = testUVWm(uvw[2], itsRangeWm);
  }
  if (flagBL) {
    // Flag entire baseline.
    std::fill_n(flagPtr, n_correlations * n_channels, true);
  } else {
    if (!itsRangeUVl.empty()) {
      // UV-distance is sqrt(u^2 + v^2).
      testUVWl(sqrt(uvdist), itsRangeUVl, flagPtr, n_correlations, baseline_id);
    }
    if (!itsRangeUl.empty()) {
      testUVWl(uvw[0], itsRangeUl, flagPtr, n_correlations, baseline_id);
    }
    if (!itsRangeVl.empty()) {
      testUVWl(uvw[1], itsRangeVl, flagPtr, n_correlations, baseline_id);
    }
    if (!itsRangeWl.empty()) {
      testUVWl(uvw[2], itsRangeWl, flagPtr, n_correlations, baseline_id);
    }
  }
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
                          bool* flagPtr, unsigned int n_correlations,
                          unsigned int baseline_id) {
  // This loop could be made more efficient if it is guaranteed that
  // itsRecWavel is in strict ascending or descending order.
  // It is expected that the nr of ranges is so small that it is not
  // worth the trouble, but it could be done if ever needed.
  for (unsigned int j = 0; j < itsRecWavel[baseline_id].size(); ++j) {
    double uvwl = uvw * itsRecWavel[baseline_id][j];
    for (size_t i = 0; i < ranges.size(); i += 2) {
      if (uvwl > ranges[i] && uvwl < ranges[i + 1]) {
        std::fill_n(flagPtr, n_correlations, true);
        break;
      }
    }
    flagPtr += n_correlations;
  }
}

std::vector<double> UVWFlagger::fillUVW(const common::ParameterSet& parset,
                                        const std::string& prefix,
                                        const std::string& name, bool square) {
  // Get possible range, minimum, and maximum.
  std::vector<std::string> uvs = parset.getStringVector(
      prefix + name + "range", std::vector<std::string>());
  double minuv = parset.getDouble(prefix + name + "min", 0.);
  double maxuv = parset.getDouble(prefix + name + "max", 0.);
  // Process the ranges.
  std::vector<double> vals;
  vals.reserve(2 * uvs.size());
  for (std::vector<std::string>::const_iterator str = uvs.begin();
       str != uvs.end(); ++str) {
    // Each range can be given as st..end or val+-halfwidth.
    // Find the .. or +- token.
    bool usepm = false;
    std::string::size_type pos;
    pos = str->find("..");
    if (pos == std::string::npos) {
      usepm = true;
      pos = str->find("+-");
      if (pos == std::string::npos)
        throw std::runtime_error("UVWFlagger " + name + "range '" + *str +
                                 "' should be range using .. or +-");
    }
    std::string str1 = str->substr(0, pos);
    std::string str2 = str->substr(pos + 2);
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
    throw std::runtime_error(
        "Up to 3 values can be given in UVWFlagger phasecenter");
  casacore::MDirection phaseCenter;
  if (itsCenter.size() == 1) {
    std::string str = boost::to_upper_copy(itsCenter[0]);
    casacore::MDirection::Types tp;
    if (!casacore::MDirection::getType(tp, str))
      throw std::runtime_error(str +
                               " is an invalid source type"
                               " in UVWFlagger phasecenter");
    phaseCenter = casacore::MDirection(tp);
  } else {
    casacore::Quantity q0;
    casacore::Quantity q1;
    if (!casacore::MVAngle::read(q0, itsCenter[0]))
      throw std::runtime_error(itsCenter[0] +
                               " is an invalid RA or longitude"
                               " in UVWFlagger phasecenter");
    if (!casacore::MVAngle::read(q1, itsCenter[1]))
      throw std::runtime_error(itsCenter[1] +
                               " is an invalid DEC or latitude"
                               " in UVWFlagger phasecenter");
    casacore::MDirection::Types type = casacore::MDirection::J2000;
    if (itsCenter.size() > 2) {
      std::string str = boost::to_upper_copy(itsCenter[2]);
      casacore::MDirection::Types tp;
      if (!casacore::MDirection::getType(tp, str))
        throw std::runtime_error(str +
                                 " is an invalid direction type in UVWFlagger"
                                 " in UVWFlagger phasecenter");
    }
    phaseCenter = casacore::MDirection(q0, q1, type);
  }
  // Create the UVW calculator.
  itsUVWCalc = std::make_unique<base::UVWCalculator>(
      phaseCenter, getInfoOut().arrayPos(), getInfoOut().antennaPos());
}

}  // namespace steps
}  // namespace dp3
