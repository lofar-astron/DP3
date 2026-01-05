// PreFlagger.cc: DPPP step class to (un)flag data on channel, baseline, time
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "PreFlagger.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <iostream>
#include <stack>

#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/RecordGram.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/BasicMath/Functors.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCEpoch.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <xtensor/xcomplex.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xview.hpp>

#include "base/DPBuffer.h"
#include "base/DPInfo.h"

#include "../common/ParameterSet.h"
#include "../common/StreamUtil.h"

using casacore::Block;
using casacore::IPosition;
using casacore::MDirection;
using casacore::MeasFrame;
using casacore::MEpoch;
using casacore::MVAngle;
using casacore::MVDirection;
using casacore::MVEpoch;
using casacore::MVTime;
using casacore::Quantity;
using casacore::Record;
using casacore::RecordGram;
using casacore::TableExprNode;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

PreFlagger::PreFlagger(const common::ParameterSet& parset,
                       const std::string& prefix)
    : itsName(prefix),
      itsMode(Mode::kSetFlag),
      itsPSet(parset, prefix),
      itsCount(0),
      itsFlagCounter(parset, prefix + "count.") {
  string mode = boost::to_lower_copy(parset.getString(prefix + "mode", "set"));
  if (mode == "clear") {
    itsMode = Mode::kClearFlag;
  } else if (mode == "setcomplement" || mode == "setother") {
    itsMode = Mode::kSetComp;
  } else if (mode == "clearcomplement" || mode == "clearother") {
    itsMode = Mode::kClearComp;
  } else {
    if (mode != "set")
      throw std::runtime_error(
          "invalid preflagger mode: "
          "only set, clear, and set/clearcomplement/other are valid");
  }
}

PreFlagger::~PreFlagger() {}

void PreFlagger::show(std::ostream& os) const {
  os << "PreFlagger " << itsName << '\n';
  os << "  mode:           ";
  switch (itsMode) {
    case Mode::kSetFlag:
      os << "set";
      break;
    case Mode::kClearFlag:
      os << "clear";
      break;
    case Mode::kSetComp:
      os << "setcomplement";
      break;
    case Mode::kClearComp:
      os << "clearcomplement";
      break;
  }
  os << '\n';
  itsPSet.show(os, false);
}

void PreFlagger::showCounts(std::ostream& os) const {
  os << '\n' << "Flags set by PreFlagger " << itsName;
  os << '\n' << "=======================" << '\n';
  itsFlagCounter.showBaseline(os, itsCount);
  itsFlagCounter.showChannel(os, itsCount);
}

void PreFlagger::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " PreFlagger " << itsName << '\n';
}

void PreFlagger::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);
  itsPSet.updateInfo(info_in);
  itsFlagCounter.init(info_in);
}

bool PreFlagger::process(std::unique_ptr<DPBuffer> buffer) {
  itsTimer.start();
  assert(buffer->GetFlags().size() != 0);
  // Do the PSet steps and combine the result with the current flags.
  // Only count if the flag changes.
  const xt::xtensor<uint8_t, 2>* const flags =
      itsPSet.process(*buffer, itsCount, xt::xtensor<bool, 1>(), itsTimer);
  switch (itsMode) {
    case Mode::kSetFlag:
      setFlags(*flags, buffer->GetFlags(), true);
      break;
    case Mode::kClearFlag:
      clearFlags(*flags, buffer->GetFlags(), true, buffer->GetData(),
                 buffer->GetWeights());
      break;
    case Mode::kSetComp:
      setFlags(*flags, buffer->GetFlags(), false);
      break;
    case Mode::kClearComp:
      clearFlags(*flags, buffer->GetFlags(), false, buffer->GetData(),
                 buffer->GetWeights());
      break;
  }
  itsTimer.stop();
  // Let the next step do its processing.
  getNextStep()->process(std::move(buffer));
  itsCount++;
  return true;
}

void PreFlagger::setFlags(const xt::xtensor<uint8_t, 2>& in,
                          base::DPBuffer::FlagsType& out, bool mode) {
  const size_t n_baselines = out.shape(0);
  const size_t n_channels = out.shape(1);
  assert(in.shape() == xt::shape({n_baselines, n_channels}));

  for (std::size_t baseline = 0; baseline < n_baselines; ++baseline) {
    for (std::size_t channel = 0; channel < n_channels; ++channel) {
      if (in(baseline, channel) == mode && !out(baseline, channel, 0)) {
        // Only 1st corr is counted.
        itsFlagCounter.incrBaseline(baseline);
        itsFlagCounter.incrChannel(channel);
        xt::view(out, baseline, channel, xt::all()).fill(true);
      }
    }
  }
}

void PreFlagger::clearFlags(const xt::xtensor<uint8_t, 2>& in,
                            base::DPBuffer::FlagsType& out, bool mode,
                            const base::DPBuffer::DataType& data,
                            const base::DPBuffer::WeightsType& weights) {
  const size_t n_baselines = out.shape(0);
  const size_t n_channels = out.shape(1);
  assert(in.shape() == xt::shape({n_baselines, n_channels}));
  assert(out.shape() == data.shape());
  assert(out.shape() == weights.shape());

  for (std::size_t baseline = 0; baseline < n_baselines; ++baseline) {
    for (std::size_t channel = 0; channel < n_channels; ++channel) {
      if (in(baseline, channel) == mode) {
        // Flags for invalid data are not cleared.
        auto data_view = xt::view(data, baseline, channel, xt::all());
        auto weights_view = xt::view(weights, baseline, channel, xt::all());
        const bool flag =
            xt::any(!xt::isfinite(data_view) || xt::equal(weights_view, 0.0f));

        if (out(baseline, channel, 0) != flag) {
          itsFlagCounter.incrBaseline(baseline);
          itsFlagCounter.incrChannel(channel);
          xt::view(out, baseline, channel, xt::all()).fill(flag);
        }
      }
    }
  }
}

void PreFlagger::finish() {
  // Let the next step finish its processing.
  getNextStep()->finish();
}

PreFlagger::PSet::PSet(const common::ParameterSet& parset,
                       const std::string& prefix)
    : itsName(prefix),
      itsFlagOnUV(false),
      itsFlagOnBL(false),
      itsFlagOnAmpl(false),
      itsFlagOnPhase(false),
      itsFlagOnReal(false),
      itsFlagOnImag(false),
      itsFlagOnAzEl(false),
      itsSelBL(parset, prefix, true) {
  // Read all possible parameters.
  itsStrTime =
      parset.getStringVector(prefix + "timeofday", std::vector<std::string>());
  itsStrLST =
      parset.getStringVector(prefix + "lst", std::vector<std::string>());
  itsStrATime =
      parset.getStringVector(prefix + "abstime", std::vector<std::string>());
  itsStrRTime =
      parset.getStringVector(prefix + "reltime", std::vector<std::string>());
  itsTimeSlot = parset.getUintVector(
      prefix + "timeslot", std::vector<unsigned int>(), true);  // expand ..
  itsStrAzim =
      parset.getStringVector(prefix + "azimuth", std::vector<std::string>());
  itsStrElev =
      parset.getStringVector(prefix + "elevation", std::vector<std::string>());
  itsMinUV = parset.getDouble(prefix + "uvmmin", -1);
  itsMaxUV = parset.getDouble(prefix + "uvmmax", -1);
  itsFlagOnUV = (itsMinUV >= 0) || (itsMaxUV > 0);
  itsStrFreq =
      parset.getStringVector(prefix + "freqrange", std::vector<std::string>());
  itsStrChan =
      parset.getStringVector(prefix + "chan", std::vector<std::string>());
  itsAmplMin = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "amplmin", std::string())),
                                -1e30, itsFlagOnAmpl);
  itsAmplMax = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "amplmax", std::string())),
                                1e30, itsFlagOnAmpl);
  itsPhaseMin = fillValuePerCorr(common::ParameterValue(parset.getString(
                                     prefix + "phasemin", std::string())),
                                 -1e30, itsFlagOnPhase);
  itsPhaseMax = fillValuePerCorr(common::ParameterValue(parset.getString(
                                     prefix + "phasemax", std::string())),
                                 1e30, itsFlagOnPhase);
  itsRealMin = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "realmin", std::string())),
                                -1e30, itsFlagOnReal);
  itsRealMax = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "realmax", std::string())),
                                1e30, itsFlagOnReal);
  itsImagMin = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "imagmin", std::string())),
                                -1e30, itsFlagOnImag);
  itsImagMax = fillValuePerCorr(common::ParameterValue(parset.getString(
                                    prefix + "imagmax", std::string())),
                                1e30, itsFlagOnImag);
  itsStrExpr = parset.getString(prefix + "expr", std::string());
  // Parse the possible pset expression and convert to RPN form.
  if (!itsStrExpr.empty()) {
    std::vector<std::string> names = exprToRpn(itsStrExpr);
    // Create PSet objects for all operands.
    itsPSets.reserve(names.size());
    for (unsigned int i = 0; i < names.size(); ++i) {
      itsPSets.push_back(
          std::make_shared<PSet>(parset, prefix + names[i] + '.'));
    }
  }
}

void PreFlagger::PSet::updateInfo(const DPInfo& info) {
  itsInfo = &info;
  // Fill the matrix with the baselines to flag.
  fillBLMatrix();
  // Handle the possible date/time parameters.
  itsTimes = fillTimes(itsStrTime, true, true);
  itsLST = fillTimes(itsStrLST, true, true);
  itsATimes = fillTimes(itsStrATime, false, false);
  itsRTimes = fillTimes(itsStrRTime, true, false);
  itsFlagOnTime = !(itsTimeSlot.empty() && itsTimes.empty() && itsLST.empty() &&
                    itsATimes.empty() && itsRTimes.empty());
  // Handle possible azimuth/elevation ranges.
  itsAzimuth = fillTimes(itsStrAzim, true, true);
  itsElevation = fillTimes(itsStrElev, true, true);
  itsFlagOnAzEl = !(itsAzimuth.empty() && itsElevation.empty());
  // Determine if to flag on UV distance.
  // If so, square the distances to avoid having to take the sqrt in flagUV.
  if (itsMinUV >= 0) {
    itsMinUV *= itsMinUV;
  }
  if (itsMaxUV > 0) {
    itsMaxUV *= itsMaxUV;
  } else {
    // Make it a very high number.
    itsMaxUV = 1e30;
  }
  if (itsMinUV >= itsMaxUV)
    throw std::runtime_error("PreFlagger uvmmin should be < uvmmax");
  // Determine if only flagging on time info is done.
  itsFlagOnTimeOnly =
      (!(itsFlagOnUV || itsFlagOnBL || itsFlagOnAzEl || itsFlagOnAmpl ||
         itsFlagOnPhase || itsFlagOnReal || itsFlagOnImag) &&
       itsPSets.empty());
  // Size the object's buffers (used in process) correctly.
  const size_t nrchan = info.nchan();
  const size_t nrbaselines = info.nbaselines();
  itsFlags.resize({nrbaselines, nrchan});
  itsMatchBL.resize({nrbaselines});
  // Determine the channels to be flagged.
  if (!(itsStrChan.empty() && itsStrFreq.empty())) {
    fillChannels(info);
    if (!itsChannels.empty()) {
      itsFlagOnTimeOnly = false;
    }
  }
  // Do the same for the child steps.
  for (unsigned int i = 0; i < itsPSets.size(); ++i) {
    itsPSets[i]->updateInfo(info);
  }
}

void PreFlagger::PSet::fillChannels(const DPInfo& info) {
  unsigned int nrchan = info.nchan();
  xt::xtensor<int, 1> selChan(std::array<size_t, 1>{nrchan});
  if (itsStrChan.empty()) {
    selChan.fill(true);
  } else {
    // Set selChan for channels not exceeding nr of channels.
    selChan.fill(false);
    Record rec;
    rec.define("nchan", nrchan);
    double result;
    for (unsigned int i = 0; i < itsStrChan.size(); ++i) {
      // Evaluate possible expressions.
      // Split the value if start..end is given.
      unsigned int startch, endch;
      std::string::size_type pos = itsStrChan[i].find("..");
      if (pos == std::string::npos) {
        TableExprNode node(RecordGram::parse(rec, itsStrChan[i]));
        node.get(rec, result);
        startch = (unsigned int)(result + 0.001);
        endch = startch;
      } else {
        if (pos == 0 || pos >= itsStrChan[i].size() - 2)
          throw std::runtime_error(
              "No start or end given in PreFlagger channel range " +
              itsStrChan[i]);
        TableExprNode node1(
            RecordGram::parse(rec, itsStrChan[i].substr(0, pos)));
        node1.get(rec, result);
        startch = (unsigned int)(result + 0.001);
        TableExprNode node2(
            RecordGram::parse(rec, itsStrChan[i].substr(pos + 2)));
        node2.get(rec, result);
        endch = (unsigned int)(result + 0.001);
        if (startch > endch)
          throw std::runtime_error("Start " + std::to_string(startch) +
                                   " must be <= end " + std::to_string(endch) +
                                   " in PreFlagger channel range " +
                                   itsStrChan[i]);
      }
      if (startch < nrchan) {
        xt::view(selChan, xt::range(startch, std::min(endch + 1, nrchan)))
            .fill(true);
      }
    }
  }
  // Now determine which channels to use from given frequency ranges.
  // AND it with the channel selection given above.
  if (!itsStrFreq.empty()) {
    selChan &= handleFreqRanges(itsInfo->chanFreqs());
  }
  // Turn the channels into a mask.
  itsChannels.clear();
  itsChanFlags.resize({nrchan});
  itsChanFlags.fill(false);
  for (unsigned int channel = 0; channel < nrchan; ++channel) {
    if (selChan[channel]) {
      itsChannels.push_back(channel);
      itsChanFlags[channel] = true;
    }
  }
}

void PreFlagger::PSet::show(std::ostream& os, bool showName) const {
  if (showName) {
    os << "  pset " << itsName << '\n';
  }
  if (!itsStrExpr.empty()) {
    os << "   expr:          " << itsStrExpr << '\n';
  }
  if (!itsStrLST.empty()) {
    os << "   lst:           " << itsStrLST << '\n';
  }
  if (!itsStrTime.empty()) {
    os << "   timeofday:     " << itsStrTime << '\n';
  }
  if (!itsStrATime.empty()) {
    os << "   abstime:       " << itsStrATime << '\n';
  }
  if (!itsStrRTime.empty()) {
    os << "   reltime:       " << itsStrRTime << '\n';
  }
  if (!itsTimeSlot.empty()) {
    os << "   timeslot:      " << itsTimeSlot << '\n';
  }
  if (itsFlagOnBL) {
    itsSelBL.show(os);
  }
  if (itsFlagOnUV) {
    if (itsMinUV >= 0) {
      os << "   uvmmin:        " << sqrt(itsMinUV) << '\n';
    } else {
      os << "   uvmmin:        " << itsMinUV << '\n';
    }
    os << "   uvmmax:        " << sqrt(itsMaxUV) << '\n';
  }
  if (itsFlagOnAzEl) {
    os << "   azimuth:       " << itsStrAzim << '\n';
    os << "   elevation:     " << itsStrElev << '\n';
  }
  if (!itsChannels.empty()) {
    os << "   channel:       " << itsStrChan << '\n';
    os << "   freqrange:     " << itsStrFreq << '\n';
    os << "    chan to flag: " << itsChannels << '\n';
  }
  if (itsFlagOnAmpl) {
    os << "   amplmin:       " << itsAmplMin << '\n';
    os << "   amplmax:       " << itsAmplMax << '\n';
  }
  if (itsFlagOnPhase) {
    os << "   phasemin:      " << itsPhaseMin << '\n';
    os << "   phasemax:      " << itsPhaseMax << '\n';
  }
  if (itsFlagOnReal) {
    os << "   realmin:       " << itsRealMin << '\n';
    os << "   realmax:       " << itsRealMax << '\n';
  }
  if (itsFlagOnImag) {
    os << "   imagmin:       " << itsImagMin << '\n';
    os << "   imagmax:       " << itsImagMax << '\n';
  }
  // Do it for the child steps.
  for (unsigned int i = 0; i < itsPSets.size(); ++i) {
    itsPSets[i]->show(os, true);
  }
}

xt::xtensor<uint8_t, 2>* PreFlagger::PSet::process(
    DPBuffer& out, unsigned int timeSlot, const xt::xtensor<bool, 1>& matchBL,
    common::NSTimer& timer) {
  // No need to process it if the time mismatches or if only time selection.
  if (itsFlagOnTime) {
    if (!matchTime(out.GetTime(), timeSlot)) {
      itsFlags.fill(false);
      return &itsFlags;
    }
  }
  if (itsFlagOnTimeOnly) {
    itsFlags.fill(itsFlagOnTime);
    return &itsFlags;
  }
  // Initialize the flags.
  itsFlags.fill(false);
  const size_t n_channels = out.GetFlags().shape(1);
  // Take over the baseline info from the parent. Default is all.
  if (matchBL.size() == 0) {
    itsMatchBL.fill(true);
  } else {
    itsMatchBL = matchBL;
  }
  // The PSet tree is a combination of ORs and ANDs.
  // Depth is AND, breadth is OR.
  // In a PSet flagging is done in two stages.
  // First it is determined which baselines are not flagged. It is kept
  // in the itsMatchBL block.
  // This is passed to the children who do their flagging. In this way
  // a child can minimize the amount of work to do.
  // The resulting flags of the children are handled according to the
  // operators in the RPN list.

  // First flag on baseline if necessary. Stop if no matches.
  if (itsFlagOnBL && !flagBL()) {
    return &itsFlags;
  }
  // Flag on UV distance if necessary.
  if (itsFlagOnUV && !flagUV(out.GetUvw())) {
    return &itsFlags;
  }
  // Flag on AzEl is necessary.
  if (itsFlagOnAzEl && !flagAzEl(out.GetTime())) {
    return &itsFlags;
  }
  // Convert each baseline flag to a flag per correlation/channel.
  uint8_t* flagPtr = itsFlags.data();
  for (unsigned int i = 0; i < itsMatchBL.size(); ++i) {
    if (itsMatchBL[i]) {
      std::fill(flagPtr, flagPtr + n_channels, itsMatchBL[i]);
    }
    flagPtr += n_channels;
  }
  // Flag on channel if necessary.
  if (!itsChannels.empty()) {
    flagChannels();
  }
  // Flag on amplitude, phase or real/imaginary if necessary.
  if (itsFlagOnAmpl) {
    flagAmpl(out.GetData());
  }
  if (itsFlagOnReal) {
    flagReal(out.GetData());
  }
  if (itsFlagOnImag) {
    flagImag(out.GetData());
  }
  if (itsFlagOnPhase) {
    flagPhase(out.GetData());
  }
  // Evaluate the PSet expression.
  // The expression is in RPN notation. A stack of array pointers is used
  // to keep track of intermediate results. The arrays (in the PSet objects)
  // are reused to AND or OR subexpressions. This can be done harmlessly
  // and saves the creation of too many arrays.
  if (!itsPSets.empty()) {
    std::stack<xt::xtensor<uint8_t, 2>*> results;
    for (int oper : itsRpn) {
      if (oper >= 0) {
        results.push(itsPSets[oper]->process(out, timeSlot, itsMatchBL, timer));
      } else if (oper == OpNot) {
        xt::xtensor<uint8_t, 2>* left = results.top();
        *left = !*left;
      } else if (oper == OpOr || oper == OpAnd) {
        xt::xtensor<uint8_t, 2>* right = results.top();
        results.pop();
        xt::xtensor<uint8_t, 2>* left = results.top();
        if (oper == OpOr) {
          *left |= *right;
        } else {
          *left &= *right;
        }
      } else {
        throw std::runtime_error("Expected operation NOT, OR or AND");
      }
    }
    // Finally AND the children's flags with the flags of this pset.
    if (results.size() != 1)
      throw std::runtime_error(
          "Something went wrong while evaluating expression: results.size() != "
          "1");
    xt::xtensor<uint8_t, 2>* mflags = results.top();
    itsFlags &= *mflags;
  }
  return &itsFlags;
}

bool PreFlagger::PSet::matchTime(double time, unsigned int timeSlot) const {
  if (!itsATimes.empty() && !matchRange(time, itsATimes)) {
    return false;
  }
  if (!itsRTimes.empty() &&
      !matchRange(time - itsInfo->startTime(), itsRTimes)) {
    return false;
  }
  if (!itsTimes.empty()) {
    MVTime mvtime(time / 86400);  // needs time in days
    double timeofday = time - int(mvtime.day()) * 86400.;
    if (!matchRange(timeofday, itsTimes)) {
      return false;
    }
  }
  if (!itsTimeSlot.empty()) {
    if (std::find(itsTimeSlot.begin(), itsTimeSlot.end(), timeSlot) ==
        itsTimeSlot.end()) {
      return false;
    }
  }
  if (!itsLST.empty()) {
    // Convert time from UTC to Local Apparent Sidereal Time.
    MeasFrame frame;
    frame.set(itsInfo->arrayPos());
    Quantity qtime(time, "s");
    MEpoch lst = MEpoch::Convert(MEpoch(MVEpoch(qtime), MEpoch::UTC),
                                 MEpoch::Ref(MEpoch::LAST, frame))();
    double lstSec = lst.getValue().get();       // in days
    lstSec -= int(lstSec);                      // time of day
    if (!matchRange(lstSec * 86400, itsLST)) {  // use seconds
      return false;
    }
  }
  return true;
}

bool PreFlagger::PSet::matchRange(double v,
                                  const std::vector<double>& ranges) const {
  for (unsigned int i = 0; i < ranges.size(); i += 2) {
    if (v > ranges[i] && v < ranges[i + 1]) {
      return true;
    }
  }
  return false;
}

bool PreFlagger::PSet::flagUV(const base::DPBuffer::UvwType& uvw) {
  bool match = false;
  unsigned int nrbl = itsMatchBL.size();
  for (unsigned int bl = 0; bl < nrbl; ++bl) {
    if (itsMatchBL[bl]) {
      // UV-distance is sqrt(u^2 + v^2).
      // The sqrt is not needed because minuv and maxuv are squared.
      const double uvdist = uvw(bl, 0) * uvw(bl, 0) + uvw(bl, 1) * uvw(bl, 1);
      if (uvdist >= itsMinUV && uvdist <= itsMaxUV) {
        // UV-dist mismatches, so do not flag baseline.
        itsMatchBL[bl] = false;
      } else {
        match = true;
      }
    }
  }
  return match;
}

bool PreFlagger::PSet::flagBL() {
  bool match = false;
  unsigned int nrbl = itsMatchBL.size();
  for (unsigned int i = 0; i < nrbl; ++i) {
    if (itsMatchBL[i]) {
      if (!itsFlagBL(itsInfo->getAnt1()[i], itsInfo->getAnt2()[i])) {
        // do not flag this baseline
        itsMatchBL[i] = false;
      } else {
        match = true;
      }
    }
  }
  return match;
}

bool PreFlagger::PSet::flagAzEl(double time) {
  bool match = false;
  unsigned int nrbl = itsMatchBL.size();
  // Calculate AzEl for each flagged antenna for this time slot.
  MeasFrame frame;
  Quantity qtime(time, "s");
  MEpoch epoch(MVEpoch(qtime), MEpoch::UTC);
  frame.set(epoch);
  MDirection::Convert converter(itsInfo->phaseCenter(),
                                MDirection::Ref(MDirection::AZEL, frame));
  unsigned int nrant = itsInfo->antennaNames().size();
  Block<bool> done(nrant, false);
  for (unsigned int i = 0; i < nrbl; ++i) {
    if (itsMatchBL[i]) {
      // If needed, check if ant1 matches AzEl criterium.
      // If not matching, itsMatchBL is cleared for this baseline and all
      // subsequent baselines with this antenna.
      std::size_t a1 = itsInfo->getAnt1()[i];
      std::size_t a2 = itsInfo->getAnt2()[i];
      if (!done[a1]) {
        frame.set(itsInfo->antennaPos()[a1]);
        testAzEl(converter, i, a1, itsInfo->getAnt1(), itsInfo->getAnt2());
        done[a1] = true;
      }
      // If needed, check if ant2 matches AzEl criterium.
      if (itsMatchBL[i] && !done[a2]) {
        frame.set(itsInfo->antennaPos()[a2]);
        testAzEl(converter, i, a2, itsInfo->getAnt1(), itsInfo->getAnt2());
        done[a2] = true;
      }
      if (itsMatchBL[i]) {
        match = true;
      }
    }
  }
  return match;
}

void PreFlagger::PSet::testAzEl(MDirection::Convert& converter,
                                unsigned int blnr, int ant,
                                const std::vector<int>& ant1,
                                const std::vector<int>& ant2) {
  // Calculate AzEl (n seconds because ranges are in seconds too).
  MVDirection mvAzel(converter().getValue());
  casacore::Vector<double> azel = mvAzel.getAngle("s").getValue();
  double az = azel[0];
  double el = azel[1];
  if (az < 0) az += 86400;
  if (el < 0) el += 86400;
  // If outside the ranges, there is no match.
  // Set no match for all baselines containing this antenna.
  // It needs to be done from this baseline on, because the earlier
  // baselines have already been handled.
  bool res = ((itsAzimuth.empty() || matchRange(az, itsAzimuth)) &&
              (itsElevation.empty() || matchRange(el, itsElevation)));
  if (!res) {
    for (unsigned int i = blnr; i < itsMatchBL.size(); ++i) {
      if (ant1[i] == ant || ant2[i] == ant) {
        itsMatchBL[i] = false;
      }
    }
  }
}

void PreFlagger::PSet::flagAmpl(const base::DPBuffer::DataType& data) {
  const std::size_t nr = data.shape(0) * data.shape(1);
  const std::size_t nrcorr = data.shape(2);

  const xt::xtensor<float, 3> amplitudes = xt::abs(data);
  const float* valPtr = amplitudes.data();
  uint8_t* flagPtr = itsFlags.data();
  for (unsigned int i = 0; i < nr; ++i) {
    bool flag = false;
    for (unsigned int j = 0; j < nrcorr; ++j) {
      if (valPtr[j] < itsAmplMin[j] || valPtr[j] > itsAmplMax[j]) {
        flag = true;
        break;
      }
    }
    if (!flag) {
      *flagPtr = false;
    }
    valPtr += nrcorr;
    ++flagPtr;
  }
}

void PreFlagger::PSet::flagPhase(const base::DPBuffer::DataType& data) {
  const std::size_t nr = data.shape(0) * data.shape(1);
  const std::size_t nrcorr = data.shape(2);

  const xt::xtensor<float, 3> phases = xt::arg(data);
  const float* valPtr = phases.data();
  uint8_t* flagPtr = itsFlags.data();
  for (unsigned int i = 0; i < nr; ++i) {
    bool flag = false;
    for (unsigned int j = 0; j < nrcorr; ++j) {
      if (valPtr[j] < itsPhaseMin[j] || valPtr[j] > itsPhaseMax[j]) {
        flag = true;
        break;
      }
    }
    if (!flag) {
      *flagPtr = false;
    }
    valPtr += nrcorr;
    ++flagPtr;
  }
}

void PreFlagger::PSet::flagReal(const base::DPBuffer::DataType& data) {
  const std::size_t nr = data.shape(0) * data.shape(1);
  const std::size_t nrcorr = data.shape(2);

  const std::complex<float>* valPtr = data.data();
  uint8_t* flagPtr = itsFlags.data();
  for (unsigned int i = 0; i < nr; ++i) {
    bool flag = false;
    for (unsigned int j = 0; j < nrcorr; ++j) {
      if (valPtr[j].real() < itsRealMin[j] ||
          valPtr[j].real() > itsRealMax[j]) {
        flag = true;
        break;
      }
    }
    if (!flag) {
      *flagPtr = false;
    }
    valPtr += nrcorr;
    ++flagPtr;
  }
}

void PreFlagger::PSet::flagImag(const base::DPBuffer::DataType& data) {
  const std::size_t nr = data.shape(0) * data.shape(1);
  const std::size_t nrcorr = data.shape(2);

  const std::complex<float>* valPtr = data.data();
  uint8_t* flagPtr = itsFlags.data();
  for (unsigned int i = 0; i < nr; ++i) {
    bool flag = false;
    for (unsigned int j = 0; j < nrcorr; ++j) {
      if (valPtr[j].imag() < itsImagMin[j] ||
          valPtr[j].imag() > itsImagMax[j]) {
        flag = true;
        break;
      }
    }
    if (!flag) {
      *flagPtr = false;
    }
    valPtr += nrcorr;
    ++flagPtr;
  }
}

void PreFlagger::PSet::flagChannels() {
  const size_t n_baselines = itsFlags.shape()[0];

  // 4x faster equivalent of:
  // itsFlags &= xt::broadcast(itsChanFlags, itsFlags.shape())
  for (size_t i = 0; i < n_baselines; ++i) {
    auto slice = xt::view(itsFlags, i, xt::all());

    std::transform(itsChanFlags.begin(), itsChanFlags.end(), slice.begin(),
                   slice.begin(), std::bit_and<uint8_t>{});
  }
}

// See http://montcs.bloomu.edu/~bobmon/Information/RPN/infix2rpn.shtml
// for the algorithm used here.
// Some code was added to check if no two subsequent operators or names
// are given.
std::vector<std::string> PreFlagger::PSet::exprToRpn(
    const std::string& origExpr) {
  // Operators & (or &&) | (or ||) and , are used as well as parentheses.
  // The operators must have a value in order of precedence, thus &&
  // has a higher precedence than || (as in C).
  string expr = boost::to_upper_copy(origExpr);
  std::stack<int> tokens;
  std::vector<std::string> names;
  unsigned int i = 0;
  bool hadName = false;  // the last token was a name.
  while (i < expr.size()) {
    int oper = 0;
    // skip whitespace.
    // Look for parenthesis or operator.
    if (expr[i] == ' ' || expr[i] == '\t') {
      i++;
    } else if (expr[i] == '(') {
      if (hadName)
        throw std::runtime_error(
            "no operator before opening parenthesis at pos. " +
            std::to_string(i) + " in expression: " + origExpr);
      oper = OpParen;
      tokens.push(oper);
      i++;
    } else if (expr[i] == '|') {
      oper = OpOr;
      i++;
      if (i < expr.size() && expr[i] == '|') i++;
    } else if (expr[i] == ',') {
      oper = OpOr;
      i++;
    } else if (expr[i] == '&') {
      oper = OpAnd;
      i++;
      if (i < expr.size() && expr[i] == '&') i++;
    } else if (expr[i] == '!') {
      oper = OpNot;
      i++;
    } else if (expr.size() - i >= 3 &&
               (expr.substr(i, 3) == "OR " || expr.substr(i, 3) == "OR\t" ||
                expr.substr(i, 3) == "OR!" || expr.substr(i, 3) == "OR(")) {
      oper = OpOr;
      i += 2;
    } else if (expr.size() - i == 2 && (expr.substr(i, 2) == "OR")) {
      oper = OpOr;
      i += 2;
    } else if (expr.size() - i >= 4 &&
               (expr.substr(i, 4) == "AND " || expr.substr(i, 4) == "AND\t" ||
                expr.substr(i, 4) == "AND!" || expr.substr(i, 4) == "AND(")) {
      oper = OpAnd;
      i += 3;
    } else if (expr.size() - i == 3 && (expr.substr(i, 3) == "AND")) {
      oper = OpAnd;
      i += 3;
    } else if (expr.size() - i >= 4 &&
               (expr.substr(i, 4) == "NOT " || expr.substr(i, 4) == "NOT\t" ||
                expr.substr(i, 4) == "NOT(")) {
      oper = OpNot;
      i += 3;
    } else if (expr.size() - i == 3 && (expr.substr(i, 3) == "NOT")) {
      oper = OpNot;
      i += 3;
    } else if (expr[i] == ')') {
      // Closing parenthesis. Push till opening parenthesis found.
      if (!hadName)
        throw std::runtime_error(
            "no set name before closing parenthesis at pos. " +
            std::to_string(i) + " in expression: " + origExpr);
      while (true) {
        if (tokens.empty())
          throw std::runtime_error("mismatched parentheses at pos. " +
                                   std::to_string(i) +
                                   " in expression: " + origExpr);
        if (tokens.top() == OpParen) {
          tokens.pop();
          break;
        }
        itsRpn.push_back(tokens.top());
        tokens.pop();
      }
      i++;
    } else {
      // No operator, thus it must be an operand (a set name).
      int st = i;
      if (hadName)
        throw std::runtime_error("No operator between set names at pos. " +
                                 std::to_string(i) +
                                 " in expression: " + origExpr);
      while (i < expr.size() && expr[i] != ' ' && expr[i] != '\t' &&
             expr[i] != '(' && expr[i] != ')' && expr[i] != '!' &&
             expr[i] != ',' && expr[i] != '&' && expr[i] != '|') {
        i++;
      }
      hadName = true;
      itsRpn.push_back(names.size());
      casacore::String setName(origExpr.substr(st, i - st));
      // Check the name is valid (no special characters).
      if (!setName.matches(casacore::RXidentifier))
        throw std::runtime_error("Invalid set name " + std::string(setName) +
                                 " used in set expression " + origExpr);
      names.push_back(setName);
    }
    if (oper < OpParen) {
      // Check if an operator was preceeded correctly.
      if (oper == OpNot) {
        if (hadName)
          throw std::runtime_error("No set name before operator ! at pos. " +
                                   std::to_string(i) +
                                   " in expression: " + origExpr);
      } else {
        if (!hadName)
          throw std::runtime_error("No set name before operator at pos. " +
                                   std::to_string(i) +
                                   " in expression: " + origExpr);
      }
      hadName = false;
      // Push till lower precedence found.
      while (!tokens.empty() && tokens.top() < oper) {
        itsRpn.push_back(tokens.top());
        tokens.pop();
      }
      tokens.push(oper);
    }
  }
  if (!hadName)
    throw std::runtime_error("no set name after last operator in expression: " +
                             origExpr);
  while (!tokens.empty()) {
    if (tokens.top() >= OpParen)
      throw std::runtime_error("mismatched parentheses in expression: " +
                               origExpr);
    itsRpn.push_back(tokens.top());
    tokens.pop();
  }
  return names;
}

std::vector<double> PreFlagger::PSet::fillTimes(
    const std::vector<std::string>& vec, bool asTime, bool canEndBeforeStart) {
  std::vector<double> result;
  result.reserve(2 * vec.size());
  // A time range can be given as time..time or time+-value.
  for (std::vector<std::string>::const_iterator str = vec.begin();
       str != vec.end(); ++str) {
    // Find the .. or +- token.
    bool usepm = false;
    std::string::size_type pos = str->find("..");
    if (pos == std::string::npos) {
      usepm = true;
      pos = str->find("+-");
      if (pos == std::string::npos)
        throw std::runtime_error("PreFlagger time range '" + *str +
                                 "' should be range using .. or +-");
    }
    // Get the time or datetime in seconds. The values must be positive.
    double v1 = getSeconds(str->substr(0, pos), asTime, false);
    double v2 = getSeconds(str->substr(pos + 2), asTime, usepm);
    if (v1 < 0 || v2 < 0)
      throw std::runtime_error("PreFlagger time range " + *str +
                               " must have positive values");
    if (usepm) {
      double pm = v2;
      v2 = v1 + pm;
      v1 -= pm;
    }
    // If time is used, values around midnight can be given.
    // Note there are 86400 seconds in a day.
    // They are split in 2 ranges.
    if (!canEndBeforeStart) {
      if (v1 >= v2)
        throw std::runtime_error("PreFlagger time range " + *str +
                                 " is invalid");
    } else {
      if (v1 < 0) {
        v1 += 86400;
      }
      if (v2 > 86400) {
        v2 -= 86400;
      }
    }
    if (v1 < v2) {
      result.push_back(v1);
      result.push_back(v2);
    } else {
      result.push_back(-1);
      result.push_back(v2);
      result.push_back(v1);
      result.push_back(86401);
    }
  }
  return result;
}

double PreFlagger::PSet::getSeconds(const std::string& str, bool asTime,
                                    bool usepm) {
  Quantity q;
  if (asTime || usepm) {
    if (!MVAngle::read(q, str, true))
      throw std::runtime_error("PreFlagger time " + str + " is invalid");
  } else {
    // It should be a proper date/time, so MVAngle::read should fail.
    if (MVAngle::read(q, str, true))
      throw std::runtime_error("PreFlagger datetime " + str +
                               " is not a proper date/time");
    if (!MVTime::read(q, str, true))
      throw std::runtime_error("PreFlagger datetime " + str + " is invalid" +
                               " is not a proper date/time");
  }
  double v = q.getValue("s");
  if (usepm) {
    if (v <= 0)
      throw std::runtime_error("Preflagger time plusminus value " + str +
                               " must be positive");
  }
  return v;
}

std::vector<float> PreFlagger::PSet::fillValuePerCorr(
    const common::ParameterValue& value, float defVal, bool& doFlag) {
  // Initialize with the default value per correlation.
  std::vector<float> result(4);
  std::fill(result.begin(), result.end(), defVal);
  if (!value.get().empty()) {
    if (value.isVector()) {
      // Defined as a vector, take the values given.
      std::vector<std::string> valstr = value.getStringVector();
      unsigned int sz = std::min(valstr.size(), result.size());
      if (sz > 0) {
        // It contains a value, so set that flagging is done.
        doFlag = true;
        for (unsigned int i = 0; i < sz; ++i) {
          if (!valstr[i].empty()) {
            result[i] = common::strToFloat(valstr[i]);
          }
        }
      }
    } else {
      // A single value means use it for all correlations.
      doFlag = true;
      std::fill(result.begin(), result.end(), value.getFloat());
    }
  }
  return result;
}

void PreFlagger::PSet::fillBLMatrix() {
  itsFlagOnBL = itsSelBL.hasSelection();
  if (itsFlagOnBL) {
    itsFlagBL.reference(itsSelBL.apply(*itsInfo));
  }
}

xt::xtensor<int, 1> PreFlagger::PSet::handleFreqRanges(
    const std::vector<double>& chanFreqs) {
  unsigned int nrchan = chanFreqs.size();
  xt::xtensor<int, 1> selChan({nrchan}, false);
  // A frequency range can be given as  value..value or value+-value.
  // Units can be given for each value; if one is given it applies to both.
  // Default unit is MHz.
  for (std::vector<std::string>::const_iterator str = itsStrFreq.begin();
       str != itsStrFreq.end(); ++str) {
    // Find the .. or +- token.
    bool usepm = false;
    std::string::size_type pos = str->find("..");
    if (pos == std::string::npos) {
      usepm = true;
      pos = str->find("+-");
      if (pos == std::string::npos)
        throw std::runtime_error("PreFlagger freqrange '" + *str +
                                 "' should be range using .. or +-");
    }
    string str1 = str->substr(0, pos);
    string str2 = str->substr(pos + 2);
    casacore::String u1, u2;
    double v1, v2;
    getValue(str1, v1, u1);
    // Default unit for 2nd value is that of 1st value.
    u2 = u1;
    getValue(str2, v2, u2);
    // If no unit, use MHz.
    if (u2.empty()) {
      u2 = "MHz";
    }
    // Default unit of 1st value is that of 2nd value.
    if (u1.empty()) {
      u1 = u2;
    }
    v1 = getFreqHz(v1, u1);
    v2 = getFreqHz(v2, u2);
    if (usepm) {
      double pm = v2;
      v2 = v1 + pm;
      v1 -= pm;
    }
    // Add any channel inside this range.
    for (unsigned int i = 0; i < chanFreqs.size(); ++i) {
      if (chanFreqs[i] > v1 && chanFreqs[i] < v2) {
        selChan[i] = true;
      }
    }
  }
  return selChan;
}

void PreFlagger::PSet::getValue(const std::string& str, double& value,
                                casacore::String& unit) {
  // See if a unit is given at the end.
  casacore::String v(str);
  // Remove possible trailing blanks.
  boost::algorithm::trim_right(v);
  casacore::Regex regex("[a-zA-Z]+$");
  std::string::size_type pos = v.index(regex);
  if (pos != casacore::String::npos) {
    unit = v.from(pos);
    v = v.before(pos);
  }
  // Set value and unit.
  value = common::strToDouble(v);
}

double PreFlagger::PSet::getFreqHz(double value, const casacore::String& unit) {
  Quantity q(value, unit);
  return q.getValue("Hz");
}

}  // namespace steps
}  // namespace dp3
