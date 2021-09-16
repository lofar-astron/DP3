// DPInfo.cc: General info about DPPP data processing attributes like averaging
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DPInfo.h"
#include "Exceptions.h"

#include "../common/Epsilon.h"

#include "../steps/InputStep.h"

#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/STLIO.h>

#include <aocommon/threadpool.h>

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace casacore;
using namespace std;

namespace dp3 {
namespace base {

DPInfo::DPInfo()
    : itsNeedVisData(false),
      itsWriteData(false),
      itsWriteFlags(false),
      itsWriteWeights(false),
      itsMetaChanged(false),
      itsNCorr(0),
      itsStartChan(0),
      itsNChan(0),
      itsChanAvg(1),
      itsNTime(0),
      itsTimeAvg({1}),
      itsStartTime(0),
      itsTimeInterval(0),
      itsPhaseCenterIsOriginal(true),
      itsBeamCorrectionMode(everybeam::CorrectionMode::kNone),
      itsNThreads(aocommon::ThreadPool::NCPUs()) {}

void DPInfo::init(unsigned int ncorr, unsigned int startChan,
                  unsigned int nchan, unsigned int ntime, double startTime,
                  double timeInterval, const string& msName,
                  const string& antennaSet) {
  itsNCorr = ncorr;
  itsStartChan = startChan;
  itsNChan = nchan;
  itsOrigNChan = nchan;
  itsNTime = ntime;
  itsStartTime = startTime;
  itsTimeInterval = timeInterval;
  itsMSName = msName;
  itsAntennaSet = antennaSet;
}

void DPInfo::set(std::vector<double>&& chan_freqs,
                 std::vector<double>&& chan_widths,
                 std::vector<double>&& resolutions,
                 std::vector<double>&& effective_bw, double ref_freq) {
  if (resolutions.empty()) {
    resolutions = chan_widths;
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  }

  if (ref_freq == 0) {
    int n = chan_freqs.size();
    // Takes mean of middle elements if n is even; takes middle if odd.
    ref_freq = 0.5 * (chan_freqs[(n - 1) / 2] + chan_freqs[n / 2]);
  }

  itsChanFreqs.clear();
  itsChanWidths.clear();
  itsResolutions.clear();
  itsEffectiveBW.clear();

  itsChanFreqs.push_back(std::move(chan_freqs));
  itsChanWidths.push_back(std::move(chan_widths));
  itsResolutions.push_back(std::move(resolutions));
  itsEffectiveBW.push_back(std::move(effective_bw));

  itsTotalBW = std::accumulate(itsEffectiveBW.front().begin(),
                               itsEffectiveBW.front().end(), 0.0);
  itsRefFreq = ref_freq;
}

void DPInfo::set(std::vector<std::vector<double>>&& chan_freqs,
                 std::vector<std::vector<double>>&& chan_widths,
                 std::vector<std::vector<double>>&& resolutions,
                 std::vector<std::vector<double>>&& effective_bw,
                 double ref_freq) {
  if (resolutions.empty()) {
    resolutions = chan_widths;
  }
  if (effective_bw.empty()) {
    effective_bw = chan_widths;
  }
  if (chan_freqs.size() != nbaselines() || chan_widths.size() != nbaselines() ||
      resolutions.size() != nbaselines() ||
      effective_bw.size() != nbaselines()) {
    throw Exception("Invalid baseline count while setting frequency info");
  }

  const double total_bw = std::accumulate(effective_bw.front().begin(),
                                          effective_bw.front().end(), 0.0);
  for (std::vector<double>& eff_bw_bl : effective_bw) {
    if (std::accumulate(eff_bw_bl.begin(), eff_bw_bl.end(), 0.0) != total_bw) {
      throw Exception("Total BW is not equal for all baselines");
    }
  }

  // Find the baseline with the most channels.
  auto comp = [](const std::vector<double>& left,
                 const std::vector<double>& right) {
    return left.size() < right.size();
  };
  auto it = std::max_element(chan_freqs.begin(), chan_freqs.end(), comp);
  itsNChan = it->size();

  if (ref_freq == 0) {
    // Takes mean of middle elements if n is even; takes middle if odd.
    ref_freq = 0.5 * ((*it)[(itsNChan - 1) / 2] + (*it)[itsNChan / 2]);
  }

  itsChanFreqs = std::move(chan_freqs);
  itsChanWidths = std::move(chan_widths);
  itsResolutions = std::move(resolutions);
  itsEffectiveBW = std::move(effective_bw);
  itsTotalBW = total_bw;
  itsRefFreq = ref_freq;
}

bool DPInfo::channelsAreRegular() const {
  if (itsChanFreqs.empty()) {
    return true;
  }

  // Check that all baselines have equal channel layouts.
  const double kTolerance = 1.0;  // Hz
  for (std::size_t bl = 1; bl < itsChanFreqs.size(); ++bl) {
    if (!common::EpsilonEqual(itsChanFreqs.front(), itsChanFreqs[bl],
                              kTolerance) ||
        !common::EpsilonEqual(itsChanWidths.front(), itsChanWidths[bl],
                              kTolerance) ||
        !common::EpsilonEqual(itsResolutions.front(), itsResolutions[bl],
                              kTolerance) ||
        !common::EpsilonEqual(itsEffectiveBW.front(), itsEffectiveBW[bl],
                              kTolerance)) {
      return false;
    }
  }

  // Check that channels are evenly spaced.
  const std::vector<double>& freqs = itsChanFreqs.front();
  const std::vector<double>& widths = itsChanWidths.front();
  if (freqs.size() > 1) {
    const double freqstep0 = freqs[1] - freqs[0];
    const double kTolerance = 1.e3;  // Compare up to 1kHz accuracy.
    for (std::size_t i = 1; i < freqs.size(); ++i) {
      if ((std::abs(freqs[i] - freqs[i - 1] - freqstep0) >= kTolerance) ||
          (std::abs(widths[i] - widths[0]) >= kTolerance)) {
        return false;
      }
    }
  }

  return true;
}

void DPInfo::set(const MPosition& arrayPos, const MDirection& phaseCenter,
                 const MDirection& delayCenter, const MDirection& tileBeamDir) {
  itsArrayPos = arrayPos;
  itsPhaseCenter = phaseCenter;
  itsDelayCenter = delayCenter;
  itsTileBeamDir = tileBeamDir;
}

void DPInfo::set(const Vector<casacore::String>& antNames,
                 const Vector<Double>& antDiam, const vector<MPosition>& antPos,
                 const Vector<Int>& ant1, const Vector<Int>& ant2) {
  if (antNames.size() != antDiam.size() || antNames.size() != antPos.size())
    throw std::invalid_argument(
        "The name, diameter and position arrays are not of the same size");
  if (ant1.size() != ant2.size())
    throw std::invalid_argument(
        "The ant1 and ant2 arrays are not of the same size");
  itsAntNames.reference(antNames);
  itsAntDiam.reference(antDiam);
  itsAntPos = antPos;
  itsAnt1 = std::vector<std::size_t>(ant1.begin(), ant1.end());
  itsAnt2 = std::vector<std::size_t>(ant2.begin(), ant2.end());
  // Set which antennae are used.
  setAntUsed();
}

void DPInfo::setAntUsed() {
  itsAntUsed.clear();
  itsAntMap.resize(itsAntNames.size());
  std::fill(itsAntMap.begin(), itsAntMap.end(), -1);
  for (unsigned int i = 0; i < itsAnt1.size(); ++i) {
    if (itsAnt1[i] >= itsAntMap.size() || itsAnt2[i] >= itsAntMap.size())
      throw std::runtime_error("Antenna map has an inconsistent size");
    itsAntMap[itsAnt1[i]] = 0;
    itsAntMap[itsAnt2[i]] = 0;
  }
  itsAntUsed.reserve(itsAntNames.size());
  for (unsigned int i = 0; i < itsAntMap.size(); ++i) {
    if (itsAntMap[i] == 0) {
      itsAntMap[i] = itsAntUsed.size();
      itsAntUsed.push_back(i);
    }
  }
}

MeasureHolder DPInfo::copyMeasure(const MeasureHolder fromMeas) {
  Record rec;
  String msg;
  if (!fromMeas.toRecord(msg, rec))
    throw std::runtime_error("Could not copy MeasureHolder record to record");
  MeasureHolder mh2;
  if (!mh2.fromRecord(msg, rec))
    throw std::runtime_error("Could not copy record to MeasureHolder");
  return mh2;
}

unsigned int DPInfo::update(unsigned int chanAvg, unsigned int timeAvg) {
  if (itsChanFreqs.size() != 1) {
    throw Exception("Averaging does not support BDA");
  }

  if (chanAvg > itsNChan) {
    chanAvg = itsNChan;
  }
  if (timeAvg > itsNTime) {
    timeAvg = itsNTime;
  }
  if (itsNChan % chanAvg != 0)
    throw Exception(
        "When averaging, nr of channels must divide integrally; "
        "itsNChan=" +
        std::to_string(itsNChan) + " chanAvg=" + std::to_string(chanAvg));
  itsChanAvg *= chanAvg;
  itsNChan = (itsNChan + chanAvg - 1) / chanAvg;
  itsTimeAvg.front() *= timeAvg;
  itsNTime = (itsNTime + timeAvg - 1) / timeAvg;
  itsTimeInterval *= timeAvg;
  std::vector<double> freqs(itsNChan);
  std::vector<double> widths(itsNChan, 0.);
  std::vector<double> resols(itsNChan, 0.);
  std::vector<double> effBWs(itsNChan, 0.);
  double totBW = 0;
  for (unsigned int i = 0; i < itsNChan; ++i) {
    freqs[i] = 0.5 * (itsChanFreqs.front()[i * chanAvg] +
                      itsChanFreqs.front()[(i + 1) * chanAvg - 1]);
    for (unsigned int j = 0; j < chanAvg; ++j) {
      widths[i] += itsChanWidths.front()[i * chanAvg + j];
      resols[i] += itsResolutions.front()[i * chanAvg + j];
      effBWs[i] += itsEffectiveBW.front()[i * chanAvg + j];
    }
    totBW += effBWs[i];
  }
  itsChanFreqs.front() = std::move(freqs);
  itsChanWidths.front() = std::move(widths);
  itsResolutions.front() = std::move(resols);
  itsEffectiveBW.front() = std::move(effBWs);
  itsTotalBW = totBW;
  return chanAvg;
}

void DPInfo::update(std::vector<unsigned int>&& timeAvg) {
  itsTimeAvg = std::move(timeAvg);
}

void DPInfo::update(unsigned int startChan, unsigned int nchan,
                    const vector<unsigned int>& baselines, bool removeAnt) {
  if (itsChanFreqs.size() != 1) {
    throw Exception("Channel selection does not support BDA");
  }
  itsStartChan = startChan;
  auto freqs_begin = itsChanFreqs.front().begin() + startChan;
  auto widths_begin = itsChanWidths.front().begin() + startChan;
  auto resol_begin = itsResolutions.front().begin() + startChan;
  auto effbw_begin = itsEffectiveBW.front().begin() + startChan;
  itsChanFreqs.front() = std::vector<double>(freqs_begin, freqs_begin + nchan);
  itsChanWidths.front() =
      std::vector<double>(widths_begin, widths_begin + nchan);
  itsResolutions.front() =
      std::vector<double>(resol_begin, resol_begin + nchan);
  itsEffectiveBW.front() =
      std::vector<double>(effbw_begin, effbw_begin + nchan);
  itsNChan = nchan;
  // Keep only selected baselines.
  if (!baselines.empty()) {
    std::vector<std::size_t> ant1(baselines.size());
    std::vector<std::size_t> ant2(baselines.size());
    for (unsigned int i = 0; i < baselines.size(); ++i) {
      ant1[i] = itsAnt1[baselines[i]];
      ant2[i] = itsAnt2[baselines[i]];
    }
    itsAnt1 = std::move(ant1);
    itsAnt2 = std::move(ant2);
    // Clear; they'll be recalculated if needed.
    itsBLength.resize(0);
    itsAutoCorrIndex.resize(0);
  }
  setAntUsed();
  // If needed, remove the stations and renumber the baselines.
  if (removeAnt) {
    removeUnusedAnt();
  }
}

void DPInfo::removeUnusedAnt() {
  if (itsAntUsed.size() < itsAntMap.size()) {
    // First remove stations.
    Vector<casacore::String> antNames(itsAntUsed.size());
    Vector<Double> antDiam(itsAntUsed.size());
    vector<MPosition> antPos;
    antPos.reserve(itsAntUsed.size());
    for (unsigned int i = 0; i < itsAntUsed.size(); ++i) {
      antNames[i] = itsAntNames[itsAntUsed[i]];
      antDiam[i] = itsAntDiam[itsAntUsed[i]];
      antPos.push_back(itsAntPos[itsAntUsed[i]]);
    }
    // Use the new vectors.
    itsAntNames.reference(antNames);
    itsAntDiam.reference(antDiam);
    itsAntPos.swap(antPos);
    // Renumber the baselines.
    for (unsigned int i = 0; i < itsAnt1.size(); ++i) {
      itsAnt1[i] = itsAntMap[itsAnt1[i]];
      itsAnt2[i] = itsAntMap[itsAnt2[i]];
    }
    // Now fill the itsAntUsed and itsAntMap vectors again.
    setAntUsed();
    // Clear; they'll be recalculated if needed.
    itsBLength.resize(0);
    itsAutoCorrIndex.resize(0);
  }
}

const vector<double>& DPInfo::getBaselineLengths() const {
  // Calculate the baseline lengths if not done yet.
  if (itsBLength.empty()) {
    // First get the antenna positions.
    const vector<MPosition>& antPos = antennaPos();
    vector<Vector<double>> antVec;
    antVec.reserve(antPos.size());
    for (vector<MPosition>::const_iterator iter = antPos.begin();
         iter != antPos.end(); ++iter) {
      // Convert to ITRF and keep as x,y,z in m.
      antVec.push_back(
          MPosition::Convert(*iter, MPosition::ITRF)().getValue().getValue());
    }
    // Fill in the length of each baseline.
    vector<double> blength;
    itsBLength.reserve(itsAnt1.size());
    for (unsigned int i = 0; i < itsAnt1.size(); ++i) {
      Array<double> diff(antVec[itsAnt2[i]] - antVec[itsAnt1[i]]);
      itsBLength.push_back(sqrt(sum(diff * diff)));
    }
  }
  return itsBLength;
}

const vector<int>& DPInfo::getAutoCorrIndex() const {
  if (itsAutoCorrIndex.empty()) {
    int nant = 1 + std::max(*std::max_element(itsAnt1.begin(), itsAnt1.end()),
                            *std::max_element(itsAnt2.begin(), itsAnt2.end()));
    itsAutoCorrIndex.resize(nant);
    std::fill(itsAutoCorrIndex.begin(), itsAutoCorrIndex.end(), -1);
    // Keep the baseline table index for the autocorrelations.
    for (unsigned int i = 0; i < itsAnt1.size(); ++i) {
      if (itsAnt1[i] == itsAnt2[i]) {
        itsAutoCorrIndex[itsAnt1[i]] = i;
      }
    }
  }
  return itsAutoCorrIndex;
}

Record DPInfo::toRecord() const {
  Record rec;
  rec.define("NeedVisData", itsNeedVisData);
  rec.define("WriteData", itsWriteData);
  rec.define("WriteFlags", itsWriteFlags);
  rec.define("WriteWeights", itsWriteWeights);
  rec.define("MetaChanged", itsMetaChanged);
  rec.define("MSName", itsMSName);
  rec.define("AntennaSet", itsAntennaSet);
  rec.define("NCorr", itsNCorr);
  rec.define("StartChan", itsStartChan);
  rec.define("OrigNChan", itsOrigNChan);
  rec.define("NChan", itsNChan);
  rec.define("ChanAvg", itsChanAvg);
  rec.define("NTime", itsNTime);
  rec.define("TimeAvg", itsTimeAvg.front());
  rec.define("StartTime", itsStartTime);
  rec.define("TimeInterval", itsTimeInterval);
  rec.define("ChanFreqs", casacore::Vector<double>(itsChanFreqs.front()));
  rec.define("ChanWidths", casacore::Vector<double>(itsChanWidths.front()));
  rec.define("Resolutions", casacore::Vector<double>(itsResolutions.front()));
  rec.define("EffectiveBW", casacore::Vector<double>(itsEffectiveBW.front()));
  rec.define("TotalBW", itsTotalBW);
  rec.define("RefFreq", itsRefFreq);
  rec.define("AntNames", itsAntNames);
  rec.define("AntDiam", itsAntDiam);
  rec.define("AntUsed", Vector<int>(itsAntUsed));
  rec.define("AntMap", Vector<int>(itsAntMap));
  rec.define("Ant1", Vector<int>(itsAnt1));
  rec.define("Ant2", Vector<int>(itsAnt2));
  rec.define("BLength", Vector<double>(itsBLength));
  rec.define("AutoCorrIndex", Vector<int>(itsAutoCorrIndex));
  return rec;
}

void DPInfo::fromRecord(const Record& rec) {
  if (rec.isDefined("NeedVisData")) {
    rec.get("NeedVisData", itsNeedVisData);
  }
  if (rec.isDefined("WriteData")) {
    rec.get("WriteData", itsWriteData);
  }
  if (rec.isDefined("WriteFlags")) {
    rec.get("WriteFlags", itsWriteFlags);
  }
  if (rec.isDefined("WriteWeights")) {
    rec.get("WriteWeights", itsWriteWeights);
  }
  if (rec.isDefined("MetaChanged")) {
    rec.get("MetaChanged", itsMetaChanged);
  }
  if (rec.isDefined("MSName")) {
    itsMSName = rec.asString("MSName");
  }
  if (rec.isDefined("AntennaSet")) {
    itsAntennaSet = rec.asString("AntennaSet");
  }
  if (rec.isDefined("NCorr")) {
    rec.get("NCorr", itsNCorr);
  }
  if (rec.isDefined("StartChan")) {
    rec.get("StartChan", itsStartChan);
  }
  if (rec.isDefined("OrigNChan")) {
    rec.get("OrigNChan", itsOrigNChan);
  }
  if (rec.isDefined("NChan")) {
    rec.get("NChan", itsNChan);
  }
  if (rec.isDefined("ChanAvg")) {
    rec.get("ChanAvg", itsChanAvg);
  }
  if (rec.isDefined("NTime")) {
    rec.get("NTime", itsNTime);
  }
  if (rec.isDefined("TimeAvg")) {
    rec.get("TimeAvg", itsTimeAvg.front());
  }
  if (rec.isDefined("StartTime")) {
    rec.get("StartTime", itsStartTime);
  }
  if (rec.isDefined("TimeInterval")) {
    rec.get("TimeInterval", itsTimeInterval);
  }
  if (rec.isDefined("ChanFreqs")) {
    itsChanFreqs.clear();
    itsChanFreqs.push_back(rec.toArrayDouble("ChanFreqs").tovector());
  }
  if (rec.isDefined("ChanWidths")) {
    itsChanWidths.clear();
    itsChanWidths.push_back(rec.toArrayDouble("ChanWidths").tovector());
  }
  if (rec.isDefined("Resolutions")) {
    itsResolutions.clear();
    itsResolutions.push_back(rec.toArrayDouble("Resolutions").tovector());
  }
  if (rec.isDefined("EffectiveBW")) {
    itsEffectiveBW.clear();
    itsEffectiveBW.push_back(rec.toArrayDouble("EffectiveBW").tovector());
  }
  if (rec.isDefined("TotalBW")) {
    rec.get("TotalBW", itsTotalBW);
  }
  if (rec.isDefined("RefFreq")) {
    rec.get("RefFreq", itsRefFreq);
  }
  if (rec.isDefined("AntNames")) {
    rec.get("AntNames", itsAntNames);
  }
  if (rec.isDefined("AntDiam")) {
    rec.get("AntDiam", itsAntDiam);
  }
  /// if (rec.isDefined ("AntUsed")) {
  /// itsAntUsed = rec.toArrayInt("AntUsed").tovector();
  ///}
  /// if (rec.isDefined ("AntMap")) {
  ///  itsAntMap = rec.toArrayInt("AntMap").tovector();
  ///}
  if (rec.isDefined("Ant1")) {
    casacore::Vector<casacore::Int> ant1 = rec.toArrayInt("Ant1");
    itsAnt1 = std::vector<std::size_t>(ant1.begin(), ant1.end());
  }
  if (rec.isDefined("Ant2")) {
    casacore::Vector<casacore::Int> ant2 = rec.toArrayInt("Ant2");
    itsAnt2 = std::vector<std::size_t>(ant2.begin(), ant2.end());
  }
  /// if (rec.isDefined ("BLength")) {
  ///  itsBLength = rec.toArrayDouble("BLength").tovector();
  ///}
  /// if (rec.isDefined ("AutoCorrIndex")) {
  ///  itsAutoCorrIndex = rec.toArrayInt("AutoCorrIndex").tovector();
  ///}
}

}  // namespace base
}  // namespace dp3
