// DPInfo.h: General info about DPPP data processing attributes like averaging
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief General info about DPPP data processing attributes like averaging
/// @author Ger van Diepen

#ifndef DPPP_DPINFO_H
#define DPPP_DPINFO_H

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Record.h>

// Include appears unused
#include <aocommon/threadpool.h>

#include <boost/algorithm/string/case_conv.hpp>

#include <EveryBeam/correctionmode.h>

namespace dp3 {
namespace steps {
class InputStep;
}

namespace base {

/// @brief General info about DPPP data processing attributes like averaging

/// This class contains the information about the number of correlations,
/// channels, baselines, and times.
/// It is initialized by the first step and updated by steps like
/// Averager that change the number of channels or times.
/// Steps can take information from it to know about shapes.

class DPInfo {
 public:
  /// Default constructor.
  DPInfo();

  /// Set the initial info from the input.
  void init(unsigned int ncorr, unsigned int startChan, unsigned int nchan,
            unsigned int ntime, double startTime, double timeInterval,
            const string& msName, const string& antennaSet);

  void setIsBDAIntervalFactorInteger(bool isIntervalInteger) {
    itsBDAIntervalFactorInteger = isIntervalInteger;
  }

  /// Set nr of channels.
  void setNChan(unsigned int nchan) { itsNChan = nchan; }

  /// Set the time interval and the number of time steps.
  void setTimeIntervalAndSteps(double timeInterval, unsigned int ntime) {
    itsTimeInterval = timeInterval;
    itsNTime = ntime;
  }

  /// Set the frequency info.
  /// An empty resolutions or effectiveBW is default to chanWidths.
  /// itsTotalBW is set to the sum of effectiveBW.
  /// If refFreq is 0, it is set to the middle of chanFreqs (mean if even).
  void set(std::vector<double>&& chanFreqs, std::vector<double>&& chanWidths,
           std::vector<double>&& resolutions = std::vector<double>(),
           std::vector<double>&& effectiveBW = std::vector<double>(),
           double refFreq = 0);

  /// Set the frequency info, using different info per baseline.
  /// An empty resolutions or effectiveBW is default to chanWidths.
  /// itsTotalBW is set to the sum of effectiveBW, which should be equal for
  /// all baselines.
  /// If refFreq is 0, it is set to the middle of chanFreqs (mean if even).
  /// of the baseline with the most channels.
  void set(std::vector<std::vector<double>>&& chanFreqs,
           std::vector<std::vector<double>>&& chanWidths,
           std::vector<std::vector<double>>&& resolutions =
               std::vector<std::vector<double>>(),
           std::vector<std::vector<double>>&& effectiveBW =
               std::vector<std::vector<double>>(),
           double refFreq = 0);

  /// Set array info.
  void set(const casacore::MPosition& arrayPos,
           const casacore::MDirection& phaseCenter,
           const casacore::MDirection& delayCenter,
           const casacore::MDirection& tileBeamDir);

  /// Set the info for the given antennae and baselines.
  void set(const casacore::Vector<casacore::String>& antNames,
           const casacore::Vector<casacore::Double>& antDiam,
           const std::vector<casacore::MPosition>& antPos,
           const casacore::Vector<casacore::Int>& ant1,
           const casacore::Vector<casacore::Int>& ant2);

  /// Set the name of the data column
  void setDataColName(const std::string& dataColName) {
    itsDataColName = dataColName;
  }

  /// Set the name of the weight column
  void setWeightColName(const std::string& weightColName) {
    itsWeightColName = weightColName;
  }

  /// Update the info for the given average factors.
  /// If chanAvg is higher than the actual nr of channels, it is reset.
  /// The same is true for timeAvg.
  /// It returns the possibly reset nr of channels to average.
  unsigned int update(unsigned int chanAvg, unsigned int timeAvg);

  /// Update the info from the given selection parameters.
  /// Optionally unused stations are really removed from the antenna lists.
  void update(unsigned int startChan, unsigned int nchan,
              const std::vector<unsigned int>& baselines, bool remove);

  /// Update the info for the given average factors.
  void update(std::vector<unsigned int>&& timeAvg);

  /// Remove unused stations from the antenna lists.
  void removeUnusedAnt();

  /// Set the phase center.
  /// If original=true, it is set to the original phase center.
  void setPhaseCenter(const casacore::MDirection& phaseCenter, bool original) {
    itsPhaseCenter = phaseCenter;
    itsPhaseCenterIsOriginal = original;
  }

  /// Get the info.
  ///@{
  const string& msName() const { return itsMSName; }
  const string& antennaSet() const { return itsAntennaSet; }
  unsigned int ncorr() const { return itsNCorr; }
  unsigned int nchan() const { return itsNChan; }
  unsigned int startchan() const { return itsStartChan; }
  unsigned int origNChan() const { return itsOrigNChan; }
  unsigned int nchanAvg() const { return itsChanAvg; }
  unsigned int nantenna() const { return itsAntNames.size(); }
  unsigned int nbaselines() const { return itsAnt1.size(); }
  unsigned int ntime() const { return itsNTime; }
  unsigned int ntimeAvg(unsigned int baseline = 0) const {
    return itsTimeAvg[baseline];
  }
  const std::vector<unsigned int>& ntimeAvgs() const { return itsTimeAvg; }
  double startTime() const { return itsStartTime; }
  double timeInterval() const { return itsTimeInterval; }
  bool isBDAIntervalFactorInteger() const {
    return itsBDAIntervalFactorInteger;
  }
  const std::vector<std::size_t>& getAnt1() const { return itsAnt1; }
  const std::vector<std::size_t>& getAnt2() const { return itsAnt2; }
  const casacore::Vector<casacore::String>& antennaNames() const {
    return itsAntNames;
  }
  const casacore::Vector<casacore::Double>& antennaDiam() const {
    return itsAntDiam;
  }
  const std::vector<casacore::MPosition>& antennaPos() const {
    return itsAntPos;
  }
  const casacore::MPosition& arrayPos() const { return itsArrayPos; }
  const casacore::MPosition arrayPosCopy() const {
    return copyMeasure(casacore::MeasureHolder(itsArrayPos)).asMPosition();
  }
  const casacore::MDirection& phaseCenter() const { return itsPhaseCenter; }
  const casacore::MDirection phaseCenterCopy() const {
    return copyMeasure(casacore::MeasureHolder(itsPhaseCenter)).asMDirection();
  }
  bool phaseCenterIsOriginal() const { return itsPhaseCenterIsOriginal; }
  const casacore::MDirection& delayCenter() const { return itsDelayCenter; }
  const casacore::MDirection delayCenterCopy() const {
    return copyMeasure(casacore::MeasureHolder(itsDelayCenter)).asMDirection();
  }
  const casacore::MDirection& tileBeamDir() const { return itsTileBeamDir; }
  const casacore::MDirection tileBeamDirCopy() const {
    return copyMeasure(casacore::MeasureHolder(itsTileBeamDir)).asMDirection();
  }
  const std::vector<double>& chanFreqs(std::size_t baseline = 0) const {
    return itsChanFreqs[baseline];
  }
  const std::vector<double>& chanWidths(std::size_t baseline = 0) const {
    return itsChanWidths[baseline];
  }
  const std::vector<double>& resolutions(std::size_t baseline = 0) const {
    return itsResolutions[baseline];
  }
  const std::vector<double>& effectiveBW(std::size_t baseline = 0) const {
    return itsEffectiveBW[baseline];
  }
  bool hasBDAChannels() const { return itsChanFreqs.size() == nbaselines(); }
  const std::string& getDataColName() const { return itsDataColName; }
  const std::string& getWeightColName() const { return itsWeightColName; }
  double totalBW() const { return itsTotalBW; }
  double refFreq() const { return itsRefFreq; }
  ///@}

  /// Get the antenna numbers actually used in the (selected) baselines.
  /// E.g. [0,2,5,6]
  const std::vector<int>& antennaUsed() const { return itsAntUsed; }

  /// Get the indices of all antennae in the used antenna vector above.
  /// -1 means that the antenna is not used.
  /// E.g. [0,-1,1,-1,-1,2,3] for the example above.
  const std::vector<int>& antennaMap() const { return itsAntMap; }

  /// Are the visibility data needed?
  bool needVisData() const { return itsNeedVisData; }
  /// Does the last step need to write data and/or flags?
  bool needWrite() const {
    return itsWriteData || itsWriteFlags || itsWriteWeights;
  }
  bool writeData() const { return itsWriteData; }
  bool writeFlags() const { return itsWriteFlags; }
  bool writeWeights() const { return itsWriteWeights; }
  /// Has the meta data been changed in a step (precluding an update)?
  bool metaChanged() const { return itsMetaChanged; }

  /// Set if visibility data needs to be read.
  void setNeedVisData() { itsNeedVisData = true; }
  /// Set if data needs to be written.
  void setWriteData() { itsWriteData = true; }
  void setWriteFlags() { itsWriteFlags = true; }
  void setWriteWeights() { itsWriteWeights = true; }
  /// Clear all write flags.
  void clearWrites() { itsWriteData = itsWriteFlags = itsWriteWeights = false; }
  /// Set change of meta data.
  void setMetaChanged() { itsMetaChanged = true; }
  void clearMetaChanged() { itsMetaChanged = false; }

  /// Get the baseline table index of the autocorrelations.
  /// A negative value means there are no autocorrelations for that antenna.
  const std::vector<int>& getAutoCorrIndex() const;

  /// Get the lengths of the baselines (in meters).
  const std::vector<double>& getBaselineLengths() const;

  void setNThreads(unsigned int nThreads) { itsNThreads = nThreads; }

  unsigned int nThreads() const { return itsNThreads; }

  /// Convert to a Record.
  /// The names of the fields in the record are the data names without 'its'.
  casacore::Record toRecord() const;

  /// Update the DPInfo object from a Record.
  /// It is possible that only a few fields are defined in the record.
  void fromRecord(const casacore::Record& rec);

  everybeam::CorrectionMode beamCorrectionMode() const {
    return itsBeamCorrectionMode;
  };
  void setBeamCorrectionMode(everybeam::CorrectionMode mode) {
    itsBeamCorrectionMode = mode;
  }

  const casacore::MDirection& beamCorrectionDir() const {
    return itsBeamCorrectionDir;
  }
  void setBeamCorrectionDir(const casacore::MDirection& dir) {
    itsBeamCorrectionDir = dir;
  }

  /// Determine if the channels have a regular layout.
  bool channelsAreRegular() const;

 private:
  /// Set which antennae are actually used.
  void setAntUsed();

  /// Creates a real copy of a casacore::Measure by exporting to a Record
  static casacore::MeasureHolder copyMeasure(
      const casacore::MeasureHolder fromMeas);

  bool itsNeedVisData;   ///< Are the visibility data needed?
  bool itsWriteData;     ///< Must the data be written?
  bool itsWriteFlags;    ///< Must the flags be written?
  bool itsWriteWeights;  ///< Must the weights be written?
  bool itsMetaChanged;   ///< Are meta data changed? (e.g., by averaging)
  string itsMSName;
  std::string itsDataColName;
  std::string itsWeightColName;
  string itsAntennaSet;
  unsigned int itsNCorr;
  unsigned int itsStartChan;
  unsigned int itsOrigNChan;
  unsigned int itsNChan;
  unsigned int itsChanAvg;
  unsigned int itsNTime;
  std::vector<unsigned int> itsTimeAvg;
  double itsStartTime;
  double itsTimeInterval;
  bool itsBDAIntervalFactorInteger;  // INTEGER_INTERVAL_FACTORS in
                                     // [BDA_TIME_AXIS}
  casacore::MDirection itsPhaseCenter;
  bool itsPhaseCenterIsOriginal;
  casacore::MDirection itsDelayCenter;
  casacore::MDirection itsTileBeamDir;
  everybeam::CorrectionMode itsBeamCorrectionMode;
  casacore::MDirection itsBeamCorrectionDir;
  casacore::MPosition itsArrayPos;
  /// Channels may differ between baselines when using BDA.
  /// If all baselines have equal channels, the outer vector only holds one
  /// inner vector. When using BDA, each baseline has its own inner vector.
  std::vector<std::vector<double>> itsChanFreqs;
  std::vector<std::vector<double>> itsChanWidths;
  std::vector<std::vector<double>> itsResolutions;
  std::vector<std::vector<double>> itsEffectiveBW;
  double itsTotalBW;
  double itsRefFreq;
  casacore::Vector<casacore::String> itsAntNames;
  casacore::Vector<casacore::Double> itsAntDiam;
  std::vector<casacore::MPosition> itsAntPos;
  std::vector<int> itsAntUsed;
  std::vector<int> itsAntMap;
  std::vector<std::size_t> itsAnt1;           ///< ant1 of all baselines
  std::vector<std::size_t> itsAnt2;           ///< ant2 of all baselines
  mutable std::vector<double> itsBLength;     ///< baseline lengths
  mutable std::vector<int> itsAutoCorrIndex;  ///< autocorr index per ant
  unsigned int itsNThreads;
};

}  // namespace base
}  // namespace dp3

#endif
