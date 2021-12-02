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
    bda_interval_factor_is_integer_ = isIntervalInteger;
  }

  /// Set nr of channels.
  void setNChan(unsigned int nchan) { n_channels_ = nchan; }

  /// Set the time interval and the number of time steps.
  void setTimeIntervalAndSteps(double timeInterval, unsigned int ntime) {
    time_interval_ = timeInterval;
    n_times_ = ntime;
  }

  /// Set the frequency info.
  /// An empty resolutions or effectiveBW is default to chanWidths.
  /// total_bandwidth_ is set to the sum of effectiveBW.
  /// If refFreq is 0, it is set to the middle of chanFreqs (mean if even).
  void set(std::vector<double>&& chanFreqs, std::vector<double>&& chanWidths,
           std::vector<double>&& resolutions = std::vector<double>(),
           std::vector<double>&& effectiveBW = std::vector<double>(),
           double refFreq = 0);

  /// Set the frequency info, using different info per baseline.
  /// An empty resolutions or effectiveBW is default to chanWidths.
  /// total_bandwidth_ is set to the sum of effectiveBW, which should be equal
  /// for all baselines. If refFreq is 0, it is set to the middle of chanFreqs
  /// (mean if even). of the baseline with the most channels.
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
    data_column_name_ = dataColName;
  }

  /// Set the name of the flag column
  void setFlagColName(const std::string& flagColName) {
    flag_column_name_ = flagColName;
  }

  /// Set the name of the weight column
  void setWeightColName(const std::string& weightColName) {
    weight_column_name_ = weightColName;
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
    phase_center_ = phaseCenter;
    phase_center_is_original_ = original;
  }

  /// Get the info.
  ///@{
  const string& msName() const { return ms_name_; }
  const string& antennaSet() const { return antenna_set_; }
  unsigned int ncorr() const { return n_correlations_; }
  unsigned int nchan() const { return n_channels_; }
  unsigned int startchan() const { return start_channel_; }
  unsigned int origNChan() const { return original_n_channels_; }
  unsigned int nchanAvg() const { return channel_averaging_factor_; }
  unsigned int nantenna() const { return antenna_names_.size(); }
  unsigned int nbaselines() const { return antenna1_.size(); }
  unsigned int ntime() const { return n_times_; }
  unsigned int ntimeAvg(unsigned int baseline = 0) const {
    return time_averaging_factors_[baseline];
  }
  const std::vector<unsigned int>& ntimeAvgs() const {
    return time_averaging_factors_;
  }
  double startTime() const { return start_time_; }
  double timeInterval() const { return time_interval_; }
  bool isBDAIntervalFactorInteger() const {
    return bda_interval_factor_is_integer_;
  }
  const std::vector<std::size_t>& getAnt1() const { return antenna1_; }
  const std::vector<std::size_t>& getAnt2() const { return antenna2_; }
  const casacore::Vector<casacore::String>& antennaNames() const {
    return antenna_names_;
  }
  const casacore::Vector<casacore::Double>& antennaDiam() const {
    return antenna_diameters_;
  }
  const std::vector<casacore::MPosition>& antennaPos() const {
    return antenna_positions_;
  }
  const casacore::MPosition& arrayPos() const { return array_position_; }
  const casacore::MPosition arrayPosCopy() const {
    return copyMeasure(casacore::MeasureHolder(array_position_)).asMPosition();
  }
  const casacore::MDirection& phaseCenter() const { return phase_center_; }
  const casacore::MDirection phaseCenterCopy() const {
    return copyMeasure(casacore::MeasureHolder(phase_center_)).asMDirection();
  }
  bool phaseCenterIsOriginal() const { return phase_center_is_original_; }
  const casacore::MDirection& delayCenter() const { return delay_center_; }
  const casacore::MDirection delayCenterCopy() const {
    return copyMeasure(casacore::MeasureHolder(delay_center_)).asMDirection();
  }
  const casacore::MDirection& tileBeamDir() const {
    return tile_beam_direction_;
  }
  const casacore::MDirection tileBeamDirCopy() const {
    return copyMeasure(casacore::MeasureHolder(tile_beam_direction_))
        .asMDirection();
  }
  const std::vector<double>& chanFreqs(std::size_t baseline = 0) const {
    return channel_frequencies_[baseline];
  }
  const std::vector<std::vector<double>>& BdaChanFreqs() const {
    return channel_frequencies_;
  }
  const std::vector<double>& chanWidths(std::size_t baseline = 0) const {
    return channel_widths_[baseline];
  }
  const std::vector<std::vector<double>>& BdaChanWidths() const {
    return channel_widths_;
  }
  const std::vector<double>& resolutions(std::size_t baseline = 0) const {
    return resolutions_[baseline];
  }
  const std::vector<double>& effectiveBW(std::size_t baseline = 0) const {
    return effective_bandwidth_[baseline];
  }
  bool hasBDAChannels() const {
    return channel_frequencies_.size() == nbaselines();
  }
  const std::string& getDataColName() const { return data_column_name_; }
  const std::string& getFlagColName() const { return flag_column_name_; }
  const std::string& getWeightColName() const { return weight_column_name_; }
  double totalBW() const { return total_bandwidth_; }
  double refFreq() const { return reference_frequency_; }
  ///@}

  /// Get the antenna numbers actually used in the (selected) baselines.
  /// E.g. [0,2,5,6]
  const std::vector<int>& antennaUsed() const { return antennas_used_; }

  /// Get the indices of all antennae in the used antenna vector above.
  /// -1 means that the antenna is not used.
  /// E.g. [0,-1,1,-1,-1,2,3] for the example above.
  const std::vector<int>& antennaMap() const { return antenna_map_; }

  /// Is the visibility data needed?
  bool needVisData() const { return need_data_; }
  /// Does the last step need to write data and/or flags?
  bool needWrite() const {
    return write_data_ || write_flags_ || write_weights_;
  }
  bool writeData() const { return write_data_; }
  bool writeFlags() const { return write_flags_; }
  bool writeWeights() const { return write_weights_; }
  /// Has the meta data been changed in a step (precluding an update)?
  bool metaChanged() const { return meta_changed_; }

  /// Set if visibility data needs to be read.
  void setNeedVisData() { need_data_ = true; }
  /// Set if data needs to be written.
  void setWriteData() { write_data_ = true; }
  void setWriteFlags() { write_flags_ = true; }
  void setWriteWeights() { write_weights_ = true; }
  /// Clear all write flags.
  void clearWrites() { write_data_ = write_flags_ = write_weights_ = false; }
  /// Set change of meta data.
  void setMetaChanged() { meta_changed_ = true; }
  void clearMetaChanged() { meta_changed_ = false; }

  /// Get the baseline table index of the autocorrelations.
  /// A negative value means there are no autocorrelations for that antenna.
  const std::vector<int>& getAutoCorrIndex() const;

  /// Get the lengths of the baselines (in meters).
  const std::vector<double>& getBaselineLengths() const;

  void setNThreads(unsigned int nThreads) { n_threads_ = nThreads; }

  unsigned int nThreads() const { return n_threads_; }

  /// Convert to a Record.
  /// The names of the fields in the record are the data names without 'its'.
  casacore::Record toRecord() const;

  /// Update the DPInfo object from a Record.
  /// It is possible that only a few fields are defined in the record.
  void fromRecord(const casacore::Record& rec);

  everybeam::CorrectionMode beamCorrectionMode() const {
    return beam_correction_mode_;
  };
  void setBeamCorrectionMode(everybeam::CorrectionMode mode) {
    beam_correction_mode_ = mode;
  }

  const casacore::MDirection& beamCorrectionDir() const {
    return beam_correction_direction_;
  }
  void setBeamCorrectionDir(const casacore::MDirection& dir) {
    beam_correction_direction_ = dir;
  }

  /// Determine if the channels have a regular layout.
  bool channelsAreRegular() const;

 private:
  /// Set which antennae are actually used.
  void setAntUsed();

  /// Creates a real copy of a casacore::Measure by exporting to a Record
  static casacore::MeasureHolder copyMeasure(
      const casacore::MeasureHolder fromMeas);

  bool need_data_;      ///< Are the visibility data needed?
  bool write_data_;     ///< Must the data be written?
  bool write_flags_;    ///< Must the flags be written?
  bool write_weights_;  ///< Must the weights be written?
  bool meta_changed_;   ///< Are meta data changed? (e.g., by averaging)
  std::string ms_name_;
  std::string data_column_name_;
  std::string flag_column_name_;
  std::string weight_column_name_;
  std::string antenna_set_;
  unsigned int n_correlations_;
  unsigned int start_channel_;
  unsigned int original_n_channels_;
  unsigned int n_channels_;
  unsigned int channel_averaging_factor_;
  unsigned int n_times_;
  std::vector<unsigned int> time_averaging_factors_;
  double start_time_;
  double time_interval_;
  bool bda_interval_factor_is_integer_;  // INTEGER_INTERVAL_FACTORS in
                                         // [BDA_TIME_AXIS}
  casacore::MDirection phase_center_;
  bool phase_center_is_original_;
  casacore::MDirection delay_center_;
  casacore::MDirection tile_beam_direction_;
  everybeam::CorrectionMode beam_correction_mode_;
  casacore::MDirection beam_correction_direction_;
  casacore::MPosition array_position_;
  /// Channels may differ between baselines when using BDA.
  /// If all baselines have equal channels, the outer vector only holds one
  /// inner vector. When using BDA, each baseline has its own inner vector.
  std::vector<std::vector<double>> channel_frequencies_;
  std::vector<std::vector<double>> channel_widths_;
  std::vector<std::vector<double>> resolutions_;
  std::vector<std::vector<double>> effective_bandwidth_;
  double total_bandwidth_;
  double reference_frequency_;
  casacore::Vector<casacore::String> antenna_names_;
  casacore::Vector<casacore::Double> antenna_diameters_;
  std::vector<casacore::MPosition> antenna_positions_;
  std::vector<int> antennas_used_;
  std::vector<int> antenna_map_;
  /// For each baseseline, the index of the first antenna.
  std::vector<std::size_t> antenna1_;
  /// For each baseseline, the index of the second antenna.
  std::vector<std::size_t> antenna2_;
  mutable std::vector<double> baseline_lengths_;
  /// For each antenna, the auto correlation index.
  mutable std::vector<int> auto_correlation_indices_;
  unsigned int n_threads_;
};

}  // namespace base
}  // namespace dp3

#endif
