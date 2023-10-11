// DPInfo.h: General info about DP3 data processing attributes like averaging
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief General info about DP3 data processing attributes like averaging
/// @author Ger van Diepen

#ifndef DP3_BASE_DPINFO_H_
#define DP3_BASE_DPINFO_H_

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MeasureHolder.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Record.h>

#include "Direction.h"

namespace dp3 {
namespace base {

/// @brief General info about DP3 data processing attributes like averaging

/// This class contains the information about the number of correlations,
/// channels, baselines, and times.
/// It is initialized by the first step and updated by steps like
/// Averager that change the number of channels or times.
/// Steps can take information from it to know about shapes.

class DPInfo {
 public:
  explicit DPInfo(unsigned int n_correlations = 0,
                  unsigned int n_original_channels = 0,
                  unsigned int start_channel = 0,
                  std::string antenna_set = std::string());

  /// Set time information and derive the number of time slots.
  /// @param first_time Centroid time of the first time slot, in mjd seconds.
  /// @param last_time Centroid time of the last time slot, in mjd seconds.
  /// @param time_interval Time interval between two time slots, in mjd seconds.
  void setTimes(double first_time, double last_time, double time_interval);

  void setMsName(const std::string& ms_name) { ms_name_ = ms_name; }

  void setMsNames(const std::string& ms_name,
                  const std::string& data_column_name,
                  const std::string& flag_column_name,
                  const std::string& weight_column_name);

  void setIsBDAIntervalFactorInteger(bool isIntervalInteger) {
    bda_interval_factor_is_integer_ = isIntervalInteger;
  }

  /// Set the time interval and the number of time steps.
  void setTimeIntervalAndSteps(double timeInterval, unsigned int ntime) {
    time_interval_ = timeInterval;
    n_times_ = ntime;
  }

  /// Set the frequency info.
  /// An empty resolutions or effectiveBW is default to chanWidths.
  /// total_bandwidth_ is set to the sum of effectiveBW.
  /// If refFreq is 0, it is set to the middle of chanFreqs (mean if even).
  /// @throw std::invalid argument When vector sizes do not match.
  void setChannels(std::vector<double>&& chanFreqs,
                   std::vector<double>&& chanWidths,
                   std::vector<double>&& resolutions = std::vector<double>(),
                   std::vector<double>&& effectiveBW = std::vector<double>(),
                   double refFreq = 0, int spectralWindow = 0);

  /// Set the frequency info, using different info per baseline.
  /// total_bandwidth_ is set to the sum of effectiveBW, which should be equal
  /// for all baselines. If refFreq is 0, it is set to the middle of chanFreqs
  /// (mean if even). of the baseline with the most channels.
  /// @param chanFreqs Channel frequencies for each baseline in Hz.
  /// @param chanWidths Channel widths for each baseline in Hz.
  /// @param resolutions Channel resolution for each baseline in Hz.
  ///        If omitted/empty, channel widths are used.
  /// @param effectiveBW Effective bandwidth for the channels in each baseline
  //         in Hz. If omitted/empty, channel widths are used.
  /// @throw std::invalid argument When a vector has an incorrect size.
  ///        The size of the vector arguments should be the number of baselines.
  ///        For each baseline, the number of elements in the vectors for that
  ///        baseline should be equal.
  void setChannels(std::vector<std::vector<double>>&& chanFreqs,
                   std::vector<std::vector<double>>&& chanWidths,
                   std::vector<std::vector<double>>&& resolutions =
                       std::vector<std::vector<double>>(),
                   std::vector<std::vector<double>>&& effectiveBW =
                       std::vector<std::vector<double>>(),
                   double refFreq = 0, int spectralWindow = 0);

  void setArrayInformation(const casacore::MPosition& arrayPos,
                           const casacore::MDirection& phaseCenter,
                           const casacore::MDirection& delayCenter,
                           const casacore::MDirection& tileBeamDir);

  /// Set the info for the given antennae and baselines.
  void setAntennas(const std::vector<std::string>& antNames,
                   const std::vector<double>& antDiam,
                   const std::vector<casacore::MPosition>& antPos,
                   const std::vector<int>& ant1, const std::vector<int>& ant2);

  /// Update the info for the given average factors.
  /// If chanAvg is higher than the actual nr of channels, it is reset.
  /// The same is true for timeAvg.
  /// It returns the possibly reset nr of channels to average.
  unsigned int update(unsigned int chanAvg, unsigned int timeAvg);

  /// Update the info from the given selection parameters.
  /// Optionally unused stations are really removed from the antenna lists.
  ///
  /// @pre The range [startChan, startChan + nchan) does not exceed the ranges
  ///      of the interal vector of the class.
  ///
  /// @param startChan The first channal to use.
  /// @param nchan The number of channels to use.
  void update(unsigned int startChan, unsigned int nchan,
              const std::vector<unsigned int>& baselines, bool remove);

  /// Update the info for the given average factors.
  void update(std::vector<unsigned int>&& timeAvg);

  /// Remove unused stations from the antenna lists.
  void removeUnusedAnt();

  void setPhaseCenter(const casacore::MDirection& phase_center) {
    phase_center_ = phase_center;
  }

  /// Get the info.
  ///@{
  const std::string& msName() const { return ms_name_; }
  const std::string& dataColumnName() const { return data_column_name_; }
  const std::string& flagColumnName() const { return flag_column_name_; }
  const std::string& weightColumnName() const { return weight_column_name_; }
  const std::string& antennaSet() const { return antenna_set_; }
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
  double startTime() const { return first_time_ - 0.5 * time_interval_; }
  double firstTime() const { return first_time_; }
  double lastTime() const { return last_time_; }
  double timeInterval() const { return time_interval_; }
  bool isBDAIntervalFactorInteger() const {
    return bda_interval_factor_is_integer_;
  }
  const std::vector<int>& getAnt1() const { return antenna1_; }
  const std::vector<int>& getAnt2() const { return antenna2_; }
  const std::vector<std::string>& antennaNames() const {
    return antenna_names_;
  }
  const std::vector<double>& antennaDiam() const { return antenna_diameters_; }
  const std::vector<casacore::MPosition>& antennaPos() const {
    return antenna_positions_;
  }
  const casacore::MPosition& arrayPos() const { return array_position_; }
  const casacore::MPosition arrayPosCopy() const {
    return copyMeasure(casacore::MeasureHolder(array_position_)).asMPosition();
  }
  const casacore::MDirection& originalPhaseCenter() const {
    return original_phase_center_;
  }
  const casacore::MDirection& phaseCenter() const { return phase_center_; }
  const casacore::MDirection phaseCenterCopy() const {
    return copyMeasure(casacore::MeasureHolder(phase_center_)).asMDirection();
  }
  /// @return The ra/dec direction corresponding to the phase center.
  Direction phaseCenterDirection() const;
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
  double totalBW() const { return total_bandwidth_; }
  double refFreq() const { return reference_frequency_; }
  int spectralWindow() const { return spectral_window_; }
  ///@}

  /// Get the antenna numbers actually used in the (selected) baselines.
  /// E.g. [0,2,5,6]
  const std::vector<int>& antennaUsed() const { return antennas_used_; }

  /// Get the indices of all antennae in the used antenna vector above.
  /// -1 means that the antenna is not used.
  /// E.g. [0,-1,1,-1,-1,2,3] for the example above.
  const std::vector<int>& antennaMap() const { return antenna_map_; }

  /// Has the meta data been changed in a step (precluding an update)?
  bool metaChanged() const { return meta_changed_; }

  /// Set change of meta data.
  void setMetaChanged() { meta_changed_ = true; }
  void clearMetaChanged() { meta_changed_ = false; }

  /// Get the baseline table index of the autocorrelations.
  /// A negative value means there are no autocorrelations for that antenna.
  const std::vector<int>& getAutoCorrIndex() const;

  /// Get the lengths of the baselines (in meters).
  const std::vector<double>& getBaselineLengths() const;

  /// Get the beam correction mode.
  /// @return An integer representation of an everybeam::CorrectionMode value.
  int beamCorrectionMode() const { return beam_correction_mode_; }

  /// Set the beam correction mode.
  /// @param mode An integer representation of an everybeam::CorrectionMode
  /// value.
  void setBeamCorrectionMode(int mode) { beam_correction_mode_ = mode; }

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

  bool meta_changed_;  ///< Are meta data changed? (e.g., by averaging)
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
  std::vector<unsigned int> time_averaging_factors_;
  double
      first_time_;    ///< Centroid time of the first time slot, in mjd seconds.
  double last_time_;  ///< Centroid time of the last time slot, in mjd seconds.
  double time_interval_;  ///< Time interval between two time slots, in mjd
                          ///< seconds.
  unsigned int n_times_;  ///< Number of time slots.
  bool bda_interval_factor_is_integer_{false};  // INTEGER_INTERVAL_FACTORS in
                                                // [BDA_TIME_AXIS}
  casacore::MDirection original_phase_center_;
  casacore::MDirection phase_center_;
  casacore::MDirection delay_center_;
  casacore::MDirection tile_beam_direction_;
  /// Correction mode for EveryBeam. Since DPInfo is part of the public DP3
  /// interface, using an integer avoids a public dependency on EveryBeam.
  int beam_correction_mode_;
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
  int spectral_window_;
  std::vector<std::string> antenna_names_;
  std::vector<double> antenna_diameters_;
  std::vector<casacore::MPosition> antenna_positions_;
  std::vector<int> antennas_used_;
  std::vector<int> antenna_map_;
  /// For each baseseline, the index of the first antenna.
  std::vector<int> antenna1_;
  /// For each baseseline, the index of the second antenna.
  std::vector<int> antenna2_;
  mutable std::vector<double> baseline_lengths_;
  /// For each antenna, the auto correlation index.
  mutable std::vector<int> auto_correlation_indices_;
};

}  // namespace base
}  // namespace dp3

#endif  // DP3_BASE_DPINFO_H_
