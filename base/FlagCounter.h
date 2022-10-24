// FlagCounter.h: Class to keep counts of nr of flagged points
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Class to keep counts of nr of flagged points
/// @author Ger van Diepen

#ifndef DPPP_FLAGCOUNTER_H
#define DPPP_FLAGCOUNTER_H

#include <casacore/casa/Arrays/Vector.h>

#include <cstdint>
#include <ostream>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace base {
class DPInfo;

/// @brief Class to keep counts of nr of flagged points

/// This class contains counts the number of flags.
/// The flags can be counted per baseline, channel, and correlation.
/// Once the counting is completed, they can be printed using the 'show'
/// functions. When printing, the baselines counts are shown per antenna.
/// Optionally the flagging percentages can be saved in a table.
/// The name of the table is the MS name suffixed by the step name and
/// '.flagxx'.

class FlagCounter {
 public:
  /// The default constructor creates an empty object. It does not save.
  FlagCounter() = default;

  /// This constructor creates an empty object.
  /// It reads info from the parset to see if percentages have to be saved.
  FlagCounter(const common::ParameterSet&, const std::string& prefix);

  /// Size all counters and initialize them to zero using the sizes
  /// from the DPInfo object.
  void init(const DPInfo& info);

  /// Increment the count per baseline.
  void incrBaseline(unsigned int bl) { base_line_counts_[bl]++; }

  /// Increment the count per channel.
  void incrChannel(unsigned int chan) { channel_counts_[chan]++; }

  /// Increment the count per correlation.
  void incrCorrelation(unsigned int corr) { correlation_counts_[corr]++; }

  /// Add the contents of that to this.
  void add(const FlagCounter& that);

  /// Get the counts.
  ///@{
  const std::vector<int64_t>& baselineCounts() const {
    return base_line_counts_;
  }
  const std::vector<int64_t>& channelCounts() const { return channel_counts_; }
  const std::vector<int64_t>& correlationCounts() const {
    return correlation_counts_;
  }
  ///@}

  /// Print the counts and optionally save percentages in a table.
  void showBaseline(std::ostream& os, int64_t ntimes) const;
  void showChannel(std::ostream& os, int64_t ntimes) const;
  void showCorrelation(std::ostream& os, int64_t ntimes) const;

  //// Print ratio of flagged visibilities per antenna.
  void showStation(std::ostream& os, int64_t ntimes) const;

  /// Show percentage with 1 decimal.
  static void showPerc1(std::ostream&, double value, double total);

  /// Show percentage with 3 decimals.
  static void showPerc3(std::ostream&, double value, double total);

 private:
  /// Save the percentages per station in a table.
  void saveStation(int64_t npoints, const casacore::Vector<int64_t>& nused,
                   const casacore::Vector<int64_t>& count) const;

  /// Save the percentages per channel.
  void saveChannel(int64_t npoints, const std::vector<int64_t>& count) const;

  const DPInfo* info_{nullptr};
  std::string save_filename_{};
  double warning_percentage_{0.0};
  bool show_fully_flagged_{false};
  bool save_{false};
  std::string path_{};
  std::string name_{};
  std::vector<int64_t> base_line_counts_{};
  std::vector<int64_t> channel_counts_{};
  std::vector<int64_t> correlation_counts_{};
};

}  // namespace base
}  // namespace dp3

#endif
