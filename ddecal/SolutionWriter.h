// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLUTIONWRITER_H
#define DP3_DDECAL_SOLUTIONWRITER_H

#include "constraints/Constraint.h"

#include "../base/CalType.h"
#include "../base/Direction.h"

#include <schaapcommon/h5parm/h5parm.h>

namespace dp3 {
namespace ddecal {

class SolutionWriter {
 public:
  /**
   * Constructor. Initializes the writer for writing to a specific file.
   * @param filename Name of H5Parm file. Will be overwritten if it exists.
   */
  explicit SolutionWriter(const std::string& filename);

  /**
   * Write main antenna properties.
   */
  void AddAntennas(
      const std::vector<std::string>& all_antenna_names,
      const std::vector<std::array<double, 3>>& all_antenna_positions);

  /**
   * (Over)write solutions to the H5Parm file.
   */
  void Write(
      const std::vector<std::vector<std::vector<std::complex<double>>>>&
          solutions,
      const std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>&
          constraint_solutions,
      double start_time, double solution_interval, base::CalType mode,
      const std::vector<std::string>& used_antenna_names,
      const std::vector<base::Direction>& source_directions,
      const std::vector<std::vector<std::string>>& directions,
      const std::vector<double>& chan_freqs,
      const std::vector<double>& chan_block_freqs, const std::string& history);

 private:
  schaapcommon::h5parm::H5Parm h5parm_;
};

}  // namespace ddecal
}  // namespace dp3

#endif