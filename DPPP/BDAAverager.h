// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Step for compressing regular data into BDA data.
/// @author Maik Nijhuis

#ifndef DPPP_BDAAVERAGER_H
#define DPPP_BDAAVERAGER_H

#include "DPStep.h"

#include <vector>

namespace DP3 {
namespace Common {
class ParameterSet;
}  // namespace Common
}  // namespace DP3

namespace DP3 {
namespace DPPP {

class BDABuffer;

class BDAAverager : public DPStep {
 public:
  BDAAverager();

  ~BDAAverager() override;

  bool process(const DPBuffer&) override;

  void finish() override;

  void show(std::ostream&) const override;

  void updateInfo(const DPInfo&) override;

 private:
  struct Baseline {
    Baseline(std::size_t factor, std::size_t n_channels,
             std::size_t n_correlations);
    std::size_t added;         ///< Number of added regular intervals.
    const std::size_t factor;  ///< Averaging factor for this baseline.
    std::vector<std::complex<float>> data;
    std::vector<float> weights;
    float summed_weight;
    double uvw[3];
  };

  rownr_t next_rownr_;
  std::size_t bda_pool_size_;
  std::unique_ptr<BDABuffer> bda_buffer_;
  std::vector<Baseline> baselines_;
};

}  // namespace DPPP
}  // namespace DP3

#endif  // DPPP_BDAAVERAGER_H
