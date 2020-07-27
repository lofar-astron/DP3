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

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <deque>
#include <string>

#include "../DPPP/DPInput.h"
#include "../DPPP/DPBuffer.h"

#include "../Common/ParameterSet.h"

#include <aocommon/lane.h>

#include <casacore/casa/Arrays/Cube.h>

extern "C" void register_interpolate();

namespace DP3 {
namespace DPPP {

class Interpolate : public DPStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Interpolate(DPInput*, const ParameterSet&, const string& prefix);

  virtual ~Interpolate() = default;

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  static DPStep::ShPtr makeStep(DPInput* input, const ParameterSet& parset,
                                const std::string& prefix);

 private:
  void interpolateTimestep(size_t index);
  void interpolateSample(size_t timestep, size_t baseline, size_t channel,
                         size_t pol);
  void sendFrontBufferToNextStep();
  void interpolationThread();

  struct Sample {
    Sample() = default;
    Sample(size_t timestep_, size_t baseline_, size_t channel_, size_t pol_)
        : timestep(timestep_),
          baseline(baseline_),
          channel(channel_),
          pol(pol_) {}
    size_t timestep;
    size_t baseline;
    size_t channel;
    size_t pol;
  };

  std::string _name;
  size_t _interpolatedPos;
  std::deque<DPBuffer> _buffers;
  size_t _windowSize;
  NSTimer _timer;
  aocommon::Lane<Sample> _lane;
  std::vector<float> _kernelLookup;
};

}  // namespace DPPP
}  // namespace DP3

#endif
