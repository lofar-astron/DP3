// Upsample.cc: DPPP step class to Upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Upsample.h"

#include <iostream>

#include "../base/UVWCalculator.h"
#include "../common/ParameterSet.h"

#include <casacore/casa/BasicMath/Math.h>       // nearAbs
#include <casacore/casa/Arrays/ArrayLogical.h>  // anyTrue

#include <boost/make_unique.hpp>

#include <iomanip>
#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using dp3::base::DPBuffer;
using dp3::base::DPInfo;

namespace dp3 {
namespace steps {

Upsample::Upsample(const common::ParameterSet& parset, const string& prefix)
    : name_(prefix),
      old_time_interval_(0),
      time_step_(parset.getInt(prefix + "timestep")),
      update_uvw_(parset.getBool(prefix + "updateuvw", false)),
      prev_buffers_(),
      buffers_(time_step_),
      first_to_flush_(0),
      uvw_calculator_(),
      timer_() {}

Upsample::~Upsample() {}

void Upsample::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);
  info().setNeedVisData();
  info().setWriteData();

  old_time_interval_ = info().timeInterval();
  info().setTimeInterval(old_time_interval_ / time_step_);

  info().setMetaChanged();
  if (update_uvw_) {
    uvw_calculator_ = boost::make_unique<base::UVWCalculator>(
        info().phaseCenter(), info().arrayPos(), info().antennaPos());
  }
}

void Upsample::show(std::ostream& os) const {
  os << "Upsample " << name_ << '\n';
  os << "  time step : " << time_step_ << '\n';
}

bool Upsample::process(const DPBuffer& bufin) {
  double time0 = bufin.getTime() - 0.5 * old_time_interval_;
  double exposure = bufin.getExposure() / time_step_;

  // Duplicate the input buffer time_step_ times
  for (unsigned int i = 0; i < time_step_; ++i) {
    buffers_[i].copy(bufin);
    // Update the time centroid and time exposure
    const double time = time0 + info().timeInterval() * (i + 0.5);
    buffers_[i].setTime(time);
    buffers_[i].setExposure(exposure);

    if (update_uvw_) {
      casacore::Matrix<double>& buffer_uvw = buffers_[i].getUVW();
      if (!buffer_uvw.empty()) {
        assert(buffer_uvw.size() == 3 * info().nbaselines());
        double* uvw_ptr = buffer_uvw.data();
        for (std::size_t bl = 0; bl < info().nbaselines(); ++bl) {
          const casacore::Vector<double> uvw = uvw_calculator_->getUVW(
              info().getAnt1()[bl], info().getAnt2()[bl], time);
          std::copy_n(uvw.data(), 3, uvw_ptr);
          uvw_ptr += 3;
        }
      }
    }
  }

  if (prev_buffers_.empty()) {
    // First time slot, ask for next time slot first
    prev_buffers_.resize(time_step_);
    for (unsigned int i = 0; i < time_step_; ++i) {
      prev_buffers_[i].copy(buffers_[i]);  // No shallow copy
    }
    return false;
  }

  // Flush the prev_buffers_. Skip parts at the beginning as determined
  // in the previous call to process. Skip parts at the end as determined
  // now. Also, determine at which step the next buffer should start to
  // flush in the next call of process.
  unsigned int curIndex = 0;  // Index in the current buffers
  for (unsigned int prevIndex = first_to_flush_; prevIndex < time_step_;
       prevIndex++) {
    first_to_flush_ = 0;  // reset for next use
    // Advance curIndex until
    // buffers[curIndex].time >= prevBuffers[prevIndex].time
    while (prev_buffers_[prevIndex].getTime() > buffers_[curIndex].getTime() &&
           !casacore::nearAbs(prev_buffers_[prevIndex].getTime(),
                              buffers_[curIndex].getTime(),
                              0.4 * info().timeInterval())) {
      curIndex++;
    }
    if (casacore::nearAbs(prev_buffers_[prevIndex].getTime(),
                          buffers_[curIndex].getTime(),
                          0.4 * info().timeInterval())) {
      // Found double buffer, choose which one to use
      // If both totally flagged, prefer prevbuffer
      if (allTrue(buffers_[curIndex].getFlags())) {
        // Use prevBuffer
        first_to_flush_ = curIndex + 1;
        getNextStep()->process(prev_buffers_[prevIndex]);
      } else {
        // Use next buffer; stop processing prevbuffer.
        // This will uncorrectly give flagged if data has been flagged
        // and a time slot has been inserted and time_step_ > 2.
        break;
      }
    } else {
      // No double buffer, just flush the prevbuffer
      getNextStep()->process(prev_buffers_[prevIndex]);
    }
  }

  prev_buffers_.swap(buffers_);  // buffers_ will be overwritten later

  return false;
}

void Upsample::finish() {
  // Flush prev_buffers_
  for (unsigned int i = first_to_flush_; i < time_step_; ++i) {
    getNextStep()->process(prev_buffers_[i]);
  }

  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
