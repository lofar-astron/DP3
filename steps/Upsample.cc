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

Upsample::Upsample(const common::ParameterSet& parset,
                   const std::string& prefix)
    : Upsample(prefix, parset.getUint(prefix + "timestep"),
               parset.getBool(prefix + "updateuvw", false)) {}

Upsample::Upsample(const std::string& name, unsigned int time_step,
                   bool update_uvw)
    : name_(name),
      time_step_(time_step),
      update_uvw_(update_uvw),
      prev_buffers_(),
      buffers_(time_step),
      first_to_flush_(0),
      uvw_calculator_(),
      timer_() {
  if (time_step_ <= 1) {
    throw std::invalid_argument("Upsample: Invalid time step value: " +
                                std::to_string(time_step_));
  }
}

Upsample::~Upsample() {}

void Upsample::updateInfo(const DPInfo& info_in) {
  Step::updateInfo(info_in);
  info().setNeedVisData();
  info().setWriteData();
  info().setTimeIntervalAndSteps(info().timeInterval() / time_step_,
                                 info().ntime() * time_step_);
  info().setMetaChanged();

  if (update_uvw_) {
    uvw_calculator_ = boost::make_unique<base::UVWCalculator>(
        info().phaseCenter(), info().arrayPos(), info().antennaPos());
  }
}

void Upsample::show(std::ostream& os) const {
  os << "Upsample " << name_ << '\n';
  os << "  time step   : " << time_step_ << '\n';
  os << "  update UVW  : " << std::boolalpha << update_uvw_ << '\n';
}

bool Upsample::process(const DPBuffer& bufin) {
  double time0 = bufin.getTime() - 0.5 * info().timeInterval() * time_step_;
  double exposure = bufin.getExposure() / time_step_;

  // Duplicate the input buffer time_step_ times
  for (unsigned int i = 0; i < time_step_; ++i) {
    buffers_[i].copy(bufin);
    // Update the time centroid and time exposure
    const double time = time0 + info().timeInterval() * (i + 0.5);
    buffers_[i].setTime(time);
    buffers_[i].setExposure(exposure);

    if (update_uvw_) {
      buffers_[i].getUVW().resize(3, info().nbaselines());
      double* uvw_ptr = buffers_[i].getUVW().data();
      for (std::size_t bl = 0; bl < info().nbaselines(); ++bl) {
        const std::array<double, 3> uvw = uvw_calculator_->getUVW(
            info().getAnt1()[bl], info().getAnt2()[bl], time);
        std::copy_n(uvw.data(), 3, uvw_ptr);
        uvw_ptr += 3;
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
