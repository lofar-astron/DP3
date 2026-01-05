// Upsample.cc: DPPP step class to Upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Upsample.h"

#include <iostream>

#include "base/UVWCalculator.h"
#include "common/ParameterSet.h"

#include <casacore/casa/BasicMath/Math.h>  // nearAbs

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
  GetWritableInfoOut().setTimeIntervalAndSteps(
      getInfoOut().timeInterval() / time_step_,
      getInfoOut().ntime() * time_step_);
  GetWritableInfoOut().setMetaChanged();

  if (update_uvw_) {
    uvw_calculator_ = std::make_unique<base::UVWCalculator>(
        getInfoOut().phaseCenter(), getInfoOut().arrayPos(),
        getInfoOut().antennaPos());
  }
}

void Upsample::show(std::ostream& os) const {
  os << "Upsample " << name_ << '\n';
  os << "  time step   : " << time_step_ << '\n';
  os << "  update UVW  : " << std::boolalpha << update_uvw_ << '\n';
}

void Upsample::UpdateTimeCentroidExposureAndUvw(
    std::unique_ptr<base::DPBuffer>& buffer, double time, double exposure) {
  buffer->SetTime(time);
  buffer->SetExposure(exposure);

  if (update_uvw_) {
    buffer->GetUvw().resize({getInfoOut().nbaselines(), 3});

    double* uvw_ptr = buffer->GetUvw().data();
    for (std::size_t bl = 0; bl < getInfoOut().nbaselines(); ++bl) {
      const std::array<double, 3> uvw = uvw_calculator_->getUVW(
          getInfoOut().getAnt1()[bl], getInfoOut().getAnt2()[bl], time);
      std::copy_n(uvw.data(), 3, uvw_ptr);
      uvw_ptr += 3;
    }
  }
}

bool Upsample::process(std::unique_ptr<base::DPBuffer> buffer) {
  double time0 =
      buffer->GetTime() - 0.5 * getInfoOut().timeInterval() * time_step_;
  double exposure = buffer->GetExposure() / time_step_;

  // Make a copy of the buffer for each time_step_
  // As an optimisation we re-use the buffer for the last time_step_
  // reducing the number of copies required by 1
  for (unsigned int i = 0; i < time_step_ - 1; ++i) {
    buffers_[i] = std::make_unique<base::DPBuffer>(*buffer);
    const double time = time0 + getInfoOut().timeInterval() * (i + 0.5);
    UpdateTimeCentroidExposureAndUvw(buffers_[i], time, exposure);
  }
  buffers_[time_step_ - 1] = std::move(buffer);
  UpdateTimeCentroidExposureAndUvw(
      buffers_[time_step_ - 1],
      time0 + getInfoOut().timeInterval() * (time_step_ - 1 + 0.5), exposure);

  if (prev_buffers_.empty()) {
    // First time slot, ask for next time slot first
    prev_buffers_.resize(time_step_);
    for (unsigned int i = 0; i < time_step_; ++i) {
      // No shallow copy
      prev_buffers_[i] = std::make_unique<base::DPBuffer>(*buffers_[i]);
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
    while (prev_buffers_[prevIndex]->GetTime() >
               buffers_[curIndex]->GetTime() &&
           !casacore::nearAbs(prev_buffers_[prevIndex]->GetTime(),
                              buffers_[curIndex]->GetTime(),
                              0.4 * getInfoOut().timeInterval())) {
      curIndex++;
    }
    if (casacore::nearAbs(prev_buffers_[prevIndex]->GetTime(),
                          buffers_[curIndex]->GetTime(),
                          0.4 * getInfoOut().timeInterval())) {
      // Found double buffer, choose which one to use
      // If both totally flagged, prefer prevbuffer
      if (xt::all(xt::equal(buffers_[curIndex]->GetFlags(), true))) {
        // Use prevBuffer
        first_to_flush_ = curIndex + 1;
        getNextStep()->process(std::move(prev_buffers_[prevIndex]));
      } else {
        // Use next buffer; stop processing prevbuffer.
        // This will uncorrectly give flagged if data has been flagged
        // and a time slot has been inserted and time_step_ > 2.
        break;
      }
    } else {
      // No double buffer, just flush the prevbuffer
      getNextStep()->process(std::move(prev_buffers_[prevIndex]));
    }
  }

  prev_buffers_.swap(buffers_);  // buffers_ will be overwritten later

  return false;
}

void Upsample::finish() {
  // Flush prev_buffers_
  for (unsigned int i = first_to_flush_; i < time_step_; ++i) {
    getNextStep()->process(std::move(prev_buffers_[i]));
  }

  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
