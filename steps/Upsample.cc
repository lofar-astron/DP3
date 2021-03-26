// Upsample.cc: DPPP step class to Upsample visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "Upsample.h"

#include <iostream>

#include "../common/ParameterSet.h"

#include <casacore/casa/BasicMath/Math.h>       // nearAbs
#include <casacore/casa/Arrays/ArrayLogical.h>  // anyTrue

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

Upsample::Upsample(InputStep*, const common::ParameterSet& parset,
                   const string& prefix)
    : itsOldTimeInterval(0),
      itsTimeStep(parset.getInt(prefix + "timestep")),
      itsFirstToFlush(0) {
  itsBuffers.resize(itsTimeStep);
}

Upsample::~Upsample() {}

void Upsample::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();

  itsOldTimeInterval = info().timeInterval();
  info().setTimeInterval(itsOldTimeInterval / itsTimeStep);

  info().setMetaChanged();
}

void Upsample::show(std::ostream& os) const {
  os << "Upsample " << itsName << '\n';
  os << "  time step : " << itsTimeStep << '\n';
}

bool Upsample::process(const DPBuffer& bufin) {
  double time0 = bufin.getTime() - 0.5 * itsOldTimeInterval;
  double exposure = bufin.getExposure() / itsTimeStep;

  // Duplicate the input buffer itsTimeStep times
  for (unsigned int i = 0; i < itsTimeStep; ++i) {
    itsBuffers[i].copy(bufin);
    // Update the time centroid and time exposure
    itsBuffers[i].setTime(time0 + info().timeInterval() * (i + 0.5));
    itsBuffers[i].setExposure(exposure);
  }

  if (itsPrevBuffers.empty()) {
    // First time slot, ask for next time slot first
    itsPrevBuffers.resize(itsTimeStep);
    for (unsigned int i = 0; i < itsTimeStep; ++i) {
      itsPrevBuffers[i].copy(itsBuffers[i]);  // No shallow copy
    }
    return false;
  }

  // Flush the itsPrevBuffers. Skip parts at the beginning as determined
  // in the previous call to process. Skip parts at the end as determined
  // now. Also, determine at which step the next buffer should start to
  // flush in the next call of process.
  unsigned int curIndex = 0;  // Index in the current buffers
  for (unsigned int prevIndex = itsFirstToFlush; prevIndex < itsTimeStep;
       prevIndex++) {
    itsFirstToFlush = 0;  // reset for next use
    // Advance curIndex until
    // buffers[curIndex].time >= prevBuffers[prevIndex].time
    while (itsPrevBuffers[prevIndex].getTime() >
               itsBuffers[curIndex].getTime() &&
           !casacore::nearAbs(itsPrevBuffers[prevIndex].getTime(),
                              itsBuffers[curIndex].getTime(),
                              0.4 * info().timeInterval())) {
      curIndex++;
    }
    if (casacore::nearAbs(itsPrevBuffers[prevIndex].getTime(),
                          itsBuffers[curIndex].getTime(),
                          0.4 * info().timeInterval())) {
      // Found double buffer, choose which one to use
      // If both totally flagged, prefer prevbuffer
      if (allTrue(itsBuffers[curIndex].getFlags())) {
        // Use prevBuffer
        itsFirstToFlush = curIndex + 1;
        getNextStep()->process(itsPrevBuffers[prevIndex]);
      } else {
        // Use next buffer; stop processing prevbuffer.
        // This will uncorrectly give flagged if data has been flagged
        // and a time slot has been inserted and itsTimeStep > 2.
        break;
      }
    } else {
      // No double buffer, just flush the prevbuffer
      getNextStep()->process(itsPrevBuffers[prevIndex]);
    }
  }

  itsPrevBuffers.swap(itsBuffers);  // itsBuffers will be overwritten later

  return false;
}

void Upsample::finish() {
  // Flush itsPrevBuffers
  for (unsigned int i = itsFirstToFlush; i < itsTimeStep; ++i) {
    getNextStep()->process(itsPrevBuffers[i]);
  }

  // Let the next steps finish.
  getNextStep()->finish();
}
}  // namespace steps
}  // namespace dp3
