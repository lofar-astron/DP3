// DPBuffer.cc: Buffer holding the data of a timeslot/band
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "DPBuffer.h"

#include <casacore/casa/version.h>

// Casacore < 3.4 does not support move semantics for casacore::Array
// and uses reference semantics in the copy constructor.
#if CASACORE_MAJOR_VERSION > 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION >= 4)
#define USE_CASACORE_MOVE_SEMANTICS
#endif

namespace dp3 {
namespace base {

DPBuffer::DPBuffer(double time, double exposure)
    : itsTime(time), itsExposure(exposure) {}

DPBuffer::DPBuffer(const DPBuffer& that) { operator=(that); }

DPBuffer::DPBuffer(DPBuffer&& that)
    : itsTime(that.itsTime),
      itsExposure(that.itsExposure),
      itsRowNrs(std::move(that.itsRowNrs)),
      itsData(std::move(that.itsData)),
      itsFlags(std::move(that.itsFlags)),
      itsUVW(std::move(that.itsUVW)),
      itsWeights(std::move(that.itsWeights)),
      itsFullResFlags(std::move(that.itsFullResFlags)) {
#ifndef USE_CASACORE_MOVE_SEMANTICS
  // The copy constructor for casacore::Array creates references. Since
  // moving a buffer does not have reference semantics, clear 'that'.
  that.itsRowNrs.assign(decltype(that.itsRowNrs)());
  that.itsData.assign(decltype(that.itsData)());
  that.itsFlags.assign(decltype(that.itsFlags)());
  that.itsUVW.assign(decltype(that.itsUVW)());
  that.itsWeights.assign(decltype(that.itsWeights)());
  that.itsFullResFlags.assign(decltype(that.itsFullResFlags)());
#endif
}

DPBuffer& DPBuffer::operator=(const DPBuffer& that) {
  if (this != &that) {
    itsTime = that.itsTime;
    itsExposure = that.itsExposure;
    itsRowNrs.reference(that.itsRowNrs);
    itsData.reference(that.itsData);
    itsFlags.reference(that.itsFlags);
    itsWeights.reference(that.itsWeights);
    itsUVW.reference(that.itsUVW);
    itsFullResFlags.reference(that.itsFullResFlags);
  }
  return *this;
}

DPBuffer& DPBuffer::operator=(DPBuffer&& that) {
  if (this != &that) {
    itsTime = that.itsTime;
    itsExposure = that.itsExposure;
    // Casacore < 3.4.0 does not support move semantics for casacore::Array.
    // The copy assignment operator for casacore::Array then creates copies.
#ifdef USE_CASACORE_MOVE_SEMANTICS
    itsRowNrs = std::move(that.itsRowNrs);
    itsData = std::move(that.itsData);
    itsFlags = std::move(that.itsFlags);
    itsWeights = std::move(that.itsWeights);
    itsUVW = std::move(that.itsUVW);
    itsFullResFlags = std::move(that.itsFullResFlags);
#else
    // Create references. Since move assignment does not use reference
    // semantics, clear 'that'.
    itsRowNrs.reference(that.itsRowNrs);
    itsData.reference(that.itsData);
    itsFlags.reference(that.itsFlags);
    itsWeights.reference(that.itsWeights);
    itsUVW.reference(that.itsUVW);
    itsFullResFlags.reference(that.itsFullResFlags);
    that.itsRowNrs.assign(decltype(that.itsRowNrs)());
    that.itsData.assign(decltype(that.itsData)());
    that.itsFlags.assign(decltype(that.itsFlags)());
    that.itsUVW.assign(decltype(that.itsUVW)());
    that.itsWeights.assign(decltype(that.itsWeights)());
    that.itsFullResFlags.assign(decltype(that.itsFullResFlags)());
#endif
  }
  return *this;
}

void DPBuffer::copy(const DPBuffer& that) {
  if (this != &that) {
    itsTime = that.itsTime;
    itsExposure = that.itsExposure;
    itsRowNrs.assign(that.itsRowNrs);
    if (!that.itsData.empty()) {
      itsData.assign(that.itsData);
    }
    if (!that.itsFlags.empty()) {
      itsFlags.assign(that.itsFlags);
    }
    if (!that.itsWeights.empty()) {
      itsWeights.assign(that.itsWeights);
    }
    if (!that.itsUVW.empty()) {
      itsUVW.assign(that.itsUVW);
    }
    if (!that.itsFullResFlags.empty()) {
      itsFullResFlags.assign(that.itsFullResFlags);
    }
  }
}

void DPBuffer::referenceFilled(const DPBuffer& that) {
  if (this != &that) {
    itsTime = that.itsTime;
    itsExposure = that.itsExposure;
    itsRowNrs.reference(that.itsRowNrs);
    if (!that.itsData.empty()) {
      itsData.reference(that.itsData);
    }
    if (!that.itsFlags.empty()) {
      itsFlags.reference(that.itsFlags);
    }
    if (!that.itsWeights.empty()) {
      itsWeights.reference(that.itsWeights);
    }
    if (!that.itsUVW.empty()) {
      itsUVW.reference(that.itsUVW);
    }
    if (!that.itsFullResFlags.empty()) {
      itsFullResFlags.reference(that.itsFullResFlags);
    }
  }
}

void DPBuffer::mergeFullResFlags(casacore::Cube<bool>& fullResFlags,
                                 const casacore::Cube<bool>& flags) {
  // Flag shape is [ncorr, newnchan, nbl].
  // FullRes shape is [orignchan, navgtime, nbl]
  // where orignchan = navgchan * newnchan.
  const casacore::IPosition& fullResShape = fullResFlags.shape();
  const casacore::IPosition& flagShape = flags.shape();
  int orignchan = fullResShape[0];
  int newnchan = flagShape[1];
  int navgchan = orignchan / newnchan;
  int navgtime = fullResShape[1];
  int nbl = fullResShape[2];
  int ncorr = flagShape[0];
  bool* fullResPtr = fullResFlags.data();
  const bool* flagPtr = flags.data();
  // Loop over all baselines and new channels.
  // Only use the first correlation in the loop.
  for (int j = 0; j < nbl; ++j) {
    for (int i = 0; i < newnchan; ++i) {
      // If ta data point is flagged, the flags in the corresponding
      // FullRes window have to be set.
      // This is needed in case a data point is further averaged.
      if (*flagPtr) {
        for (int i = 0; i < navgtime; ++i) {
          std::fill(fullResPtr, fullResPtr + navgchan, true);
          fullResPtr += orignchan;
        }
        fullResPtr -= navgtime * orignchan;
      }
      flagPtr += ncorr;
      fullResPtr += navgchan;
    }
    // Set pointer to next baseline.
    fullResPtr += (navgtime - 1) * orignchan;
  }
}

}  // namespace base
}  // namespace dp3
