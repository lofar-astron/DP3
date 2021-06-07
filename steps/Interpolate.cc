// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Interpolate.h"

#include "../base/DPBuffer.h"
#include "../base/DPInfo.h"
#include "../base/DP3.h"

#include "../common/buffered_lane.h"
#include "../common/ParameterSet.h"
#include "../common/StringTools.h"

#include <casacore/casa/Arrays/ArrayMath.h>

#include <iostream>
#include <iomanip>
#include <thread>

using casacore::IPosition;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::base::FlagCounter;

namespace dp3 {
namespace steps {

Interpolate::Interpolate(InputStep* /*input*/,
                         const common::ParameterSet& parset,
                         const string& prefix)
    : _name(prefix),
      _interpolatedPos(0),
      _windowSize(parset.getUint(prefix + "windowsize", 15)) {
  if (_windowSize % 2 != 1)
    throw std::runtime_error(
        "Window size of Interpolate action should be an odd number");

  _kernelLookup.reserve(_windowSize * _windowSize);
  for (int t = 0; t != int(_windowSize); ++t) {
    int y = t - int(_windowSize / 2);
    for (int ch = 0; ch != int(_windowSize); ++ch) {
      int x = ch - int(_windowSize / 2);
      double windowDist = double(x * x + y * y);
      // Gaussian function with sigma = 1
      // (evaluated with double prec, then converted to floats)
      double w = std::exp(windowDist * -0.5);
      _kernelLookup.emplace_back(w);
    }
  }
}

void Interpolate::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;
  info().setNeedVisData();
  info().setWriteData();
  info().setWriteFlags();
}

void Interpolate::show(std::ostream& os) const {
  os << "Interpolate " << _name << '\n';
  os << "  windowsize:     " << _windowSize << '\n';
}

void Interpolate::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  FlagCounter::showPerc1(os, _timer.getElapsed(), duration);
  os << " Interpolate " << _name << '\n';
}

bool Interpolate::process(const DPBuffer& buf) {
  _timer.start();
  // Collect the data in buffers.
  _buffers.emplace_back();
  _buffers.back().copy(buf);
  // If we have a full window of data, interpolate everything
  // up to the middle of the window
  if (_buffers.size() >= _windowSize) {
    size_t mid = _windowSize / 2;
    while (_interpolatedPos <= mid) {
      interpolateTimestep(_interpolatedPos);
      ++_interpolatedPos;
    }
    // Buffers are only pushed to the next step when they are completely
    // out of the window. This is because flags need to be set to false,
    // however the flag information of the entire window is needed during
    // interpolation, so these can only be set to false after processing.
    sendFrontBufferToNextStep();
  }
  _timer.stop();
  return true;
}

void Interpolate::sendFrontBufferToNextStep() {
  casacore::IPosition shp = _buffers.front().getData().shape();
  size_t nPol = shp[0], nChan = shp[1], nBl = shp[2], n = nPol * nChan * nBl;
  // Set all flags to false
  bool* flags = _buffers.front().getFlags().data();
  casacore::Complex* data = _buffers.front().getData().data();
  std::fill(flags, flags + n, false);
  // Flag NaN values (values for which the entire window was flagged on input)
  for (size_t i = 0; i != n; ++i) {
    if (!std::isfinite(data[i].real()) || !std::isfinite(data[i].imag())) {
      // The datum value is also set to 0, because NaNs sometimes give problems
      // in certain software, even when they are flagged (e.g. in Sagecal).
      data[i] = 0.0;
      flags[i] = true;
    }
  }

  _timer.stop();
  getNextStep()->process(_buffers.front());
  _timer.start();

  _buffers.pop_front();
  --_interpolatedPos;
}

void Interpolate::finish() {
  _timer.start();

  // Interpolate everything up to the end of the window
  while (_interpolatedPos < _buffers.size()) {
    interpolateTimestep(_interpolatedPos);
    ++_interpolatedPos;
  }
  while (!_buffers.empty()) {
    sendFrontBufferToNextStep();
  }

  _timer.stop();

  getNextStep()->finish();
}

#define BUFFER_SIZE 1024

void Interpolate::interpolateTimestep(size_t index) {
  const IPosition shp = _buffers.front().getData().shape();
  const size_t nPol = shp[0], nChan = shp[1], nPerBl = nPol * nChan,
               nBl = shp[2];

  std::vector<std::thread> threads;
  size_t nthreads = std::min<size_t>(getInfo().nThreads(), 8);
  _lane.resize(nthreads * BUFFER_SIZE);
  common::lane_write_buffer<Sample> buflane(&_lane, BUFFER_SIZE);
  threads.reserve(nthreads);
  for (size_t i = 0; i != nthreads; ++i)
    threads.emplace_back(&Interpolate::interpolationThread, this);

  std::vector<casacore::Complex> dataBlock;
  for (size_t bl = 0; bl < nBl; ++bl) {
    bool* flags = _buffers[index].getFlags().data() + bl * nPerBl;
    for (size_t ch = 0; ch != nChan; ++ch) {
      for (size_t p = 0; p != nPol; ++p) {
        if (*flags) {
          buflane.emplace(index, bl, ch, p);
        }
        ++flags;
      }
    }
  }
  buflane.write_end();

  for (std::thread& t : threads) t.join();
}

void Interpolate::interpolationThread() {
  common::lane_read_buffer<Sample> buflane(&_lane, BUFFER_SIZE);
  Sample sample;
  while (buflane.read(sample)) {
    interpolateSample(sample.timestep, sample.baseline, sample.channel,
                      sample.pol);
  }
}

void Interpolate::interpolateSample(size_t timestep, size_t baseline,
                                    size_t channel, size_t pol) {
  const casacore::IPosition shp = _buffers.front().getData().shape();
  const size_t nPol = shp[0], nChan = shp[1],
               timestepBegin = (timestep > _windowSize / 2)
                                   ? (timestep - _windowSize / 2)
                                   : 0,
               timestepEnd =
                   std::min(timestep + _windowSize / 2 + 1, _buffers.size()),
               channelBegin = (channel > _windowSize / 2)
                                  ? (channel - _windowSize / 2)
                                  : 0,
               channelEnd = std::min(channel + _windowSize / 2 + 1, nChan);

  std::complex<float> valueSum = 0.0;
  float windowSum = 0.0;

  for (size_t t = timestepBegin; t != timestepEnd; ++t) {
    casacore::Complex* data = _buffers[t].getData().data() +
                              (baseline * nChan + channelBegin) * nPol + pol;
    const bool* flags = _buffers[t].getFlags().data() +
                        (baseline * nChan + channelBegin) * nPol + pol;
    const float* row =
        &_kernelLookup[_windowSize * (t + int(_windowSize / 2) - timestep)];
    for (size_t ch = channelBegin; ch != channelEnd; ++ch) {
      if (!*flags) {
        int x = ch + int(_windowSize / 2) - channel;
        float w = row[x];
        valueSum += *data * w;
        windowSum += w;
      }

      data += nPol;
      flags += nPol;
    }
  }
  // This write is multithreaded, but is allowed because this value is never
  // read from in the loops above (because flagged values are skipped).
  casacore::Complex& value =
      _buffers[timestep]
          .getData()
          .data()[(baseline * nChan + channel) * nPol + pol];
  if (windowSum != 0.0)
    value = valueSum / windowSum;
  else
    value = casacore::Complex(std::numeric_limits<float>::quiet_NaN(),
                              std::numeric_limits<float>::quiet_NaN());
}

}  // namespace steps
}  // namespace dp3
