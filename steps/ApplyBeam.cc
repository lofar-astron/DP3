// ApplyBeam.cc: DPPP step class to ApplyBeam visibilities
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema

#include "../common/ParameterSet.h"
#include "../common/Timer.h"
#include "../common/StreamUtil.h"
#include "../common/StringTools.h"

#include "ApplyBeam.h"
#include "ApplyCal.h"
// for matrix inversion
#include <dp3/base/DPInfo.h>
#include "../base/FlagCounter.h"
#include "../base/Telescope.h"

#include <EveryBeam/telescope/phasedarray.h>
#include <EveryBeam/pointresponse/pointresponse.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>

#include <aocommon/matrix2x2diag.h>
#include <aocommon/threadpool.h>

#include <cassert>
#include <ostream>
#include <stddef.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using casacore::MDirection;
using casacore::MEpoch;
using casacore::MVAngle;
using casacore::MVEpoch;
using casacore::Quantity;

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::operator<<;

namespace dp3 {
namespace steps {

ApplyBeam::ApplyBeam(const common::ParameterSet& parset, const string& prefix,
                     bool substep)
    : itsName(prefix),
      itsUpdateWeights(parset.getBool(prefix + "updateweights", false)),
      itsDirectionStr(parset.getStringVector(prefix + "direction",
                                             std::vector<std::string>())),
      itsUseChannelFreq(parset.getBool(prefix + "usechannelfreq", true)),
      itsMode(everybeam::ParseCorrectionMode(
          parset.getString(prefix + "beammode", "default"))),
      itsModeAtStart(everybeam::CorrectionMode::kNone),
      itsDebugLevel(parset.getInt(prefix + "debuglevel", 0)) {
  // only read 'invert' parset key if it is a separate step
  // if applybeam is called from gaincal/predict, the invert key should always
  // be false
  if (substep) {
    itsInvert = false;
  } else {
    itsInvert = parset.getBool(prefix + "invert", true);
  }

  string element_model = boost::to_lower_copy(
      parset.getString(prefix + "elementmodel", "hamaker"));
  if (element_model == "hamaker") {
    itsElementResponseModel = everybeam::ElementResponseModel::kHamaker;
  } else if (element_model == "lobes") {
    itsElementResponseModel = everybeam::ElementResponseModel::kLOBES;
  } else if (element_model == "oskar") {
    itsElementResponseModel =
        everybeam::ElementResponseModel::kOSKARSphericalWave;
  } else if (element_model == "oskardipole") {
    itsElementResponseModel = everybeam::ElementResponseModel::kOSKARDipole;
  } else {
    throw std::runtime_error(
        "Elementmodel should be HAMAKER, LOBES, OSKAR or OSKARDIPOLE");
  }
}

void ApplyBeam::updateInfo(const DPInfo& infoIn) {
  info() = infoIn;

  // Parse direction parset value
  if (itsDirectionStr.empty())
    itsDirection = info().phaseCenter();
  else {
    if (itsDirectionStr.size() != 2)
      throw std::runtime_error(
          "2 values must be given in direction option of ApplyBeam");
    casacore::MDirection phaseCenter;
    Quantity q0, q1;
    if (!MVAngle::read(q0, itsDirectionStr[0]))
      throw std::runtime_error(
          itsDirectionStr[0] +
          " is an invalid RA or longitude in ApplyBeam direction");
    if (!MVAngle::read(q1, itsDirectionStr[1]))
      throw std::runtime_error(
          itsDirectionStr[1] +
          " is an invalid DEC or latitude in ApplyBeam direction");
    MDirection::Types type = MDirection::J2000;
    itsDirection = MDirection(q0, q1, type);
  }

  if (itsInvert) {
    itsModeAtStart =
        static_cast<everybeam::CorrectionMode>(info().beamCorrectionMode());
    itsDirectionAtStart = info().beamCorrectionDir();
    info().setBeamCorrectionMode(static_cast<int>(itsMode));
    info().setBeamCorrectionDir(itsDirection);
  } else {
    const auto mode =
        static_cast<everybeam::CorrectionMode>(info().beamCorrectionMode());
    if (mode == everybeam::CorrectionMode::kNone)
      throw std::runtime_error(
          "In applying the beam (with invert=false): the metadata of this "
          "observation indicate that the beam has not yet been applied");
    if (mode != itsMode)
      throw std::runtime_error(
          std::string("applybeam step with invert=false has incorrect mode: "
                      "input has ") +
          everybeam::ToString(mode) + ", requested to correct for " +
          everybeam::ToString(itsMode));
    const double ra1 = info().beamCorrectionDir().getValue().getValue()[0];
    const double dec1 = info().beamCorrectionDir().getValue().getValue()[1];
    const double ra2 = itsDirection.getValue().getValue()[0];
    const double dec2 = itsDirection.getValue().getValue()[1];
    const double raDist = std::fabs(ra1 - ra2);
    const double decDist = std::fabs(dec1 - dec2);
    if (raDist > 1e-9 || decDist > 1e-9) {
      std::ostringstream str;
      str << "applybeam step with invert=false has incorrect direction: input "
             "is for "
          << info().beamCorrectionDir() << ", output is for " << itsDirection;
      throw std::runtime_error(str.str());
    }
    info().setBeamCorrectionMode(
        static_cast<int>(everybeam::CorrectionMode::kNone));
  }

  const size_t nSt = info().nantenna();
  const size_t nCh = info().nchan();

  const size_t nThreads = aocommon::ThreadPool::GetInstance().NThreads();
  itsBeamValues.resize(nThreads);

  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  itsMeasConverters.resize(nThreads);
  itsMeasFrames.resize(nThreads);
  telescopes_.resize(nThreads);

  for (size_t thread = 0; thread < nThreads; ++thread) {
    itsBeamValues[thread].resize(nSt * nCh);
    itsMeasFrames[thread].set(info().arrayPosCopy());
    itsMeasFrames[thread].set(
        MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
    itsMeasConverters[thread].set(
        MDirection::J2000,
        MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
    telescopes_[thread] = base::GetTelescope(
        info().msName(), itsElementResponseModel, itsUseChannelFreq);
  }
}

void ApplyBeam::show(std::ostream& os) const {
  os << "ApplyBeam " << itsName << '\n'
     << "  mode:              " << everybeam::ToString(itsMode) << '\n'
     << "  use channelfreq:   " << std::boolalpha << itsUseChannelFreq << '\n'
     << "  direction:         " << itsDirectionStr << '\n'
     << "  invert:            " << std::boolalpha << itsInvert << '\n'
     << "  update weights:    " << std::boolalpha << itsUpdateWeights << '\n';
  if (itsInvert) {
    if (itsModeAtStart != everybeam::CorrectionMode::kNone)
      os << "  input data has already a beam correction applied: will be "
            "undone.\n";
    else
      os << "  input data has no beam correction applied.\n";
  }
}

void ApplyBeam::showTimings(std::ostream& os, double duration) const {
  os << "  ";
  base::FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " ApplyBeam " << itsName << '\n';
}

bool ApplyBeam::processMultithreaded(std::unique_ptr<base::DPBuffer> buffer,
                                     size_t thread) {
  itsTimer.start();

  std::complex<float>* data = buffer->GetData().data();

  float* weight = buffer->GetWeights().data();

  const double time = buffer->GetTime();

  // Set up directions for beam evaluation
  everybeam::vector3r_t srcdir;

  /**
   * I'm not sure this is correct the way it is. These loops
   * seem to initialize variables that are never used in a
   * multi-threaded way, and if they were used from multiple
   * threads, it would imply process() is called multiple times,
   * and hence this initialization is already subject to a race
   * condition... ???
   * André, 2018-10-07
   */
  bool undoInputBeam =
      itsInvert && itsModeAtStart != everybeam::CorrectionMode::kNone;
  const size_t nThreads = aocommon::ThreadPool::GetInstance().NThreads();
  for (size_t threadIter = 0; threadIter < nThreads; ++threadIter) {
    itsMeasFrames[threadIter].resetEpoch(
        MEpoch(MVEpoch(time / 86400), MEpoch::UTC));
    // Do a conversion on all threads, because converters are not
    // thread safe and apparently need to be used at least once
    if (undoInputBeam)
      srcdir = dir2Itrf(itsDirectionAtStart, itsMeasConverters[threadIter]);
    else
      srcdir = dir2Itrf(itsDirection, itsMeasConverters[threadIter]);
  }

  if (undoInputBeam) {
    // A beam was previously applied to this MS, and a different direction
    // was asked this time. 'Undo' applying the input beam.
    // TODO itsElementResponseModel should be read from the measurement set
    // instead of assumed to be the same from the target beam.
    applyBeam(info(), time, data, weight, srcdir, telescopes_[thread].get(),
              itsBeamValues[thread], false, itsModeAtStart, itsUpdateWeights);
    srcdir = dir2Itrf(itsDirection, itsMeasConverters[thread]);
  }

  applyBeam(info(), time, data, weight, srcdir, telescopes_[thread].get(),
            itsBeamValues[thread], itsInvert, itsMode, itsUpdateWeights);

  itsTimer.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

everybeam::vector3r_t ApplyBeam::dir2Itrf(const MDirection& dir,
                                          MDirection::Convert& measConverter) {
  const MDirection& itrfDir = measConverter(dir);
  const casacore::Vector<double>& itrf = itrfDir.getValue().getValue();
  return {itrf[0], itrf[1], itrf[2]};
}

void ApplyBeam::finish() {
  // Let the next steps finish.
  getNextStep()->finish();
}

// applyBeam is templated on the type of the data, could be complex<double> or
// complex<float>
template <typename T>
void ApplyBeam::applyBeam(const DPInfo& info, double time, T* data0,
                          float* weight0, const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          std::vector<aocommon::MC2x2>& beamValues, bool invert,
                          everybeam::CorrectionMode mode, bool doUpdateWeights,
                          std::mutex* mutex) {
  // Get the beam values for each station.
  const size_t nCh = info.chanFreqs().size();
  const size_t nSt = beamValues.size() / nCh;
  const size_t nBl = info.nbaselines();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());

  // Apply the beam values of both stations to the ApplyBeamed data.
  for (size_t ch = 0; ch < nCh; ++ch) {
    switch (mode) {
      case everybeam::CorrectionMode::kFull:
      case everybeam::CorrectionMode::kElement:
        // Fill beamValues for channel ch
        for (size_t st = 0; st < nSt; ++st) {
          beamValues[nCh * st + ch] = point_response->Response(
              mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);
          if (invert) {
            // Terminate if the matrix is not invertible.
            [[maybe_unused]] bool status = beamValues[nCh * st + ch].Invert();
            assert(status);
          }
        }
        break;
      case everybeam::CorrectionMode::kArrayFactor: {
        aocommon::MC2x2 af_tmp;
        // Fill beamValues for channel ch
        for (size_t st = 0; st < nSt; ++st) {
          af_tmp = point_response->Response(
              mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);

          if (invert) {
            af_tmp[0] = 1. / af_tmp[0];
            af_tmp[3] = 1. / af_tmp[3];
          }
          beamValues[nCh * st + ch] = af_tmp;
        }
        break;
      }
      case everybeam::CorrectionMode::kNone:  // this should not happen
        for (size_t st = 0; st < nSt; ++st) {
          beamValues[nCh * st + ch] = aocommon::MC2x2::Unity();
        }
        break;
    }

    // Apply beam for channel ch on all baselines
    // For mode=ARRAY_FACTOR, too much work is done here because we know
    // that r and l are diagonal
    for (size_t bl = 0; bl < nBl; ++bl) {
      T* data = data0 + bl * 4 * nCh + ch * 4;
      const aocommon::MC2x2F mat(data);
      //
      const aocommon::MC2x2F left(
          beamValues[nCh * info.getAnt1()[bl] + ch].Data());
      const aocommon::MC2x2F right(
          beamValues[nCh * info.getAnt2()[bl] + ch].Data());
      const aocommon::MC2x2F result = left.Multiply(mat).MultiplyHerm(right);
      result.AssignTo(data);
      if (doUpdateWeights) {
        ApplyCal::ApplyWeights(left.Data(), right.Data(),
                               weight0 + bl * 4 * nCh + ch * 4);
      }
    }
  }
}

template void ApplyBeam::applyBeam(
    const DPInfo& info, double time, std::complex<double>* data0,
    float* weight0, const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<aocommon::MC2x2>& beamValues, bool invert,
    everybeam::CorrectionMode mode, bool doUpdateWeights, std::mutex* mutex);

void ApplyBeam::applyBeam(const DPInfo& info, double time,
                          std::complex<double>* data0, float* weight0,
                          const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          std::vector<aocommon::MC2x2>& beam_values,
                          const std::pair<size_t, size_t>& baseline_range,
                          const std::pair<size_t, size_t>& station_range,
                          aocommon::Barrier& barrier, bool invert,
                          everybeam::CorrectionMode mode,
                          bool do_update_weights, std::mutex* mutex) {
  // Get the beam values for each station.
  const size_t n_channels = info.chanFreqs().size();
  const size_t n_stations = beam_values.size() / n_channels;

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());

  // Apply the beam values of both stations to the ApplyBeamed data.

  if (station_range.first < n_stations) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      switch (mode) {
        case everybeam::CorrectionMode::kFull:
        case everybeam::CorrectionMode::kElement:
          // Fill beam_values for channel ch,
          // only for stations used by this thread
          for (size_t st = station_range.first; st < station_range.second;
               ++st) {
            beam_values[n_channels * st + ch] = point_response->Response(
                mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);
            if (invert) {
              // Terminate if the matrix is not invertible.
              [[maybe_unused]] bool status =
                  beam_values[n_channels * st + ch].Invert();
              assert(status);
            }
          }
          break;
        case everybeam::CorrectionMode::kArrayFactor: {
          aocommon::MC2x2 af_tmp;
          // Fill beam_values for channel ch
          // only for stations used by this thread
          for (size_t st = station_range.first; st < station_range.second;
               ++st) {
            af_tmp = point_response->Response(
                mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);

            if (invert) {
              af_tmp[0] = 1. / af_tmp[0];
              af_tmp[3] = 1. / af_tmp[3];
            }
            beam_values[n_channels * st + ch] = af_tmp;
          }
          break;
        }
        case everybeam::CorrectionMode::kNone:  // this should not happen
          for (size_t st = station_range.first; st < station_range.second;
               ++st) {
            beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
          }
          break;
      }
    }
  }

  barrier.wait();
  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Apply beam for channel ch on the baselines handeled by this thread
    // For mode=ARRAY_FACTOR, too much work is done here because we know
    // that r and l are diagonal
    for (size_t bl = baseline_range.first; bl < baseline_range.second; ++bl) {
      std::complex<double>* data = data0 + bl * 4 * n_channels + ch * 4;
      const aocommon::MC2x2F mat(data);
      //
      const aocommon::MC2x2F left(
          beam_values[n_channels * info.getAnt1()[bl] + ch].Data());
      const aocommon::MC2x2F right(
          beam_values[n_channels * info.getAnt2()[bl] + ch].Data());
      const aocommon::MC2x2F result = left.Multiply(mat).MultiplyHerm(right);
      result.AssignTo(data);
      if (do_update_weights) {
        ApplyCal::ApplyWeights(left.Data(), right.Data(),
                               weight0 + bl * 4 * n_channels + ch * 4);
      }
    }
  }

  barrier.wait();
}

template <typename T>
void ApplyBeam::applyBeamStokesIArrayFactor(
    const DPInfo& info, double time, T* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<everybeam::complex_t>& beamValues, bool invert,
    everybeam::CorrectionMode mode, std::mutex* mutex) {
  assert(mode == everybeam::CorrectionMode::kArrayFactor);
  // Get the beam values for each station.
  const size_t nCh = info.chanFreqs().size();
  const size_t nSt = beamValues.size() / nCh;
  const size_t nBl = info.nbaselines();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());

  // Apply the beam values of both stations to the ApplyBeamed data.
  for (size_t ch = 0; ch < nCh; ++ch) {
    // Fill beamValues for channel ch
    for (size_t st = 0; st < nSt; ++st) {
      beamValues[nCh * st + ch] = point_response->Response(
          everybeam::BeamMode::kArrayFactor, station_indices[st],
          info.chanFreqs()[ch], srcdir, mutex)[0];
      if (invert) {
        beamValues[nCh * st + ch] = 1. / beamValues[nCh * st + ch];
      }
    }

    // Apply beam for channel ch on all baselines
    for (size_t bl = 0; bl < nBl; ++bl) {
      T* data = data0 + bl * nCh + ch;
      everybeam::complex_t* left = &(beamValues[nCh * info.getAnt1()[bl]]);
      everybeam::complex_t* right = &(beamValues[nCh * info.getAnt2()[bl]]);
      data[0] = left[ch] * std::complex<double>(data[0]) * conj(right[ch]);

      // TODO: update weights?
    }
  }
}

template void ApplyBeam::applyBeamStokesIArrayFactor(
    const DPInfo& info, double time, std::complex<double>* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<everybeam::complex_t>& beamValues, bool invert,
    everybeam::CorrectionMode mode, std::mutex* mutex);

void ApplyBeam::applyBeamStokesIArrayFactor(
    const DPInfo& info, double time, std::complex<double>* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<everybeam::complex_t>& beam_values,
    const std::pair<size_t, size_t>& baseline_range,
    const std::pair<size_t, size_t>& station_range, aocommon::Barrier& barrier,
    bool invert, everybeam::CorrectionMode mode, std::mutex* mutex) {
  assert(mode == everybeam::CorrectionMode::kArrayFactor);
  // Get the beam values for each station.
  const size_t n_channels = info.chanFreqs().size();
  const size_t n_stations = beam_values.size() / n_channels;

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());

  // Apply the beam values of both stations to the ApplyBeamed data.
  if (station_range.first < n_stations) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Fill beam_values for channel ch
      // only for stations used by this thread
      for (size_t st = station_range.first; st < station_range.second; ++st) {
        beam_values[n_channels * st + ch] = point_response->Response(
            everybeam::BeamMode::kArrayFactor, station_indices[st],
            info.chanFreqs()[ch], srcdir, mutex)[0];
        if (invert) {
          beam_values[n_channels * st + ch] =
              1. / beam_values[n_channels * st + ch];
        }
      }
    }
  }
  barrier.wait();
  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Apply beam for channel ch on the baselines handeled by this thread
    for (size_t bl = baseline_range.first; bl < baseline_range.second; ++bl) {
      std::complex<double>* data = data0 + bl * n_channels + ch;
      everybeam::complex_t* left =
          &(beam_values[n_channels * info.getAnt1()[bl]]);
      everybeam::complex_t* right =
          &(beam_values[n_channels * info.getAnt2()[bl]]);
      data[0] = left[ch] * std::complex<double>(data[0]) * conj(right[ch]);

      // TODO: update weights?
    }
  }
  barrier.wait();
}

}  // namespace steps
}  // namespace dp3
