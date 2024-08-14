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
#include <aocommon/logger.h>

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

namespace {

// ApplyBeam is templated on the type of the data, could be complex<double>
// or complex<float>
template <typename T>
void ApplyBeamToData(const DPInfo& info, const size_t n_stations, T* data0,
                     float* weight0, std::vector<aocommon::MC2x2>& beam_values,
                     bool doUpdateWeights) {
  /*
    Applies the beam to each baseline and each frequency of the
    model patch
  */
  const size_t n_channels = info.chanFreqs().size();
  const size_t n_baselines = info.nbaselines();

  // Apply beam for channel ch on all baselines
  // For mode=ARRAY_FACTOR, too much work is done here
  // because we know that r and l are diagonal
  for (size_t bl = 0; bl < n_baselines; ++bl) {
    // If the beam is the same for all stations (i.e. when n_stations = 1),
    // all baselines will have the same beam values
    size_t index_left = (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
    size_t index_right =
        (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);
    for (size_t ch = 0; ch < n_channels; ++ch) {
      T* data = data0 + bl * 4 * n_channels + ch * 4;
      const aocommon::MC2x2F mat(data);

      const aocommon::MC2x2F left(beam_values[index_left + ch]);
      const aocommon::MC2x2F right(beam_values[index_right + ch]);
      const aocommon::MC2x2F result = left * mat.MultiplyHerm(right);
      result.AssignTo(data);
      if (doUpdateWeights) {
        dp3::steps::ApplyCal::ApplyWeights(
            left, right, weight0 + bl * 4 * n_channels + ch * 4);
      }
    }
  }
}

template <typename T>
void ApplyBeamToDataAndAdd(const DPInfo& info, const size_t n_stations,
                           T* data0, T* model_data, float* weight0,
                           std::vector<aocommon::MC2x2>& beam_values,
                           bool doUpdateWeights) {
  /*
    Applies the beam to each baseline and each frequency of the
    model patch and sum the contribution to the model data
  */
  const size_t n_channels = info.chanFreqs().size();
  const size_t n_baselines = info.nbaselines();

  // Apply beam for channel ch on all baselines
  // For mode=ARRAY_FACTOR, too much work is done here
  // because we know that r and l are diagonal
  for (size_t bl = 0; bl < n_baselines; ++bl) {
    // If the beam is the same for all stations (i.e. when n_stations = 1),
    // all baselines will have the same beam values
    const size_t index_left =
        (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
    const size_t index_right =
        (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);
    for (size_t ch = 0; ch < n_channels; ++ch) {
      T* data_ptr = data0 + bl * 4 * n_channels + ch * 4;
      T* model_data_ptr = model_data + bl * 4 * n_channels + ch * 4;
      const aocommon::MC2x2F data(data_ptr);
      const aocommon::MC2x2F model_data(model_data_ptr);

      const aocommon::MC2x2F left(beam_values[index_left + ch]);
      const aocommon::MC2x2F right(beam_values[index_right + ch]);
      const aocommon::MC2x2F result =
          model_data + left * data.MultiplyHerm(right);
      result.AssignTo(model_data_ptr);
      if (doUpdateWeights) {
        dp3::steps::ApplyCal::ApplyWeights(
            left, right, weight0 + bl * 4 * n_channels + ch * 4);
      }
    }
  }
}

size_t ComputeBeam(const DPInfo& info, double time,
                   const everybeam::vector3r_t& srcdir,
                   const everybeam::telescope::Telescope* telescope,
                   std::vector<aocommon::MC2x2>& beam_values, bool invert,
                   everybeam::CorrectionMode mode, std::mutex* mutex,
                   const std::vector<size_t>& skip_station_indices) {
  /*
    Compute the beam values for each station in a specific direction
    and store them into beam_values

    For convenience it returns the number of stations the beam
    was computed for.
  */
  const size_t n_channels = info.chanFreqs().size();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      dp3::base::SelectStationIndices(*telescope, info.antennaNames());
  const size_t n_stations = station_indices.size();

  // Apply the beam values of both stations to the ApplyBeamed data.
  for (size_t ch = 0; ch < n_channels; ++ch) {
    switch (mode) {
      case everybeam::CorrectionMode::kFull:
      case everybeam::CorrectionMode::kElement:
        // Fill beam_values for channel ch
        for (size_t st = 0; st < n_stations; ++st) {
          if (std::find(skip_station_indices.begin(),
                        skip_station_indices.end(),
                        st) != skip_station_indices.end()) {
            beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
          } else {
            beam_values[n_channels * st + ch] = point_response->Response(
                mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);

            if (invert) {
              // Terminate if the matrix is not invertible.
              [[maybe_unused]] bool status =
                  beam_values[n_channels * st + ch].Invert();
              assert(status);
            }
          }
        }
        break;
      case everybeam::CorrectionMode::kArrayFactor: {
        aocommon::MC2x2 af_tmp;
        for (size_t st = 0; st < n_stations; ++st) {
          if (std::find(skip_station_indices.begin(),
                        skip_station_indices.end(),
                        st) != skip_station_indices.end()) {
            beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
          } else {
            af_tmp = point_response->Response(
                mode, station_indices[st], info.chanFreqs()[ch], srcdir, mutex);

            if (invert) {
              af_tmp[0] = 1. / af_tmp[0];
              af_tmp[3] = 1. / af_tmp[3];
            }
            beam_values[n_channels * st + ch] = af_tmp;
          }
        }
        break;
      }
      case everybeam::CorrectionMode::kNone:  // this should not happen
        for (size_t st = 0; st < n_stations; ++st) {
          beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
        }
        break;
    }
  }
  return n_stations;
}
}  // namespace

namespace dp3 {
namespace steps {

ApplyBeam::ApplyBeam(const common::ParameterSet& parset, const string& prefix,
                     bool substep)
    : itsName(prefix),
      itsUpdateWeights(parset.getBool(prefix + "updateweights", false)),
      itsDirectionStr(parset.getStringVector(prefix + "direction",
                                             std::vector<std::string>())),
      itsUseChannelFreq(parset.getBool(prefix + "usechannelfreq", true)),
      itsSkipStationNames(parset.getStringVector(prefix + "skipstations",
                                                 std::vector<std::string>())),
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

  const size_t n_stations = info().nantenna();
  const size_t n_channels = info().nchan();

  const size_t nThreads = aocommon::ThreadPool::GetInstance().NThreads();
  itsBeamValues.resize(nThreads);

  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  itsMeasConverters.resize(nThreads);
  itsMeasFrames.resize(nThreads);
  telescopes_.resize(nThreads);

  for (size_t thread = 0; thread < nThreads; ++thread) {
    itsBeamValues[thread].resize(n_stations * n_channels);
    itsMeasFrames[thread].set(info().arrayPosCopy());
    itsMeasFrames[thread].set(
        MEpoch(MVEpoch(info().startTime() / 86400), MEpoch::UTC));
    itsMeasConverters[thread].set(
        MDirection::J2000,
        MDirection::Ref(MDirection::ITRF, itsMeasFrames[thread]));
    telescopes_[thread] = base::GetTelescope(
        info().msName(), itsElementResponseModel, itsUseChannelFreq);

    telescopes_[thread]->SetTime(info().startTime());
  }

  if (!itsSkipStationNames.empty()) {
    // Needs loop over itsSkipStationNames because SelectStationIndices
    // assumes some order. By giving it a length-one vector the order is as it
    // assumes (because there is only one way to order a vector of length one.
    for (std::string& skipStationName : itsSkipStationNames) {
      std::vector<size_t> station_indices = base::SelectStationIndices(
          *(telescopes_[0]), std::vector<std::string>{skipStationName});
      itsSkipStationIndices.emplace_back(std::move(station_indices[0]));
    }
  }
}

void ApplyBeam::show(std::ostream& os) const {
  os << "ApplyBeam " << itsName << '\n'
     << "  mode:              " << everybeam::ToString(itsMode) << '\n'
     << "  use channelfreq:   " << std::boolalpha << itsUseChannelFreq << '\n'
     << "  direction:         " << itsDirectionStr << '\n'
     << "  invert:            " << std::boolalpha << itsInvert << '\n'
     << "  update weights:    " << std::boolalpha << itsUpdateWeights << '\n';
  if (!itsSkipStationNames.empty()) {
    os << "  not applying for:  [";
    for (const std::string& name : itsSkipStationNames) {
      os << name << (name == itsSkipStationNames.back() ? "]\n" : ",");
    }
  }
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

  telescopes_[thread]->SetTime(time);

  if (undoInputBeam) {
    // A beam was previously applied to this MS, and a different direction
    // was asked this time. 'Undo' applying the input beam.
    // TODO itsElementResponseModel should be read from the measurement set
    // instead of assumed to be the same from the target beam.
    applyBeam(info(), time, data, weight, srcdir, telescopes_[thread].get(),
              itsBeamValues[thread], false, itsModeAtStart, itsUpdateWeights,
              nullptr, itsSkipStationIndices);
    srcdir = dir2Itrf(itsDirection, itsMeasConverters[thread]);
  }

  applyBeam(info(), time, data, weight, srcdir, telescopes_[thread].get(),
            itsBeamValues[thread], itsInvert, itsMode, itsUpdateWeights,
            nullptr, itsSkipStationIndices);

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

// applyBeam is templated on the type of the data, could be complex<double>
// or complex<float>
template <typename T>
void ApplyBeam::applyBeam(const DPInfo& info, double time, T* data0,
                          float* weight0, const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          std::vector<aocommon::MC2x2>& beam_values,
                          bool invert, everybeam::CorrectionMode mode,
                          bool doUpdateWeights, std::mutex* mutex,
                          const std::vector<size_t>& skip_station_indices) {
  const size_t n_stations =
      ComputeBeam(info, time, srcdir, telescope, beam_values, invert, mode,
                  mutex, skip_station_indices);
  ApplyBeamToData(info, n_stations, data0, weight0, beam_values,
                  doUpdateWeights);
}

// applyBeam is templated on the type of the data, could be complex<double>
// or complex<float>
template <typename T>
void ApplyBeam::ApplyBeamAndAddToModel(
    const DPInfo& info, double time, T* data0, T* model_data, float* weight0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<aocommon::MC2x2>& beam_values, bool invert,
    everybeam::CorrectionMode mode, bool doUpdateWeights, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices) {
  const size_t n_stations =
      ComputeBeam(info, time, srcdir, telescope, beam_values, invert, mode,
                  mutex, skip_station_indices);

  ApplyBeamToDataAndAdd(info, n_stations, data0, model_data, weight0,
                        beam_values, doUpdateWeights);
}

template void ApplyBeam::ApplyBeamAndAddToModel(
    const DPInfo& info, double time, std::complex<double>* data0,
    std::complex<double>* model_data, float* weight0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<aocommon::MC2x2>& beam_values, bool invert,
    everybeam::CorrectionMode mode, bool doUpdateWeights, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices);

template void ApplyBeam::ApplyBeamAndAddToModel(
    const DPInfo& info, double time, std::complex<float>* data0,
    std::complex<float>* model_data, float* weight0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<aocommon::MC2x2>& beam_values, bool invert,
    everybeam::CorrectionMode mode, bool doUpdateWeights, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices);

template void ApplyBeam::applyBeam(
    const DPInfo& info, double time, std::complex<double>* data0,
    float* weight0, const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<aocommon::MC2x2>& beam_values, bool invert,
    everybeam::CorrectionMode mode, bool doUpdateWeights, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices);

void ApplyBeam::applyBeam(const DPInfo& info, double time,
                          std::complex<double>* data0, float* weight0,
                          const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          std::vector<aocommon::MC2x2>& beam_values,
                          const std::pair<size_t, size_t>& baseline_range,
                          const std::pair<size_t, size_t>& station_range,
                          aocommon::Barrier& barrier, bool invert,
                          everybeam::CorrectionMode mode,
                          bool do_update_weights, std::mutex* mutex,
                          const std::vector<size_t>& skip_station_indices) {
  // Get the beam values for each station.
  const size_t n_channels = info.chanFreqs().size();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());
  const size_t n_stations = station_indices.size();

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
            if (std::find(skip_station_indices.begin(),
                          skip_station_indices.end(),
                          st) != skip_station_indices.end()) {
              beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
            } else {
              beam_values[n_channels * st + ch] =
                  point_response->Response(mode, station_indices[st],
                                           info.chanFreqs()[ch], srcdir, mutex);
              if (invert) {
                // Terminate if the matrix is not invertible.
                [[maybe_unused]] bool status =
                    beam_values[n_channels * st + ch].Invert();
                assert(status);
              }
            }
          }
          break;
        case everybeam::CorrectionMode::kArrayFactor: {
          aocommon::MC2x2 af_tmp;
          // Fill beam_values for channel ch
          // only for stations used by this thread
          for (size_t st = station_range.first; st < station_range.second;
               ++st) {
            if (std::find(skip_station_indices.begin(),
                          skip_station_indices.end(),
                          st) != skip_station_indices.end()) {
              beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
            } else {
              af_tmp =
                  point_response->Response(mode, station_indices[st],
                                           info.chanFreqs()[ch], srcdir, mutex);

              if (invert) {
                af_tmp[0] = 1. / af_tmp[0];
                af_tmp[3] = 1. / af_tmp[3];
              }
              beam_values[n_channels * st + ch] = af_tmp;
            }
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
    // Apply beam for channel ch on the baselines handled by this thread
    // For mode=ARRAY_FACTOR, too much work is done here because we know
    // that r and l are diagonal
    for (size_t bl = baseline_range.first; bl < baseline_range.second; ++bl) {
      std::complex<double>* data = data0 + bl * 4 * n_channels + ch * 4;
      const aocommon::MC2x2F mat(data);
      size_t index_left =
          (n_stations == 1 ? ch : n_channels * info.getAnt1()[bl] + ch);
      size_t index_right =
          (n_stations == 1 ? ch : n_channels * info.getAnt2()[bl] + ch);
      const aocommon::MC2x2F left(beam_values[index_left]);
      const aocommon::MC2x2F right(beam_values[index_right]);

      const aocommon::MC2x2F result = left * mat.MultiplyHerm(right);
      result.AssignTo(data);
      if (do_update_weights) {
        ApplyCal::ApplyWeights(left, right,
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
    std::vector<everybeam::complex_t>& beam_values, bool invert,
    everybeam::CorrectionMode mode, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices) {
  assert(mode == everybeam::CorrectionMode::kArrayFactor);
  // Get the beam values for each station.
  const size_t n_channels = info.chanFreqs().size();
  const size_t n_baselines = info.nbaselines();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());
  const size_t n_stations = station_indices.size();

  // Apply the beam values of both stations to the ApplyBeamed data.
  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Fill beam_values for channel ch
    for (size_t st = 0; st < n_stations; ++st) {
      if (std::find(skip_station_indices.begin(), skip_station_indices.end(),
                    st) != skip_station_indices.end()) {
        beam_values[n_channels * st + ch] = 1.0;
      } else {
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

  // Apply beam for channel ch on all baselines
  for (size_t bl = 0; bl < n_baselines; ++bl) {
    size_t index_left = (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
    size_t index_right =
        (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);

    for (size_t ch = 0; ch < n_channels; ++ch) {
      T* data = data0 + bl * n_channels + ch;
      everybeam::complex_t* left = &(beam_values[index_left]);
      everybeam::complex_t* right = &(beam_values[index_right]);
      data[0] = left[ch] * std::complex<double>(data[0]) * conj(right[ch]);

      // TODO: update weights?
    }
  }
}

template void ApplyBeam::applyBeamStokesIArrayFactor(
    const DPInfo& info, double time, std::complex<double>* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<everybeam::complex_t>& beam_values, bool invert,
    everybeam::CorrectionMode mode, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices);

void ApplyBeam::applyBeamStokesIArrayFactor(
    const DPInfo& info, double time, std::complex<double>* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::vector<everybeam::complex_t>& beam_values,
    const std::pair<size_t, size_t>& baseline_range,
    const std::pair<size_t, size_t>& station_range, aocommon::Barrier& barrier,
    bool invert, everybeam::CorrectionMode mode, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices) {
  assert(mode == everybeam::CorrectionMode::kArrayFactor);
  // Get the beam values for each station.
  const size_t n_channels = info.chanFreqs().size();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());
  const size_t n_stations = station_indices.size();

  // Apply the beam values of both stations to the ApplyBeamed data.
  if (station_range.first < n_stations) {
    for (size_t st = station_range.first; st < station_range.second; ++st) {
      if (std::find(skip_station_indices.begin(), skip_station_indices.end(),
                    st) != skip_station_indices.end()) {
        for (size_t ch = 0; ch < n_channels; ++ch) {
          beam_values[n_channels * st + ch] = 1.0;
        }
      } else {
        for (size_t ch = 0; ch < n_channels; ++ch) {
          // Fill beam_values for channel ch
          // only for stations used by this thread
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
  }
  barrier.wait();
  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Apply beam for channel ch on the baselines handeled by this thread
    for (size_t bl = baseline_range.first; bl < baseline_range.second; ++bl) {
      std::complex<double>* data = data0 + bl * n_channels + ch;
      size_t index_left =
          (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
      size_t index_right =
          (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);
      everybeam::complex_t* left = &(beam_values[index_left]);
      everybeam::complex_t* right = &(beam_values[index_right]);
      data[0] = left[ch] * std::complex<double>(data[0]) * conj(right[ch]);

      // TODO: update weights?
    }
  }
  barrier.wait();
}

}  // namespace steps
}  // namespace dp3
