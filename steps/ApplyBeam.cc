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
#include "base/DPInfo.h"
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
                     float* weight0, aocommon::MC2x2* beam_values,
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

}  // namespace

namespace dp3 {
namespace steps {

size_t ComputeBeam(const base::DPInfo& info, double time,
                   const everybeam::vector3r_t& srcdir,
                   const everybeam::telescope::Telescope* telescope,
                   aocommon::MC2x2* beam_values, bool invert,
                   everybeam::BeamMode mode, std::mutex* mutex,
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
  switch (mode) {
    case everybeam::BeamMode::kFull:
    case everybeam::BeamMode::kElement:
      // Fill beam_values for channel ch
      for (size_t st = 0; st < n_stations; ++st) {
        if (std::find(skip_station_indices.begin(), skip_station_indices.end(),
                      st) != skip_station_indices.end()) {
          for (size_t ch = 0; ch < n_channels; ++ch)
            beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
        } else {
          point_response->Response(&beam_values[n_channels * st], mode,
                                   station_indices[st], info.chanFreqs(),
                                   srcdir, mutex);
          if (invert) {
            for (size_t ch = 0; ch < n_channels; ++ch) {
              // Terminate if the matrix is not invertible.
              [[maybe_unused]] const bool status =
                  beam_values[n_channels * st + ch].Invert();
              assert(status);
            }
          }
        }
      }
      break;
    case everybeam::BeamMode::kArrayFactor: {
      for (size_t st = 0; st < n_stations; ++st) {
        if (std::find(skip_station_indices.begin(), skip_station_indices.end(),
                      st) != skip_station_indices.end()) {
          for (size_t ch = 0; ch < n_channels; ++ch)
            beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
        } else {
          point_response->Response(&beam_values[n_channels * st], mode,
                                   station_indices[st], info.chanFreqs(),
                                   srcdir, mutex);

          for (size_t ch = 0; ch < n_channels; ++ch) {
            if (invert) {
              const aocommon::MC2x2 af_tmp = beam_values[n_channels * st + ch];
              beam_values[n_channels * st + ch] = aocommon::MC2x2(
                  1.0 / af_tmp.Get(0), 0.0, 0.0, 1.0 / af_tmp.Get(3));
            }
          }
        }
      }
      break;
    }
    case everybeam::BeamMode::kNone:  // this should not happen
      for (size_t st = 0; st < n_stations; ++st) {
        for (size_t ch = 0; ch < n_channels; ++ch) {
          beam_values[n_channels * st + ch] = aocommon::MC2x2::Unity();
        }
      }
      break;
  }
  return n_stations;
}

void ApplyBeamToDataAndAdd(
    const DPInfo& info, size_t n_stations,
    const aocommon::xt::UTensor<std::complex<double>, 3>& data,
    aocommon::xt::UTensor<std::complex<double>, 3>& model_data,
    const aocommon::MC2x2* beam_values) {
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
      std::complex<double>* model_data_pointer = &model_data(bl, ch, 0);
      // Using (model_)data2x2 avoids aliases with the (model_)data arguments.
      const aocommon::MC2x2F data2x2(&data(bl, ch, 0));
      const aocommon::MC2x2F model_data2x2(model_data_pointer);

      const aocommon::MC2x2F left(beam_values[index_left + ch]);
      const aocommon::MC2x2F right(beam_values[index_right + ch]);
      const aocommon::MC2x2F result =
          model_data2x2 + left * data2x2.MultiplyHerm(right);
      result.AssignTo(model_data_pointer);
    }
  }
}

ApplyBeam::ApplyBeam(const common::ParameterSet& parset,
                     const std::string& prefix, bool substep)
    : itsName(prefix),
      itsUpdateWeights(parset.getBool(prefix + "updateweights", false)),
      itsDirectionStr(parset.getStringVector(prefix + "direction",
                                             std::vector<std::string>())),
      itsUseChannelFreq(parset.getBool(prefix + "usechannelfreq", true)),
      itsSkipStationNames(parset.getStringVector(prefix + "skipstations",
                                                 std::vector<std::string>())),
      itsMode(everybeam::ParseBeamMode(
          parset.getString(prefix + "beammode", "default"))),
      coefficients_path_(parset.getString(prefix + "coefficients_path", "")),
      itsDebugLevel(parset.getInt(prefix + "debuglevel", 0)),
      use_model_data_(parset.getBool(prefix + "usemodeldata", false)) {
  // only read 'invert' parset key if it is a separate step
  // if applybeam is called from gaincal/predict, the invert key should always
  // be false
  if (substep) {
    itsInvert = false;
  } else {
    itsInvert = parset.getBool(prefix + "invert", true);
  }

  std::string element_model =
      parset.getString(prefix + "elementmodel", "default");
  itsElementResponseModel =
      everybeam::ElementResponseModelFromString(element_model);
}

void ApplyBeam::updateInfo(const DPInfo& infoIn) {
  Step::updateInfo(infoIn);

  // If ApplyBeam was requested to apply the beam to model data,
  // check whether there is model data to apply the beam to.
  if (use_model_data_ && getInfoOut().GetDirections().empty()) {
    throw std::runtime_error(
        "ApplyBeam's option 'usemodeldata' is set to true, \n"
        "but the beam can not be applied to model data, \n"
        "because no model data with direction meta data is found.");
  }

  // Parse direction parset value
  if (itsDirectionStr.empty())
    itsDirection = getInfoOut().phaseCenter();
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
        static_cast<everybeam::BeamMode>(getInfoOut().beamCorrectionMode());
    itsDirectionAtStart = getInfoOut().beamCorrectionDir();
    GetWritableInfoOut().setBeamCorrectionMode(static_cast<int>(itsMode));
    GetWritableInfoOut().setBeamCorrectionDir(itsDirection);
  } else if (!use_model_data_) {
    const auto mode =
        static_cast<everybeam::BeamMode>(getInfoOut().beamCorrectionMode());
    if (mode == everybeam::BeamMode::kNone)
      throw std::runtime_error(
          "In applying the beam (with invert=false): the metadata of this "
          "observation indicate that the beam has not yet been applied");
    if (mode != itsMode)
      throw std::runtime_error(
          std::string("applybeam step with invert=false has incorrect mode: "
                      "input has ") +
          everybeam::ToString(mode) + ", requested to correct for " +
          everybeam::ToString(itsMode));
    const double ra1 =
        getInfoOut().beamCorrectionDir().getValue().getValue()[0];
    const double dec1 =
        getInfoOut().beamCorrectionDir().getValue().getValue()[1];
    const double ra2 = itsDirection.getValue().getValue()[0];
    const double dec2 = itsDirection.getValue().getValue()[1];
    const double raDist = std::fabs(ra1 - ra2);
    const double decDist = std::fabs(dec1 - dec2);
    if (raDist > 1e-9 || decDist > 1e-9) {
      std::ostringstream str;
      str << "applybeam step with invert=false has incorrect direction: input "
             "is for "
          << getInfoOut().beamCorrectionDir() << ", output is for "
          << itsDirection;
      throw std::runtime_error(str.str());
    }
    GetWritableInfoOut().setBeamCorrectionMode(
        static_cast<int>(everybeam::BeamMode::kNone));
  }

  const size_t n_stations = getInfoOut().nantenna();
  const size_t n_channels = getInfoOut().nchan();

  // Create the Measure ITRF conversion info given the array position.
  // The time and direction are filled in later.
  beam_values_.resize(n_stations * n_channels);
  measure_frame_.set(getInfoOut().arrayPosCopy());
  measure_frame_.set(
      MEpoch(MVEpoch(getInfoOut().startTime() / 86400), MEpoch::UTC));
  measure_converter_.set(MDirection::J2000,
                         MDirection::Ref(MDirection::ITRF, measure_frame_));
  telescope_ =
      base::GetTelescope(getInfoOut().msName(), itsElementResponseModel,
                         itsUseChannelFreq, coefficients_path_);
  telescope_->SetTime(getInfoOut().startTime());

  if (!itsSkipStationNames.empty()) {
    // Needs loop over itsSkipStationNames because SelectStationIndices
    // assumes some order. By giving it a length-one vector the order is as it
    // assumes (because there is only one way to order a vector of length one.
    for (std::string& skipStationName : itsSkipStationNames) {
      std::vector<size_t> station_indices = base::SelectStationIndices(
          *telescope_, std::vector<std::string>{skipStationName});
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
    if (itsModeAtStart != everybeam::BeamMode::kNone)
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

bool ApplyBeam::ProcessModelData(std::unique_ptr<base::DPBuffer> buffer) {
  itsTimer.start();

  const std::map<std::string, dp3::base::Direction>& directions =
      getInfoOut().GetDirections();
  const double time = buffer->GetTime();
  telescope_->SetTime(time);
  measure_frame_.resetEpoch(MEpoch(MVEpoch(time / 86400), MEpoch::UTC));

  for (const auto& [direction_name, direction] : directions) {
    std::complex<float>* data = buffer->GetData(direction_name).data();
    casacore::MDirection direction_j2000(
        casacore::Quantity(direction.ra, "rad"),
        casacore::Quantity(direction.dec, "rad"), MDirection::J2000);
    everybeam::vector3r_t direction_itrf =
        dir2Itrf(direction_j2000, measure_converter_);

    const size_t n_stations =
        ComputeBeam(getInfoOut(), time, direction_itrf, telescope_.get(),
                    beam_values_.data(), itsInvert, itsMode, nullptr,
                    itsSkipStationIndices);
    ApplyBeamToData(getInfoOut(), n_stations, data, nullptr,
                    beam_values_.data(), false);
  }

  itsTimer.stop();
  getNextStep()->process(std::move(buffer));
  return false;
}

bool ApplyBeam::ProcessData(std::unique_ptr<base::DPBuffer> buffer) {
  itsTimer.start();

  std::complex<float>* data = buffer->GetData().data();
  float* weight = buffer->GetWeights().data();
  const double time = buffer->GetTime();
  const bool undoInputBeam =
      itsInvert && itsModeAtStart != everybeam::BeamMode::kNone;
  measure_frame_.resetEpoch(MEpoch(MVEpoch(time / 86400), MEpoch::UTC));

  telescope_->SetTime(time);

  if (undoInputBeam) {
    // A beam was previously applied to this MS, and a different direction
    // was asked this time. 'Undo' applying the input beam.
    // TODO itsElementResponseModel should be read from the measurement set
    // instead of assumed to be the same from the target beam.
    const everybeam::vector3r_t srcdir =
        dir2Itrf(itsDirectionAtStart, measure_converter_);
    const size_t n_stations = ComputeBeam(
        getInfoOut(), time, srcdir, telescope_.get(), beam_values_.data(),
        false, itsModeAtStart, nullptr, itsSkipStationIndices);
    ApplyBeamToData(getInfoOut(), n_stations, data, weight, beam_values_.data(),
                    itsUpdateWeights);
  }

  const everybeam::vector3r_t srcdir =
      dir2Itrf(itsDirection, measure_converter_);
  const size_t n_stations = ComputeBeam(
      getInfoOut(), time, srcdir, telescope_.get(), beam_values_.data(),
      itsInvert, itsMode, nullptr, itsSkipStationIndices);
  ApplyBeamToData(getInfoOut(), n_stations, data, weight, beam_values_.data(),
                  itsUpdateWeights);

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

void ApplyBeam::ApplyBaselineBasedBeam(
    const DPInfo& info, double time, std::complex<double>* data0,
    float* weight0, const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    aocommon::MC2x2* beam_values,
    const std::pair<size_t, size_t>& baseline_range,
    const std::pair<size_t, size_t>& station_range, aocommon::Barrier& barrier,
    bool invert, everybeam::BeamMode mode, bool do_update_weights,
    std::mutex* mutex, const std::vector<size_t>& skip_station_indices) {
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
        case everybeam::BeamMode::kFull:
        case everybeam::BeamMode::kElement:
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
        case everybeam::BeamMode::kArrayFactor: {
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
                af_tmp = aocommon::MC2x2(1.0 / af_tmp.Get(0), 0.0, 0.0,
                                         1.0 / af_tmp.Get(3));
              }
              beam_values[n_channels * st + ch] = af_tmp;
            }
          }
          break;
        }
        case everybeam::BeamMode::kNone:  // this should not happen
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

size_t ComputeArrayFactor(const DPInfo& info, double time,
                          const everybeam::vector3r_t& srcdir,
                          const everybeam::telescope::Telescope* telescope,
                          std::complex<double>* beam_values, bool invert,
                          std::mutex* mutex,
                          const std::vector<size_t>& skip_station_indices) {
  const size_t n_channels = info.chanFreqs().size();

  std::unique_ptr<everybeam::pointresponse::PointResponse> point_response =
      telescope->GetPointResponse(time);

  const std::vector<size_t> station_indices =
      base::SelectStationIndices(*telescope, info.antennaNames());
  const size_t n_stations = station_indices.size();

  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Fill beam_values for channel ch
    for (size_t st = 0; st < n_stations; ++st) {
      std::complex<double>& value = beam_values[n_channels * st + ch];
      if (std::find(skip_station_indices.begin(), skip_station_indices.end(),
                    st) != skip_station_indices.end()) {
        value = 1.0;
      } else {
        value = point_response
                    ->Response(everybeam::BeamMode::kArrayFactor,
                               station_indices[st], info.chanFreqs()[ch],
                               srcdir, mutex)
                    .Get(0);
        if (invert) {
          value = 1.0 / value;
        }
      }
    }
  }
  return n_stations;
}

void ApplyArrayFactorAndAdd(
    const DPInfo& info, size_t n_stations,
    const aocommon::xt::UTensor<std::complex<double>, 3>& data,
    aocommon::xt::UTensor<std::complex<double>, 3>& model_data,
    const std::complex<double>* beam_values) {
  // Apply beam for channel ch on all baselines
  const size_t n_baselines = info.nbaselines();
  const size_t n_channels = info.chanFreqs().size();
  for (size_t bl = 0; bl < n_baselines; ++bl) {
    const size_t index_left =
        (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
    const size_t index_right =
        (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);

    const std::complex<double>* left = &beam_values[index_left];
    const std::complex<double>* right = &beam_values[index_right];

    // Using pointers ensures that indexing is as fast as possible in the loop
    // below, which is on a performance-critical path.
    const std::complex<double>* data_pointer = &data(bl, 0, 0);
    std::complex<double>* model_data_pointer = &model_data(bl, 0, 0);
    for (size_t ch = 0; ch < n_channels; ++ch) {
      model_data_pointer[ch] +=
          left[ch] * data_pointer[ch] * std::conj(right[ch]);

      // TODO: update weights?
    }
  }
}

void ApplyBeam::ApplyBaselineBasedArrayFactor(
    const DPInfo& info, double time, std::complex<double>* data0,
    const everybeam::vector3r_t& srcdir,
    const everybeam::telescope::Telescope* telescope,
    std::complex<double>* beam_values,
    const std::pair<size_t, size_t>& baseline_range,
    const std::pair<size_t, size_t>& station_range, aocommon::Barrier& barrier,
    bool invert, everybeam::BeamMode mode, std::mutex* mutex,
    const std::vector<size_t>& skip_station_indices) {
  assert(mode == everybeam::BeamMode::kArrayFactor);
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
          beam_values[n_channels * st + ch] =
              point_response
                  ->Response(everybeam::BeamMode::kArrayFactor,
                             station_indices[st], info.chanFreqs()[ch], srcdir,
                             mutex)
                  .Get(0);
          if (invert) {
            beam_values[n_channels * st + ch] =
                1.0 / beam_values[n_channels * st + ch];
          }
        }
      }
    }
  }
  barrier.wait();
  for (size_t ch = 0; ch < n_channels; ++ch) {
    // Apply beam for channel ch on the baselines handled by this thread
    for (size_t bl = baseline_range.first; bl < baseline_range.second; ++bl) {
      std::complex<double>* data = data0 + bl * n_channels + ch;
      size_t index_left =
          (n_stations == 1 ? 0 : n_channels * info.getAnt1()[bl]);
      size_t index_right =
          (n_stations == 1 ? 0 : n_channels * info.getAnt2()[bl]);
      std::complex<double>* left = &(beam_values[index_left]);
      std::complex<double>* right = &(beam_values[index_right]);
      data[0] = left[ch] * std::complex<double>(data[0]) * conj(right[ch]);

      // TODO: update weights?
    }
  }
  barrier.wait();
}

}  // namespace steps
}  // namespace dp3
